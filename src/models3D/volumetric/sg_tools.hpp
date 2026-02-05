//
//  sg_tools.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 5/6/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  tools for structured-grid based model

#ifndef sg_tools_hpp
#define sg_tools_hpp

#include "inparam.hpp"
#include "geodesy.hpp"
#include "StructuredGrid.hpp"
#include <iomanip> // required by gcc

namespace sg_tools {
    ////////////////////////// inparam //////////////////////////
    // horizontal coords
    inline
    void inparamHorizontal(const InparamYAML &gm, const std::string &rootc,
                           const std::string &modelName,
                           const std::string &className,
                           bool &sourceCentered, bool &xy, bool &ellipticity) {
        // type
        const std::string &type = gm.get<std::string>(rootc + ":horizontal");
        if (type == "DISTANCE_AZIMUTH") {
            sourceCentered = true;
            xy = false;
            ellipticity = false;
        } else if (type == "XY_CARTESIAN") {
            if (!geodesy::isCartesian()) {
                throw std::runtime_error("sg_tools::inparamHorizontal || "
                                         "XY_CARTESIAN is incompatible with a "
                                         "spherical mesh.");
            }
            sourceCentered = true;
            xy = true;
            ellipticity = false;
        } else if (type == "LATITUDE_LONGITUDE") {
            sourceCentered = false;
            xy = false;
            ellipticity = gm.get<bool>(rootc + ":ellipticity");
        } else {
            throw std::runtime_error
            ("sg_tools::inparamHorizontal || "
             "Parameter coordinates:horizontal must be DISTANCE_AZIMUTH, || "
             "XY_CARTESIAN or LATITUDE_LONGITUDE."
             " || Model name: " + modelName +
             " || Class name: " + className);
        }
    }
    
    // vertical coords
    inline void inparamVertical(const InparamYAML &gm, const std::string &rootc,
                                const std::string &modelName,
                                const std::string &className,
                                bool &useDepth, bool &depthSolid) {
        // type
        const std::string &type = gm.get<std::string>(rootc + ":vertical");
        if (type == "DEPTH") {
            useDepth = true;
            depthSolid = gm.get<bool>(rootc + ":depth_below_solid_surface");
        } else if (type == "RADIUS") {
            useDepth = false;
            depthSolid = false;
        } else {
            throw std::runtime_error
            ("sg_tools::inparamVertical || "
             "Parameter coordinates:vertical must be either RADIUS or DEPTH."
             " || Model name: " + modelName +
             " || Class name: " + className);
        }
    }
    
    // nc variables and ranks
    template <int NDIM>
    void inparamVarRank(const InparamYAML &gm, const std::string &rootc,
                        const std::string &modelName,
                        const std::string &className,
                        std::array<std::string, NDIM> &crdVarNames,
                        std::array<int, NDIM> &shuffleData) {
        const std::vector<std::string> &vars =
        gm.getVector<std::string>(rootc + ":nc_variables");
        const std::vector<int> &rank =
        gm.getVector<int>(rootc + ":data_rank");
        if (vars.size() != NDIM || rank.size() != NDIM) {
            throw std::runtime_error
            ("sg_tools::inparamVarRank || "
             "The size of coordinates:nc_variables and "
             "coordinates:data_rank must be " + bstring::toString(NDIM) + "."
             " || Model name: " + modelName +
             " || Class name: " + className);
        }
        std::copy(vars.begin(), vars.end(), crdVarNames.begin());
        std::copy(rank.begin(), rank.end(), shuffleData.begin());
    }
    
    // coord units
    inline void inparamUnits(const InparamYAML &gm, const std::string &rootc,
                             bool xy, double &lengthUnit, double &angleUnit) {
        lengthUnit = gm.getWithOptions<double>(rootc + ":length_unit", {
            {"m", 1.}, {"km", 1e3}});
        if (!xy) {
            angleUnit = gm.getWithLimits<double>(rootc + ":angle_unit", {
                {"degree", numerical::dDegree}, {"radian", 1.}});
        }
    }
    
    
    ////////////////////////// construct //////////////////////////
    // coord units
    template <class Grid>
    void constructUnits(Grid &grid, bool sourceCentered, bool xy,
                        bool vertical, double lengthUnit, double angleUnit) {
        // transform lambda expr
        auto vectorTimesFactor = [](std::vector<double> &vec, double factor) {
            std::transform(vec.begin(), vec.end(), vec.begin(),
                           [factor](double x) -> double {return factor * x;});};
        
        // horizontal
        if (sourceCentered) {
            if (geodesy::isCartesian()) {
                if (xy) {
                    // x, y
                    vectorTimesFactor(grid.getGridCoords()[0], lengthUnit);
                    vectorTimesFactor(grid.getGridCoords()[1], lengthUnit);
                } else {
                    // dist, azim
                    vectorTimesFactor(grid.getGridCoords()[0], lengthUnit);
                    vectorTimesFactor(grid.getGridCoords()[1], angleUnit);
                }
            } else {
                // dist, azim
                vectorTimesFactor(grid.getGridCoords()[0], angleUnit);
                vectorTimesFactor(grid.getGridCoords()[1], angleUnit);
            }
        } else {
            // lat, lon
            double angle = angleUnit / numerical::dDegree;
            vectorTimesFactor(grid.getGridCoords()[0], angle);
            vectorTimesFactor(grid.getGridCoords()[1], angle);
        }
        
        // vertical
        if (vertical) {
            vectorTimesFactor(grid.getGridCoords()[2], lengthUnit);
        }
    }
    
    // longitude range
    template <class Grid>
    bool constructLon360(const Grid &grid, const std::string &modelName) {
        const std::vector<double> &lon = grid.getGridCoords()[1];
        if (lon.front() < 0. - numerical::dEpsilon &&
            lon.back() > 180. + numerical::dEpsilon) {
            throw std::runtime_error
            ("sg_tools::constructLon360 || "
             "Longitude range must be either [-180, 180] or [0, 360]. || "
             "Model name = " + modelName);
        }
        return (lon.back() > 180. + numerical::dEpsilon);
    }
    
    
    ////////////////////////// verbose //////////////////////////
    // head
    inline std::string verboseHead(const std::string &modelName,
                                   const std::string &className,
                                   const std::string &fileName) {
        std::stringstream ss;
        ss << bstring::boxSubTitle(0, modelName + " ", '~');
        ss << bstring::boxEquals(2, 11, "class name", className);
        ss << bstring::boxEquals(2, 11, "NetCDF file", fileName);
        return ss.str();
    }
    
    // coords
    inline std::string
    verboseCoords(bool sourceCentered, bool xy, bool vertical, bool useDepth,
                  const std::vector<std::string> &crdVarNames,
                  const std::vector<double> &crdMin,
                  const std::vector<double> &crdMax) {
        std::stringstream ss;
        ss << bstring::boxSubTitle(2, "Coordinates");
        std::vector<std::string> crdNames;
        if (sourceCentered) {
            if (xy) {
                crdNames.push_back("x");
                crdNames.push_back("y");
            } else {
                crdNames.push_back("distance");
                crdNames.push_back("azimuth");
            }
        } else {
            crdNames.push_back("latitude");
            crdNames.push_back("longitude");
        }
        if (vertical) {
            crdNames.push_back(useDepth ? "depth" : "radius");
        }
        // width
        int widthCrd = std::max(vector_tools::maxLength(crdNames), 5);
        int widthVar = std::max(vector_tools::maxLength(crdVarNames), 6);
        // table
        ss << "      " << std::left;
        ss << std::setw(widthCrd) << "COORD" << " | ";
        ss << std::setw(widthVar) << "NC-VAR" << " | SCOPE" << "\n";
        for (int idim = 0; idim < crdNames.size(); idim++) {
            ss << "    * " << std::left;
            ss << std::setw(widthCrd) << crdNames[idim] << " | ";
            ss << std::setw(widthVar) << crdVarNames[idim] << " | ";
            ss << bstring::range(crdMin[idim], crdMax[idim]) << "\n";
        }
        return ss.str();
    }
}

#endif /* sg_tools_hpp */
