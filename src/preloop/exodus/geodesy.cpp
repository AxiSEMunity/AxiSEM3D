//
//  geodesy.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/16/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  geodesy of the model and transformations between these CSs:
//  * geographic: (lat, lon), north pole on top
//  * geocentric: (theta, phi), north pole on top
//  * source-centered: (distance, azimuth), source on top

#include "geodesy.hpp"
#include "inparam.hpp"
#include "ExodusMesh.hpp"
#include "mpi.hpp"

namespace geodesy {
    // internal data
    namespace internal {
        // Cartesain
        bool iCartesian = false;
        
        // radius
        double iOuterRadius = -1.;
        double iOuterSolidRadius = -1.;
        double iInnerRadius = -1.;
        
        // ellipticity
        double iOuterFlattening = -1.;
        std::vector<double> iEllipR;
        std::vector<double> iEllipF;
        
        // geographic
        // (lat, lon, r) on positive z-axis
        eigen::DRow3 iLatLonRadiusAxisZ = eigen::DRow3::Zero();
        // Q matrix from source-centered to geographic
        eigen::DMat33 iSrc2GeoQ = eigen::DMat33::Zero();
    }
    
    // setup
    void setup(ExodusMesh &exMesh) {
        using namespace numerical;
        // Cartesain
        internal::iCartesian = exMesh.isCartesian();
        
        // radius
        internal::iOuterRadius = exMesh.getMeshSurface();
        internal::iInnerRadius = exMesh.getMeshBottom();
        internal::iOuterSolidRadius = exMesh.getSolidSurface();
        
        ////////////////// ellipticity //////////////////
        // flattening on surface
        internal::iOuterFlattening = inparam::gInparamModel.
        getWithOptions<double>("geodesy:flattening_on_surface", {
            {"SPHERE", 0.},
            {"WGS84", 1. / 298.257223563},
            {"GRS80", 1. / 298.257222101},
            {"SPECFEM3D_GLOBE", 1. / 299.8}});
        // no ellipticity for Cartesian
        if (internal::iOuterFlattening < dEpsilon ||
            internal::iCartesian) {
            internal::iOuterFlattening = 0.;
        }
        
        // ellipticity curve (only on root in ExodusMesh)
        if (mpi::root()) {
            const eigen::DMatXX_RM &ellip = exMesh.getEllipticityCurve();
            internal::iEllipR.resize(ellip.cols());
            internal::iEllipF.resize(ellip.cols());
            for (int col = 0; col < ellip.cols(); col++) {
                internal::iEllipR[col] = ellip(0, col);
                internal::iEllipF[col] = ellip(1, col) *
                internal::iOuterFlattening;
            }
        }
        mpi::bcast(internal::iEllipR);
        mpi::bcast(internal::iEllipF);
        
        ////////////////// geographic //////////////////
        // (lat, lon, r) on positive z-axis
        const std::string &geographicZ = inparam::gInparamModel.
        get<std::string>("geodesy:lat_lon_north_pole_mesh");
        if (geographicZ == "SOURCE") {
            // read from inparam.source.yaml
            const std::string &sourceName = inparam::gInparamSource.
            get<std::string>("list_of_sources:{0}");
            const std::string &rootl = ("list_of_sources:[0]:" + sourceName +
                                        ":location");
            try {
                internal::iLatLonRadiusAxisZ(0) = inparam::gInparamSource.
                getWithBounds(rootl + ":latitude_longitude:[0]", -90., 90.);
                internal::iLatLonRadiusAxisZ(1) = inparam::gInparamSource.
                getWithBounds(rootl + ":latitude_longitude:[1]", -180., 180.);
            } catch (...) {
                throw std::runtime_error
                ("geodesy::setup || "
                 "Error determining the geographic location of the first || "
                 "source in list_of_sources in inparam.source.yaml; || "
                 "location:latitude_longitude must be presented.");
            }
            if (inparam::gInparamSource.contains(rootl + ":depth")) {
                double depth = inparam::gInparamSource.
                getWithBounds(rootl + ":depth", 0., getOuterRadius());
                if (inparam::gInparamSource.get<bool>
                    (rootl + ":depth_below_solid_surface")) {
                    internal::iLatLonRadiusAxisZ(2) =
                    getOuterSolidRadius() - depth;
                } else {
                    internal::iLatLonRadiusAxisZ(2) =
                    getOuterRadius() - depth;
                }
            } else  {
                internal::iLatLonRadiusAxisZ(2) = inparam::gInparamSource.
                getWithBounds(rootl + ":radius", 0., getOuterRadius());
            }
        } else {
            // read numbers
            const std::string &key = "geodesy:lat_lon_north_pole_mesh";
            internal::iLatLonRadiusAxisZ(0) =
            inparam::gInparamModel.getWithBounds(key + ":[0]", -90., 90.);
            internal::iLatLonRadiusAxisZ(1) =
            inparam::gInparamModel.getWithBounds(key + ":[1]", -180., 180.);
            // assume surface
            internal::iLatLonRadiusAxisZ(2) = getOuterRadius();
        }
        
        // Q matrix from source-centered to geographic
        eigen::DRow3 srctpr = llr2tpr(internal::iLatLonRadiusAxisZ, true);
        double theta = srctpr(0);
        double phi = srctpr(1);
        internal::iSrc2GeoQ(0, 0) = cos(theta) * cos(phi);
        internal::iSrc2GeoQ(0, 1) = -sin(phi);
        internal::iSrc2GeoQ(0, 2) = sin(theta) * cos(phi);
        internal::iSrc2GeoQ(1, 0) = cos(theta) * sin(phi);
        internal::iSrc2GeoQ(1, 1) = cos(phi);
        internal::iSrc2GeoQ(1, 2) = sin(theta) * sin(phi);
        internal::iSrc2GeoQ(2, 0) = -sin(theta);
        internal::iSrc2GeoQ(2, 1) = 0.;
        internal::iSrc2GeoQ(2, 2) = cos(theta);
        internal::iSrc2GeoQ.transposeInPlace();
    }
    
    // verbose
    std::string verbose(const eigen::DColX &discontinuities) {
        std::stringstream ss;
        ss << bstring::boxTitle("Geodesy");
        // Cartesain
        ss << bstring::boxEquals(0, 21, "model in Cartesian",
                                 internal::iCartesian);
        
        // radius
        ss << bstring::boxEquals(0, 21, "surface radius",
                                 internal::iOuterRadius);
        ss << bstring::boxEquals(0, 21, "solid surface radius",
                                 internal::iOuterSolidRadius);
        ss << bstring::boxEquals(0, 21, "bottom radius",
                                 internal::iInnerRadius);
        
        // (lat, lon, r) on positive z-axis
        ss << "a reference geographic location on the positive z-axis\n";
        ss << bstring::boxEquals(2, 19, "latitude",
                                 internal::iLatLonRadiusAxisZ[0]);
        ss << bstring::boxEquals(2, 19, "longitude",
                                 internal::iLatLonRadiusAxisZ[1]);
        ss << bstring::boxEquals(2, 19, "radius",
                                 internal::iLatLonRadiusAxisZ[2]);
        
        // ellipticity
        if (internal::iCartesian) {
            ss << "* No ellipticity correction for Cartesian models.\n";
        } else {
            if (internal::iOuterFlattening > numerical::dEpsilon) {
                ss << bstring::boxEquals(0, 21, "flattening on surface",
                                         internal::iOuterFlattening);
                ///////// f at discontinuities /////////
                ss << bstring::boxSubTitle(0, "flattening at discontinuities");
                // compute f
                const eigen::DColX &fDisc = computeFlattening(discontinuities);
                // width
                int width = (int)bstring::toString(internal::iOuterFlattening *
                                                   numerical::dPi / 3.).size();
                // verbose f from surface to center
                for (int idisc = (int)fDisc.size() - 1; idisc >=0; idisc--) {
                    ss << bstring::boxEquals(2, width + 6, "f = " +
                                             bstring::toString(fDisc[idisc]),
                                             discontinuities[idisc], "at");
                }
                ss << ("* Correcting for ellipticity in latitude-θ "
                       "conversions.\n");
            } else {
                ss << "* Perfect sphere, no ellipticity correction.\n";
            }
        }
        ss << bstring::boxBaseline() << "\n\n";
        return ss.str();
    }
}
