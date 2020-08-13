//
//  ExodusMesh.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 9/2/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  Exodus mesh created by salvus mesher

#include "ExodusMesh.hpp"

// tools
#include "timer.hpp"
#include "mpi.hpp"
#include "io.hpp"
#include "bstring.hpp"
#include "eigen_tools.hpp"
#include "Mapping.hpp"

// input
#include "NetCDF_Reader.hpp"

// nr
#include "NrField.hpp"


///////////////// construction /////////////////
// constructor
ExodusMesh::ExodusMesh(const std::string &meshFile):
mFileName(meshFile) {
    // memory info
    double memSup = 0.;
    double memAll = 0.;
    
    // opne file on root
    timer::gPreloopTimer.begin("Opening Exodus mesh file");
    timer::gPreloopTimer.message(io::popInputDir(mFileName));
    NetCDF_Reader reader;
    if (mpi::root()) {
        reader.open(io::popInputDir(mFileName));
    }
    timer::gPreloopTimer.ended("Opening Exodus mesh file");
    
    // read and broadcast
    timer::gPreloopTimer.begin("Reading and broadcasting mesh");
    
    // global variables and records
    timer::gPreloopTimer.begin("Global variables");
    readBcastGlobal(reader, memSup, memAll);
    timer::gPreloopTimer.ended("Global variables");
    
    // connectivity
    timer::gPreloopTimer.begin("Connectivity");
    readBcastConnectivity(reader, memSup, memAll);
    timer::gPreloopTimer.ended("Connectivity");
    
    // coords
    timer::gPreloopTimer.begin("Coordinates");
    readBcastCoordinates(reader, memSup, memAll);
    timer::gPreloopTimer.ended("Coordinates");
    
    // side sets
    timer::gPreloopTimer.begin("Side sets");
    readBcastSideSets(reader, memSup, memAll);
    timer::gPreloopTimer.ended("Side sets");
    
    // elemental variables
    timer::gPreloopTimer.begin("Element variables");
    readBcastElemental(reader, memSup, memAll);
    timer::gPreloopTimer.ended("Element variables");
    
    // radial variables
    timer::gPreloopTimer.begin("Radial variables");
    readBcastRadial(reader, memSup, memAll);
    timer::gPreloopTimer.ended("Radial variables");
    
    // ellipticity
    timer::gPreloopTimer.begin("Ellipticity");
    readEllipticity(reader, memSup, memAll);
    timer::gPreloopTimer.ended("Ellipticity");
    
    // report total memory
    std::stringstream ss;
    ss << "*** Non-scalable memory allocation on super processors: ";
    ss << memSup << " GB ***";
    timer::gPreloopTimer.message(ss.str());
    ss.str("");
    ss << "*** Non-scalable memory allocation on all processors: ";
    ss << memAll << " GB ***";
    timer::gPreloopTimer.message(ss.str());
    
    // finish
    timer::gPreloopTimer.ended("Reading and broadcasting mesh");
}

// verbose
std::string ExodusMesh::verbose() const {
    if (!mpi::root()) {
        // only root has complete data
        return "";
    }
    
    using namespace bstring;
    std::stringstream ss;
    ss << boxTitle("Exodus Mesh");
    
    // overview
    ss << boxSubTitle(0, "Overview");
    ss << boxEquals(2, 20, "Exodus file", mFileName);
    ss << boxEquals(2, 20, "mesh CS", mGlobalRecords.at("crdsys"));
    ss << boxEquals(2, 20, "model name", mGlobalRecords.at("model"));
    ss << boxEquals(2, 20, "isotropic", isIsotropic());
    ss << boxEquals(2, 20, "attenuation", hasAttenuation());
    ss << boxEquals(2, 20, "storage type", mElementNodesStorage ?
                    "element_nodes" : "elements");
    ss << boxEquals(2, 20, "# discontinuities", mDiscontinuities.size());
    
    // geometry
    ss << boxSubTitle(0, "Mesh geometry");
    ss << boxEquals(2, 20, "# elements", getNumQuads());
    ss << boxEquals(2, 20, "# nodes", getNumNodes());
    const eigen::DMatX2_RM &nodalCoords = mySuperOnly().mNodalCoords;
    if (isCartesian()) {
        const std::string &rangeS = range(nodalCoords.col(0).minCoeff(),
                                          nodalCoords.col(0).maxCoeff());
        const std::string &rangeZ = range(nodalCoords.col(1).minCoeff(),
                                          nodalCoords.col(1).maxCoeff());
        ss << boxEquals(2, 20, "range of s-axis", rangeS);
        ss << boxEquals(2, 20, "range of z-axis", rangeZ);
    } else {
        const eigen::DMatX2_RM &rt = geodesy::sz2rtheta(nodalCoords, true);
        const std::string &rangeR = range(rt.col(0).minCoeff(),
                                          rt.col(0).maxCoeff());
        const std::string &rangeT = range(rt.col(1).minCoeff(),
                                          rt.col(1).maxCoeff());
        ss << boxEquals(2, 20, "range of r-axis", rangeR);
        ss << boxEquals(2, 21, "range of θ-axis", rangeT);
    }
    
    // global
    ss << boxSubTitle(0, "Global variables");
    ss << boxEquals(2, 20, mGlobalVariables);
    ss << boxEquals(2, 20, mGlobalRecords);
    
    // element geometry
    ss << boxSubTitle(0, "Element geometry");
    const eigen::IColX &geometryType = mySuperOnly().mGeometryType;
    std::map<std::string, int> geomTypeMap = {
        {"spherical", geometryType.cwiseEqual(0).cast<int>().sum()},
        {"linear", geometryType.cwiseEqual(1).cast<int>().sum()},
        {"semi-spherical", geometryType.cwiseEqual(2).cast<int>().sum()}};
    ss << boxEquals(2, 20, geomTypeMap, ":");
    
    // solid-fluid
    ss << boxSubTitle(0, "Element medium");
    const eigen::IColX &isElementFluid = mySuperOnly().mIsElementFluid;
    std::map<std::string, int> solidFluidMap = {
        {"solid", getNumQuads() - isElementFluid.sum()},
        {"fluid", isElementFluid.sum()}};
    ss << boxEquals(2, 20, solidFluidMap, ":");
    
    // radial
    ss << boxSubTitle(0, "Material properties");
    for (auto it = mRadialVariables.begin();
         it != mRadialVariables.end(); ++it) {
        const std::string &prange = range(it->second.minCoeff(),
                                          it->second.maxCoeff());
        ss << boxEquals(2, 20, it->first, prange, "∈");
    }
    
    // side sets
    ss << boxSubTitle(0, "Side sets");
    for (auto it = mSideSets.begin(); it != mSideSets.end(); ++it) {
        ss << boxEquals(2, 20, it->first, toString(it->second.size())
                        + " sides", ":");
    }
    
    // miscellaneous (currently only ellipticity)
    ss << boxSubTitle(0, "Miscellaneous");
    ss << boxEquals(2, 20, "ellipticity curve",
                    toString(mEllipticityCurve.cols()) + " knots", ":");
    
    // mesh generation cmdline
    ss << boxBaseline('-');
    ss << "Mesh generation command line:\n>> ";
    std::vector<std::string> subcmd = split(mCmdMeshGen, "\n");
    for (int icmd = 0; icmd < subcmd.size() - 1; icmd++) {
        ss << subcmd[icmd] + "\n   ";
    }
    ss << subcmd.back() + "\n";
    
    // finish
    ss << boxBaseline() << "\n\n";
    return ss.str();
}

// global variables
void ExodusMesh::readBcastGlobal(const NetCDF_Reader &reader,
                                 double &memSup, double &memAll) {
    // read
    std::vector<std::string> glbVarNames, glbRecords;
    eigen::DRowX glbVarValues;
    if (mpi::root()) {
        reader.readString("name_glo_var", glbVarNames);
        reader.readString("info_records", glbRecords);
        reader.readMatrix("vals_glo_var", glbVarValues);
    }
    
    // broadcast
    mpi::bcast(glbVarNames);
    mpi::bcast(glbRecords);
    mpi::bcastEigen(glbVarValues);
    
    // cast variables to map
    for (int igv = 0; igv < glbVarNames.size(); igv++) {
        // "dt" could be misleading
        if (glbVarNames[igv] == "dt") {
            glbVarNames[igv] = "dt (nPol = 1)";
        }
        // make attenuation coefficients readable
        if (glbVarNames[igv].find("w_") == 0 ||
            glbVarNames[igv].find("y_") == 0 ||
            glbVarNames[igv].find("f_") == 0) {
            glbVarNames[igv] = "attenuation_" + glbVarNames[igv];
        }
        mGlobalVariables.insert({glbVarNames[igv], glbVarValues(igv)});
    }
    
    // cast records to map
    const std::vector<std::string> keys = {"crdsys", "model"};
    for (int igr = 0; igr < glbRecords.size(); igr++) {
        // split by "="
        std::vector<std::string> substrs = bstring::split(glbRecords[igr], "=");
        // only include required keys
        if (std::find(keys.begin(), keys.end(), substrs[0]) != keys.end()) {
            mGlobalRecords.insert({substrs[0], substrs[1]});
        }
        // mesh generation cmd
        if (substrs[0].find("cmdl") != std::string::npos) {
            mCmdMeshGen += substrs[1] + "\n";
        }
    }
}

// connectivity
void ExodusMesh::readBcastConnectivity(const NetCDF_Reader &reader,
                                       double &memSup, double &memAll) {
    mpi::enterSuper();
    if (mpi::super()) {
        // read
        if (mpi::root()) {
            reader.readMatrix("connect1", mySuperOnly().mConnectivity);
            // let node index start from 0
            mySuperOnly().mConnectivity.array() -= 1;
        }
        
        // broadcast
        mpi::bcastEigen(mySuperOnly().mConnectivity);
        // memory info
        timer::gPreloopTimer.message
        (eigen_tools::memoryInfo(mySuperOnly().mConnectivity,
                                 "connectivity (super-only)", memSup));
    }
    mpi::enterWorld();
}

// coordinates
void ExodusMesh::readBcastCoordinates(const NetCDF_Reader &reader,
                                      double &memSup, double &memAll) {
    mpi::enterSuper();
    if (mpi::super()) {
        // read
        if (mpi::root()) {
            // read
            eigen::DColX x, y;
            reader.readMatrix("coordx", x);
            reader.readMatrix("coordy", y);
            // cast to matrix
            mySuperOnly().mNodalCoords.resize(x.rows(), 2);
            mySuperOnly().mNodalCoords.col(0) = x;
            mySuperOnly().mNodalCoords.col(1) = y;
        }
        
        // broadcast
        mpi::bcastEigen(mySuperOnly().mNodalCoords);
        // memory info
        timer::gPreloopTimer.message
        (eigen_tools::memoryInfo(mySuperOnly().mNodalCoords,
                                 "coordinates (super-only)", memSup));
    }
    mpi::enterWorld();
}

// side sets
void ExodusMesh::readBcastSideSets(const NetCDF_Reader &reader,
                                   double &memSup, double &memAll) {
    // names
    std::vector<std::string> ssNames;
    if (mpi::root()) {
        reader.readString("ss_names", ssNames);
    }
    mpi::bcast(ssNames);
    
    // values
    for (int iss = 0; iss < ssNames.size(); iss++) {
        // read
        eigen::IColX elems, sides;
        if (mpi::root()) {
            const std::string &istr = bstring::toString(iss + 1);
            // elem
            reader.readMatrix("elem_ss" + istr, elems);
            // side
            reader.readMatrix("side_ss" + istr, sides);
            // let side index start from 0
            elems.array() -= 1;
            sides.array() -= 1;
        }
        
        // broadcast
        mpi::bcastEigen(elems);
        mpi::bcastEigen(sides);
        // memory info
        timer::gPreloopTimer.message
        (eigen_tools::memoryInfo(elems, ssNames[iss], memAll, 2));
        
        // cast to map
        std::map<int, int> mapES;
        for (int ie = 0; ie < elems.rows(); ie++) {
            mapES.insert({elems(ie), sides(ie)});
        }
        mSideSets.insert({ssNames[iss], mapES});
    }
    
    // side set keys
    mKeyLeftSS = isCartesian() ? "x0" : "t0";
    mKeyRightSS = isCartesian() ? "x1" : "t1";
    mKeyBottomSS = isCartesian() ? "y0" : "r0";
    mKeyTopSS = isCartesian() ? "y1" : "r1";
    
    // add minimum edge length to global variable
    double minEdge = std::numeric_limits<double>::max();
    if (mpi::root()) {
        for (int iquad = 0; iquad < getNumQuads(); iquad++) {
            // only consider axial elements
            if (getAxialSide(iquad) >= 0) {
                for (int pA = 0; pA < 4; pA++) {
                    int pB = (pA == 3) ? 0 : pA + 1;
                    const auto &szA = mySuperOnly().mNodalCoords.row
                    (mySuperOnly().mConnectivity(iquad, pA));
                    const auto &szB = mySuperOnly().mNodalCoords.row
                    (mySuperOnly().mConnectivity(iquad, pB));
                    minEdge = std::min(minEdge, (szA - szB).norm());
                }
            }
        }
    }
    mpi::bcast(minEdge);
    mGlobalVariables.insert({"min_edge_length", minEdge});
    mGlobalVariables.insert({"dist_tolerance", minEdge / 1000.});
}

// elemental variables
void ExodusMesh::readBcastElemental(const NetCDF_Reader &reader,
                                    double &memSup, double &memAll) {
    mpi::enterSuper();
    if (mpi::super()) {
        // read
        if (mpi::root()) {
            // get exodus keys
            std::vector<std::string> evNames;
            reader.readString("name_elem_var", evNames);
            
            // geometry type
            std::stringstream exKey;
            exKey << "vals_elem_var" <<
            (int)(std::find(evNames.begin(), evNames.end(), "element_type") -
                  evNames.begin() + 1) << "eb1";
            eigen::DRowX dbuffer;
            reader.readMatrix(exKey.str(), dbuffer);
            mySuperOnly().mGeometryType = dbuffer.array().round().cast<int>();
            
            // solid-fluid type
            exKey.str("");
            exKey << "vals_elem_var" <<
            (int)(std::find(evNames.begin(), evNames.end(), "fluid") -
                  evNames.begin() + 1) << "eb1";
            reader.readMatrix(exKey.str(), dbuffer);
            mySuperOnly().mIsElementFluid = dbuffer.array().round().cast<int>();
        }
        
        // broadcast
        mpi::bcastEigen(mySuperOnly().mGeometryType);
        mpi::bcastEigen(mySuperOnly().mIsElementFluid);
        // memory
        timer::gPreloopTimer.message
        (eigen_tools::memoryInfo(mySuperOnly().mGeometryType,
                                 "element_type (super-only)", memSup));
        timer::gPreloopTimer.message
        (eigen_tools::memoryInfo(mySuperOnly().mIsElementFluid,
                                 "fluid (super-only)", memSup));
    }
    mpi::enterWorld();
}

// radial variables
void ExodusMesh::readBcastRadial(const NetCDF_Reader &reader,
                                 double &memSup, double &memAll) {
    //////////////////////// name ////////////////////////
    // names of elemental variables
    std::vector<std::string> evNames;
    if (mpi::root()) {
        reader.readString("name_elem_var", evNames);
    }
    mpi::bcast(evNames);
    
    // name map for quick access
    std::map<std::string, std::string> evNameMap;
    for (int iev = 0; iev < evNames.size(); iev++) {
        std::stringstream exKey;
        exKey << "vals_elem_var" << iev + 1 << "eb1";
        evNameMap.insert({evNames[iev], exKey.str()});
    }
    
    // storage type
    mElementNodesStorage = (evNameMap.find("RHO_0") != evNameMap.end());
    
    // names of radial variables
    std::vector<std::string> rvNames;
    // isIsotropic() is not callable here
    bool isotropic = (evNameMap.find("ETA_0") == evNameMap.end() &&
                      evNameMap.find("ETA") == evNameMap.end());
    if (isotropic) {
        rvNames = std::vector<std::string>({"RHO", "VP", "VS"});
    } else {
        rvNames = std::vector<std::string>({"RHO", "VPV", "VPH",
            "VSV", "VSH", "ETA"});
    }
    if (hasAttenuation()) {
        rvNames.push_back("QMU");
        rvNames.push_back("QKAPPA");
    }
    
    
    //////////////////////// values ////////////////////////
    // dist tolerance of the mesh
    double distTol = mGlobalVariables.at("dist_tolerance");
    
    // values
    eigen::DMatXX rvValues;
    if (mpi::root()) {
        //////////////////////// depth profile ////////////////////////
        // coords of axial points
        std::vector<double> &coords = mRadialCoords;
        // quad and node tags of axial points
        std::vector<std::pair<int, int>> quad_nodes;
        for (int iquad = 0; iquad < getNumQuads(); iquad++) {
            // only consider axial elements
            int axisSide = getAxialSide(iquad);
            if (axisSide >= 0) {
                // coords
                int otherNode = (axisSide == 3) ? 0 : axisSide + 1;
                double z1 = mySuperOnly().mNodalCoords
                (mySuperOnly().mConnectivity(iquad, axisSide), 1);
                double z2 = mySuperOnly().mNodalCoords
                (mySuperOnly().mConnectivity(iquad, otherNode), 1);
                // negative values not needed
                if (std::min(z1, z2) < 0.) {
                    continue;
                }
                // move both ends inward by tolerance
                z1 += (z2 - z1) / std::abs(z2 - z1) * distTol;
                z2 -= (z2 - z1) / std::abs(z2 - z1) * distTol;
                // add to profile, keep sorted
                auto ir1 = coords.insert
                (std::upper_bound(coords.begin(), coords.end(), z1), z1);
                quad_nodes.insert(quad_nodes.begin() + (ir1 - coords.begin()),
                                  {iquad, axisSide});
                auto ir2 = coords.insert
                (std::upper_bound(coords.begin(), coords.end(), z2), z2);
                quad_nodes.insert(quad_nodes.begin() + (ir2 - coords.begin()),
                                  {iquad, otherNode});
            }
        }
        
        
        //////////////////////// depth values ////////////////////////
        rvValues = eigen::DMatXX::Zero(coords.size(), rvNames.size());
        for (int iname = 0; iname < rvNames.size(); iname++) {
            const std::string &rvName = rvNames[iname];
            if (mElementNodesStorage) {
                // storage = element_nodes
                std::array<eigen::DRowX, 4> buf;
                reader.readMatrix(evNameMap.at(rvName + "_0"), buf[0]);
                reader.readMatrix(evNameMap.at(rvName + "_1"), buf[1]);
                reader.readMatrix(evNameMap.at(rvName + "_2"), buf[2]);
                reader.readMatrix(evNameMap.at(rvName + "_3"), buf[3]);
                for (int idep = 0; idep < coords.size(); idep++) {
                    rvValues(idep, iname) = buf[quad_nodes[idep].second]
                    (quad_nodes[idep].first);
                }
            } else {
                // storage = elements
                eigen::DColX buf;
                reader.readMatrix(evNameMap.at(rvName), buf);
                for (int idep = 0; idep < coords.size(); idep++) {
                    rvValues(idep, iname) = buf(quad_nodes[idep].first);
                }
            }
        }
        
        
        //////////////////////// discontinuities ////////////////////////
        // read
        eigen::DColX allDiscs;
        reader.readMatrix("discontinuities", allDiscs);
        
        // remove those out of model range
        std::vector<double> myDiscs;
        for (int idisc = 0; idisc < allDiscs.size(); idisc++) {
            if (allDiscs(idisc) > coords.front() - 2. * distTol &&
                allDiscs(idisc) < coords.back() + 2. * distTol) {
                myDiscs.push_back(allDiscs(idisc));
            }
        }
        mDiscontinuities = Eigen::Map<eigen::DColX>
        (myDiscs.data(), myDiscs.size());
        
        // verify discontinuities
        int disc_found = 0;
        for (int idep = 0; idep < coords.size() - 1; idep++) {
            if (coords[idep + 1] - coords[idep] > distTol * 4) {
                // this gap spans an element
                // do nothing
                continue;
            }
            // check gap between discontinuities
            bool gapIsDisc = false;
            for (int idisc = 0; idisc < mDiscontinuities.size(); idisc++) {
                if (mDiscontinuities[idisc] < coords[idep + 1] &&
                    mDiscontinuities[idisc] > coords[idep]) {
                    // this gap spans a discontinuity
                    // keep it as is
                    gapIsDisc = true;
                    disc_found++;
                    break;
                }
            }
            // this gap is fake
            if (!gapIsDisc) {
                // average the upper and the lower values
                // such difference results from "elements" storage
                rvValues.row(idep) = (rvValues.row(idep) +
                                      rvValues.row(idep + 1)) * .5;
                rvValues.row(idep + 1) = rvValues.row(idep);
            }
        }
        // the top and bottom won't be found
        if (disc_found != mDiscontinuities.size() - 2) {
            throw std::runtime_error("ExodusMesh::readBcastRadial || "
                                     "Error verifying mesh discontinuities.");
        }
    }
    
    // broadcast
    mpi::bcastEigen(mDiscontinuities);
    mpi::bcast(mRadialCoords);
    mpi::bcastEigen(rvValues);
    // memory infor
    timer::gPreloopTimer.message
    (eigen_tools::memoryInfo(mDiscontinuities, "discontinuities", memAll));
    timer::gPreloopTimer.message
    (eigen_tools::memoryInfo
     (Eigen::Map<eigen::DColX>(mRadialCoords.data(), mRadialCoords.size()),
      "anchoring radii", memAll));
    
    // cast to map
    for (int irv = 0; irv < rvNames.size(); irv++) {
        mRadialVariables.insert({rvNames[irv], rvValues.col(irv)});
        // memory infor
        timer::gPreloopTimer.message
        (eigen_tools::memoryInfo(rvValues.col(irv), rvNames[irv], memAll));
    }
    
    // rotate nodes of axial elements such that side 3 lies on the axis
    // NOTE: doing this after forming the radial variables; otherwise the
    //       elemental variables must also be rotated
    for (auto ita = mSideSets.at(mKeyLeftSS).begin();
         ita != mSideSets.at(mKeyLeftSS).end(); ita++) {
        int iquad = ita->first;
        int axialSide = ita->second;
        // connectivity
        mpi::enterSuper();
        if (mpi::super()) {
            // cannot use const reference here
            eigen::IMatX4_RM connect = mySuperOnly().mConnectivity.row(iquad);
            for (int ivt = 0; ivt < 4; ivt++) {
                mySuperOnly().mConnectivity(iquad, ivt) =
                connect(Mapping::cycle4(ivt + axialSide - 3));
            }
        }
        mpi::enterWorld();
        // side sets
        for (auto it = mSideSets.begin(); it != mSideSets.end(); ++it) {
            if (getSide(it->first, iquad) != -1) {
                it->second.at(iquad) =
                Mapping::cycle4(it->second.at(iquad) - axialSide + 3);
            }
        }
    }
}

// ellipticity
void ExodusMesh::readEllipticity(const NetCDF_Reader &reader,
                                 double &memSup, double &memAll) {
    if (mpi::root()) {
        if (!isCartesian()) {
            // read
            reader.readMatrix("ellipticity", mEllipticityCurve);
        } else {
            // constant 1 between [0, 1]
            mEllipticityCurve = eigen::DMatXX_RM::Ones(2, 2);
            mEllipticityCurve(0, 0) = 0.;
        }
        // to absolute
        mEllipticityCurve.row(0) *= getMeshSurface();
        // message (broadcast later in geodesy)
        timer::gPreloopTimer.message
        (eigen_tools::memoryInfo(mEllipticityCurve, "ellipticity", memAll));
    }
}


///////////////// nr field /////////////////
// form Nr at nodes
void ExodusMesh::formNrAtNodes(const NrField &nrField,
                               bool boundByInplane, bool useLuckyNumbers) {
    // partitioning task on ranks
    int nNodePP = getNumNodes() / mpi::nproc();
    int nNodeOnMe = nNodePP;
    if (mpi::rank() == mpi::nproc() - 1) {
        nNodeOnMe += getNumNodes() % mpi::nproc();
    }
    const eigen::IColX &idR = eigen::IColX::LinSpaced
    (nNodeOnMe, mpi::rank() * nNodePP, mpi::rank() * nNodePP + nNodeOnMe - 1);
    
    // compute Nr on ranks
    timer::gPreloopTimer.begin("Computing Nr on ranks");
    const eigen::DMatX2_RM &crdMe = mySuperOnly().mNodalCoords(idR, Eigen::all);
    eigen::IColX nrMe = nrField.getNrAtPoints(crdMe);
    // check Nr >= 1
    Eigen::Index loc = -1;
    int minNr = nrMe.minCoeff(&loc);
    if (minNr < 1) {
        std::stringstream ss;
        ss << "ExodusMesh::formNrAtNodes || Nr must be positive. || ";
        ss << "Nr = " << minNr << " || Location (s, z) = ";
        ss << bstring::range(crdMe.row(loc)(0), crdMe.row(loc)(1), '(', ')');
        throw std::runtime_error(ss.str());
    }
    timer::gPreloopTimer.ended("Computing Nr on ranks");
    
    // bound by inplane resolution
    if (boundByInplane) {
        timer::gPreloopTimer.begin("Bounding Nr by inplane resolution");
        // node-quad reference list
        timer::gPreloopTimer.begin("Generating node-quad reference list");
        std::vector<std::vector<int>> refQuads;
        refQuads.resize(getNumNodes());
        for (int iquad = 0; iquad < getNumQuads(); iquad++) {
            refQuads[mySuperOnly().mConnectivity(iquad, 0)].push_back(iquad);
            refQuads[mySuperOnly().mConnectivity(iquad, 1)].push_back(iquad);
            refQuads[mySuperOnly().mConnectivity(iquad, 2)].push_back(iquad);
            refQuads[mySuperOnly().mConnectivity(iquad, 3)].push_back(iquad);
        }
        timer::gPreloopTimer.ended("Generating node-quad reference list");
        
        // average gll spacing
        timer::gPreloopTimer.begin("Computing average GLL spacing");
        for (int inode = 0; inode < nNodeOnMe; inode++) {
            double aveGLLSpacing = 0.;
            int inodeG = inode + mpi::rank() * nNodePP;
            for (int iquad: refQuads[inodeG]) {
                const eigen::DRow2 &sz0 = mySuperOnly().mNodalCoords
                .row(mySuperOnly().mConnectivity(iquad, 0));
                const eigen::DRow2 &sz1 = mySuperOnly().mNodalCoords
                .row(mySuperOnly().mConnectivity(iquad, 1));
                const eigen::DRow2 &sz2 = mySuperOnly().mNodalCoords
                .row(mySuperOnly().mConnectivity(iquad, 2));
                const eigen::DRow2 &sz3 = mySuperOnly().mNodalCoords
                .row(mySuperOnly().mConnectivity(iquad, 3));
                aveGLLSpacing += (sz0 - sz1).norm();
                aveGLLSpacing += (sz1 - sz2).norm();
                aveGLLSpacing += (sz2 - sz3).norm();
                aveGLLSpacing += (sz3 - sz0).norm();
            }
            aveGLLSpacing /= refQuads[inodeG].size() * 4 * spectral::nPol;
            // upper bound of Nr
            double circ = 2. * numerical::dPi * crdMe(inode, 0);
            int upperNr = std::max((int)round(circ / aveGLLSpacing), 1);
            // for axial nodes, 0 is not enough
            if (crdMe(inode, 0) < getGlobalVariable("dist_tolerance")) {
                // Q: Use 3 or 5 here? This is very tricky.
                // A: An axial moment tensor requires nu=2 (nr=5), but source
                //    is added at element level and handled by non-axial nodes;
                //    therefore, 3 is enough for an axial GLL point, as used in
                //    the old code. However, here we implement nr on the nodes,
                //    and thus 5 must be used on the axis so that the
                //    "interpolated" nr on a non-axial GLL point can reach 5.
                upperNr = 5;
            }
            nrMe(inode) = std::min(nrMe(inode), upperNr);
        }
        timer::gPreloopTimer.ended("Computing average GLL spacing");
        timer::gPreloopTimer.ended("Bounding Nr by inplane resolution");
    }
    
    // lucky numbers
    if (useLuckyNumbers) {
        timer::gPreloopTimer.begin("Rounding Nr up to FFTW lucky numbers");
        nrMe = nrMe.unaryExpr([](int nr) {
            return NrField::nextLuckyNumber(nr);
        });
        timer::gPreloopTimer.ended("Rounding Nr up to FFTW lucky numbers");
    }
    
    // assemble Nr on ranks
    timer::gPreloopTimer.begin("Assembling Nr on ranks");
    // allocate
    mySuperOnly().mNodalNr = eigen::IColX::Zero(getNumNodes());
    // copy to whole
    mySuperOnly().mNodalNr(idR, Eigen::all) = nrMe;
    // assemble by summation
    mpi::sumEigen(mySuperOnly().mNodalNr);
    // verbose memory
    timer::gPreloopTimer.message
    (eigen_tools::memoryInfo(mySuperOnly().mNodalNr, "nodal Nr (super-only)"));
    timer::gPreloopTimer.ended("Assembling Nr on ranks");
}

// verbose Nr
std::string ExodusMesh::
verboseNr(bool boundByInplane, bool useLuckyNumbers) const {
    std::stringstream ss;
    if (mpi::root()) {
        const eigen::IColX &nodalNr = mySuperOnly().mNodalNr;
        ss << bstring::boxTitle("Computed Nr on Mesh");
        ss << bstring::boxEquals(0, 6, "min Nr", nodalNr.minCoeff());
        ss << bstring::boxEquals(0, 6, "max Nr", nodalNr.maxCoeff());
        long sum = nodalNr.cast<long>().sum();
        int mean = (int)round(1. * sum / nodalNr.rows());
        ss << bstring::boxEquals(0, 6, "ave Nr", mean);
        ss << bstring::boxEquals(0, 6, "sum Nr", sum);
        if (boundByInplane) {
            ss << "* Nr has been limited by inplane resolution.\n";
        }
        if (useLuckyNumbers) {
            ss << "* Nr has been rounded up to FFTW lucky numbers.\n";
        }
        ss << bstring::boxBaseline() << "\n\n";
    }
    return ss.str();
}

// set access to super-only
ExodusMesh::ExodusSuperOnly &ExodusMesh::mySuperOnly() {
    if (!mpi::super()) {
        throw std::runtime_error("ExodusMesh::mySuperOnly || "
                                 "Access to super-only variables from "
                                 "a non-super rank.");
    }
    return mSuperOnly__AccessOnlyBy__mySuperOnly;
}

// get access to super-only
const ExodusMesh::ExodusSuperOnly &ExodusMesh::mySuperOnly() const {
    if (!mpi::super()) {
        throw std::runtime_error("ExodusMesh::mySuperOnly || "
                                 "Access to super-only variables from "
                                 "a non-super rank.");
    }
    return mSuperOnly__AccessOnlyBy__mySuperOnly;
}

// get boundary info for ABC
// input: "RIGHT", "BOTTOM" or "TOP"
// output: 1) whether mesh contains this boundary (return)
//         2) location of this boundary (outer)
//         3) mesh span in the normal direction of this boundary (meshSpan)
bool ExodusMesh::
boundaryInfoABC(const std::string &boundaryKey,
                double &outer, double &meshSpan) const {
    // side name
    std::string ssName;
    if (boundaryKey == "RIGHT") {
        ssName = mKeyRightSS;
    } else if (boundaryKey == "BOTTOM") {
        ssName = mKeyBottomSS;
    } else if (boundaryKey == "TOP") {
        ssName = mKeyTopSS;
    } else {
        throw std::runtime_error("ExodusMesh::boundaryInfoABC || "
                                 "Invalid boundary key: " + boundaryKey);
    }
    
    // mesh does not contain this boundary
    if (mSideSets.find(ssName) == mSideSets.end()) {
        return false;
    }
    
    // outer and mesh span
    if (boundaryKey == "BOTTOM") {
        outer = geodesy::getInnerRadius();
        // note the sign here
        meshSpan = geodesy::getInnerRadius() - geodesy::getOuterRadius();
    } else if (boundaryKey == "TOP") {
        outer = geodesy::getOuterRadius();
        meshSpan = geodesy::getOuterRadius() - geodesy::getInnerRadius();
    } else { // "RIGHT"
        if (mpi::root()) {
            if (geodesy::isCartesian()) {
                outer = mySuperOnly().mNodalCoords.col(0).maxCoeff();
            } else {
                const eigen::DMatX2_RM &rt =
                geodesy::sz2rtheta(mySuperOnly().mNodalCoords, true);
                outer = rt.col(1).maxCoeff();
            }
        }
        mpi::bcast(outer);
        meshSpan = outer;
    }
    return true;
}
