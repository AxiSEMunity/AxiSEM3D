//
//  ElementOutput.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 7/30/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  element output: generator of ElementOpGroup in core

#include "ElementOutput.hpp"
#include "SE_Model.hpp"
#include "vicinity.hpp"
#include "inparam.hpp"
#include "mpi.hpp"
#include "timer.hpp"
#include "io.hpp"

// in core
#include "ElementIO_NetCDF.hpp"
#include "ElementIO_ParNetCDF.hpp"
#include "ElementOpGroup.hpp"
#include "Domain.hpp"

using namespace bstring;

// build from inparam
std::shared_ptr<const ElementOutput>
ElementOutput::buildInparam(int gindex, double dt,
                            double tmin_simu, double tmax_simu) {
    // group name and key in inparam
    const InparamYAML &gm = inparam::gInparamOutput;
    std::string root = "list_of_element_groups";
    const std::string &groupName = gm.get<std::string>
    (root + ":{" + toString(gindex) + "}");
    root += ":[" + toString(gindex) + "]:" + groupName;
    
    /////// elements ///////
    const std::string &roote = root + ":elements";
    double minR = 0., maxR = 0.;
    double minZ = 0., maxZ = 0.;
    minR = gm.get<double>(roote + ":horizontal_range:[0]");
    maxR = gm.get<double>(roote + ":horizontal_range:[1]");
    minZ = gm.get<double>(roote + ":vertical_range:[0]");
    maxZ = gm.get<double>(roote + ":vertical_range:[1]");
    
    /////// inplane ///////
    // edge
    const std::string &rooti = root + ":inplane";
    int edgeDim = gm.getWithLimits<int>(rooti + ":edge_dimension", {
        {"HORIZONTAL", 0}, {"VERTICAL", 1}, {"BOTH", -1}});
    double edgeCoord = 0.;
    if (edgeDim >= 0) {
        edgeCoord = gm.get<double>(rooti + ":edge_position");
    }
    // points
    std::vector<int> ipols;
    if (gm.get<std::string>(rooti + ":GLL_points_one_edge") == "FULL") {
        ipols.assign(spectral::nPED, 0);
        for (int ipol = 0; ipol < spectral::nPED; ipol++) {
            ipols[ipol] = ipol;
        }
    } else {
        ipols = gm.getVector<int>(rooti + ":GLL_points_one_edge");
    }
    // sort and check
    vector_tools::sortUnique(ipols);
    if (ipols.front() < 0 || ipols.back() > spectral::nPol) {
        throw std::runtime_error
        ("ElementOutput::buildInparam || elements:GLL_points_one_edge must be "
         "ranged between 0 and " + bstring::toString(spectral::nPol) + ".");
    }
    if (ipols.size() == 0) {
        throw std::runtime_error
        ("ElementOutput::buildInparam || elements:GLL_points_one_edge cannot "
         "be empty.");
    }
    
    /////// azimuthal ///////
    const std::string &roota = root + ":azimuthal";
    std::vector<double> phis, lats, lons;
    // phi
    int nphi = gm.get<int>(roota + ":phi_list:[?]");
    for (int iphi = 0; iphi < nphi; iphi++) {
        double phi =
        gm.getWithBounds(roota + ":phi_list:[" + toString(iphi) + "]",
                         0., 2. * numerical::dPi);
        phis.push_back(phi);
    }
    // lat lon
    int nll = gm.get<int>(roota + ":lat_lon_list:[?]");
    for (int ill = 0; ill < nll; ill++) {
        const std::string &key =
        roota + ":lat_lon_list:[" + toString(ill) + "]";
        lats.push_back(gm.getWithBounds(key + ":[0]", -90., 90.));
        lons.push_back(gm.getWithBounds(key + ":[1]", -90., 90.));
    }
    // na space
    int naSpace = gm.getWithBounds<int>(roota + ":na_space", 1);
    
    // fields
    const std::string &rootw = root + ":wavefields";
    // ENZ and xyz make no sense
    channel::WavefieldCS wcs = gm.getWithLimits<channel::WavefieldCS>
    (rootw + ":coordinate_frame", {
        {"spz", channel::WavefieldCS::spz},
        {"RTZ", channel::WavefieldCS::RTZ}});
    bool fluid = gm.getWithLimits<bool>(rootw + ":medium", {
        {"SOLID", false}, {"FLUID", true}});
    const std::vector<std::string> userChannels =
    gm.getVector<std::string>(rootw + ":channels");
    
    // temporal
    const std::string &roott = root + ":temporal";
    // parse sampling period
    const std::string &strSP = gm.get<std::string>(roott + ":sampling_period");
    double samplingPeriod = 0.;
    if (strSP == "DT") {
        samplingPeriod = dt;
    } else if (strSP.substr(0, 3) == "DTx") {
        int n = cast<int>(strSP.substr(3), "StationOutput::buildInparam");
        samplingPeriod = dt * n;
    } else {
        samplingPeriod = cast<double>(strSP, "StationOutput::buildInparam");
    }
    // time window
    double tmin = tmin_simu;
    double tmax = tmax_simu;
    if (gm.get<std::string>(roott + ":time_window") != "FULL") {
        tmin = std::max(tmin_simu, gm.get<double>(roott + ":time_window:[0]"));
        tmax = std::min(tmax_simu, gm.get<double>(roott + ":time_window:[1]"));
    }
    
    // file options
    const std::string &rootf = root + ":file_options";
    // buffer size
    int bufferSize = gm.getWithBounds(rootf + ":buffer_size", 1);
    // flush
    bool flush = gm.get<bool>(rootf + ":flush");
    
    // construct and return
    return std::make_shared
    <const ElementOutput>(groupName, minR, maxR, minZ, maxZ,
                          edgeDim, edgeCoord, ipols,
                          phis, lats, lons, naSpace,
                          wcs, fluid, userChannels,
                          samplingPeriod, tmin, tmax,
                          bufferSize, flush);
}

// verbose
std::string ElementOutput::
verbose(double dt, int numRecordSteps, int npnts,
        int nphis, int numElements) const {
    // title
    std::stringstream ss;
    ss << boxSubTitle(0, mGroupName + " ", '~');
    int width = 24;
    
    //////// elements ////////
    ss << boxSubTitle(2, "Elements");
    ss << boxEquals(4, width, "horizontal range", range(mMinR, mMaxR));
    ss << boxEquals(4, width, "vertical range", range(mMinZ, mMaxZ));
    
    //////// inplane ////////
    ss << boxSubTitle(2, "Inplane sampling");
    if (mEdgeDim == 0) {
        ss << boxEquals(4, width, "edge dimension", "horizontal");
        ss << boxEquals(4, width, "edge position", mEdgeCoord);
        ss << boxEquals(4, width, "# element edges found", numElements);
    } else if (mEdgeDim == 1) {
        ss << boxEquals(4, width, "edge dimension", "vertical");
        ss << boxEquals(4, width, "edge position", mEdgeCoord);
        ss << boxEquals(4, width, "# element edges found", numElements);
    } else {
        ss << boxEquals(4, width, "edge dimension", "both");
        ss << boxEquals(4, width, "# elements found", numElements);
    }
    ss << boxEquals(4, width, "GLL points on one edge", mIPols, "=", true);
    ss << boxEquals(4, width, "# GLL points per element", npnts);
    
    //////// azimuthal ////////
    ss << boxSubTitle(2, "Azimuthal sampling");
    if (nphis == 0) {
        ss << "    * Recording all Fourier series coefficients.\n";
        ss << boxEquals(6, width - 2, "grid space for Nr-storage", mNaSpace);
    } else {
        ss << "    * Recording on " << nphis << " azimuthal slices:\n";
        for (int iphi = 0; iphi < mPhis.size(); iphi++) {
            ss << boxEquals(6, width - 2, toString(iphi + 1) + ") azimuth",
                            mPhis[iphi]);
        }
        for (int ilat = 0; ilat < mLats.size(); ilat++) {
            const std::string &key =
            toString(mPhis.size() + ilat + 1) + ") latitude, longitude";
            const std::string &val =
            toString(mLats[ilat]) + ", " + toString(mLons[ilat]);
            ss << boxEquals(6, width - 2, key, val);
        }
    }
    
    //////// fields ////////
    ss << boxSubTitle(2, "Wavefields");
    ss << boxEquals(4, width, "output coordinate system",
                    channel::WavefieldCS_Str.at(mWCS));
    ss << boxEquals(4, width, "medium type", mFluid ? "fluid" : "solid");
    // channel
    std::vector<std::string> stdChs;
    if (mFluid) {
        channel::fluid::ChannelOptions chops(mWCS, mUserChannels);
        for (int ch: chops.mStdChannels) {
            stdChs.push_back(std::get<0>(channel::fluid::gChannelMap.at(ch)));
        }
    } else {
        channel::solid::ChannelOptions chops(mWCS, mUserChannels);
        for (int ch: chops.mStdChannels) {
            stdChs.push_back(std::get<0>(channel::solid::gChannelMap.at(ch)));
        }
    }
    ss << boxEquals(4, width, "user-specified channels", mUserChannels,
                    "=", true);
    ss << boxEquals(4, width, "standardized channels", stdChs, "=", true);
    
    //////// temporal ////////
    ss << boxSubTitle(2, "Temporal");
    // smapling interval
    int sampleIntv = (int)(mSamplingPeriod / dt);
    // mSamplingPeriod < dt
    if (sampleIntv == 0) {
        sampleIntv = 1;
    }
    ss << "    sampling period\n";
    ss << boxEquals(6, width - 2, "user-specified", mSamplingPeriod);
    // check overflow with very large mSamplingPeriod
    if (sampleIntv < 0) {
        ss << boxEquals(6, width - 1, "rounded to Δt", "Inf (1 sample)");
        ss << boxEquals(4, width, "# time steps per sample", "Inf (1 sample)");
    } else {
        ss << boxEquals(6, width - 1, "rounded to Δt", sampleIntv * dt);
        ss << boxEquals(4, width, "# time steps per sample", sampleIntv);
    }
    ss << boxEquals(4, width, "time window", range(mTmin, mTmax));
    
    //////// file options ////////
    ss << boxSubTitle(2, "File options");
    ss << boxEquals(4, width, "buffer size", mBufferSize);
    ss << boxEquals(4, width, "flush file", mFlush);
    
    //////// dimensions ////////
    ss << boxSubTitle(2, "Dimensions");
    if (mEdgeDim >= 0) {
        ss << boxEquals(4, width, "# element edges", numElements);
    } else {
        ss << boxEquals(4, width, "# elements", numElements);
    }
    if (nphis == 0) {
        ss << boxEquals(4, width, "# Fourier series", "variable");
    } else {
        ss << boxEquals(4, width, "# azimuthal slices", nphis);
    }
    ss << boxEquals(4, width, "# GLL points", npnts);
    ss << boxEquals(4, width, "# channels", stdChs.size());
    ss << boxEquals(4, width, "# time steps", numRecordSteps);
    return ss.str();
    return ss.str();
}

// release ElementOp to domain
void ElementOutput::release(const SE_Model &sem, Domain &domain, double dt,
                            double tmin_simu, double tmax_simu,
                            int nTotalSteps, double distTol) {
    // group count
    int groupCount =
    inparam::gInparamOutput.get<int>("list_of_element_groups:[?]");
    
    // verbose
    std::stringstream ss;
    if (io::gVerbose != io::VerboseLevel::None && mpi::root()) {
        ss << boxTitle("Element groups");
    }
    
    // go over groups
    std::set<std::string> groupNames;
    for (int gindex = 0; gindex < groupCount; gindex++) {
        timer::gPreloopTimer.begin("Element group " + toString(gindex));
        
        //////////// inparam ////////////
        timer::gPreloopTimer.begin("Building element group from inparam");
        std::shared_ptr<const ElementOutput> elgrp =
        buildInparam(gindex, dt, tmin_simu, tmax_simu);
        timer::gPreloopTimer.message("group name = " + elgrp->mGroupName);
        if (groupNames.find(elgrp->mGroupName) != groupNames.end()) {
            throw std::runtime_error("ElementOutput::release || "
                                     "Element group name must be unique. || "
                                     "Duplicated element group name: " +
                                     elgrp->mGroupName);
        }
        groupNames.insert(elgrp->mGroupName);
        timer::gPreloopTimer.ended("Building element group from inparam");
        
        //////////// select elements ////////////
        timer::gPreloopTimer.begin("Selecting elements in mesh");
        std::vector<int> opQuads;
        const std::vector<Quad> &quads = sem.getQuads();
        for (int iquad = 0; iquad < quads.size(); iquad++) {
            // check coordinate ranges
            const eigen::DCol2 &sz = quads[iquad].getNodalSZ().rowwise().mean();
            const eigen::DCol2 &RZ = geodesy::isCartesian() ? sz :
            geodesy::sz2rtheta(sz, false, 0, 1, 1, 0);
            if (elgrp->mMinR < RZ(0) && RZ(0) < elgrp->mMaxR &&
                elgrp->mMinZ < RZ(1) && RZ(1) < elgrp->mMaxZ &&
                quads[iquad].fluid() == elgrp->mFluid) {
                opQuads.push_back(iquad);
            }
        }
        timer::gPreloopTimer.ended("Selecting elements in mesh");
        
        //////////// edge and npnts ////////////
        timer::gPreloopTimer.begin("Handling inplane sampling");
        std::vector<int> opQuadsUse;
        std::vector<std::vector<int>> ipntsElem;
        std::vector<int> ipntsElemUniform;
        if (elgrp->mEdgeDim < 0) {
            // both edges
            opQuadsUse = opQuads;
            for (int ipol: elgrp->mIPols) {
                for (int jpol: elgrp->mIPols) {
                    int ipnt = ipol * spectral::nPED + jpol;
                    ipntsElemUniform.push_back(ipnt);
                }
            }
        } else {
            // one edge
            for (int qTag: opQuads) {
                int edge = -1;
                // compute coordinates
                const eigen::DMat24 &sz = quads[qTag].getNodalSZ();
                const eigen::DMat24 &RZ = geodesy::isCartesian() ? sz :
                geodesy::sz2rtheta(sz, false, 0, 1, 1, 0);
                // coordinate to check
                eigen::DRow4 x;
                double xmin = 0, xmax = 0;
                if (elgrp->mEdgeDim == 0 && !geodesy::isCartesian()) {
                    // theta to arc length
                    double rmean = RZ.row(1).mean();
                    x = RZ.row(elgrp->mEdgeDim) * rmean;
                    xmin = elgrp->mEdgeCoord * rmean - distTol;
                    xmax = elgrp->mEdgeCoord * rmean + distTol;
                } else {
                    x = RZ.row(elgrp->mEdgeDim);
                    xmin = elgrp->mEdgeCoord - distTol;
                    xmax = elgrp->mEdgeCoord + distTol;
                }
                // check coordinate
                for (int iedge = 0; iedge < 4; iedge++) {
                    int jedge = (iedge == 3) ? 0 : iedge + 1;
                    if (xmin < x(iedge) && x(iedge)< xmax &&
                        xmin < x(jedge) && x(jedge)< xmax) {
                        edge = iedge;
                        break;
                    }
                }
                // found an element on edge
                if (edge >= 0) {
                    opQuadsUse.push_back(qTag);
                    std::vector<int> ipntsEdge;
                    for (int ipol: elgrp->mIPols) {
                        int ipnt = vicinity::constants::gEdgeIPnt[edge][ipol];
                        ipntsEdge.push_back(ipnt);
                    }
                    ipntsElem.push_back(ipntsEdge);
                }
            }
        }
        timer::gPreloopTimer.ended("Handling inplane sampling");
        
        //////////// azimuth ////////////
        timer::gPreloopTimer.begin("Handling azimuthal sampling");
        // given phi
        std::vector<double> phisToUse = elgrp->mPhis;
        // given lat, lon
        int nll = (int)elgrp->mLats.size();
        eigen::DMatXX llr(nll, 3);
        llr.col(0) = Eigen::Map<const eigen::DColX>(elgrp->mLats.data(), nll);
        llr.col(1) = Eigen::Map<const eigen::DColX>(elgrp->mLons.data(), nll);
        llr.col(2).fill(geodesy::getOuterRadius());
        const eigen::DMatXX &spz = geodesy::llr2spz(llr, true);
        for (int iphi = 0; iphi < nll; iphi++) {
            phisToUse.push_back(spz(iphi, 1));
        }
        timer::gPreloopTimer.ended("Handling azimuthal sampling");
        
        //////////// release ////////////
        timer::gPreloopTimer.begin("Releasing element group");
        // smapling interval
        int sampleIntv = (int)(elgrp->mSamplingPeriod / dt);
        // mSamplingPeriod < dt
        if (sampleIntv == 0) {
            sampleIntv = 1;
        }
        // check overflow with very large mSamplingPeriod
        if (sampleIntv < 0) {
            sampleIntv = std::numeric_limits<int>::max();
        }
        
        // total number of steps recorded
        int nRecSteps = 0;
        for (int tstep = 0; tstep < nTotalSteps; tstep++) {
            if (tstep % sampleIntv > 0) {
                continue;
            }
            double t = tmin_simu + tstep * dt;
            if (t < elgrp->mTmin - dt / 2. || t > elgrp->mTmax + dt / 2.) {
                continue;
            }
            nRecSteps++;
        }
        if (nRecSteps == 0) {
            throw std::runtime_error
            ("ElementOutput::release || "
             "No time step to record in the given time window. ||"
             "Element group name: " + elgrp->mGroupName);
        }
        
        // io
        std::unique_ptr<ElementIO> elementIO = nullptr;
#ifdef _USE_PARALLEL_NETCDF
        elementIO = std::make_unique<ElementIO_ParNetCDF>();
#else
        elementIO = std::make_unique<ElementIO_NetCDF>();
#endif
        // set flush
        elementIO->setFlush(elgrp->mFlush);
        
        // create both solid and fluid groups
        std::unique_ptr<ElementOpGroup<ElementOpFluid>> EGF = nullptr;
        std::unique_ptr<ElementOpGroup<ElementOpSolid>> EGS = nullptr;
        // but only initialize one (elementIO can be used only once)
        int npnts = (elgrp->mEdgeDim < 0 ?
                     (int)elgrp->mIPols.size() * (int)elgrp->mIPols.size() :
                     (int)elgrp->mIPols.size());
        if (elgrp->mFluid) {
            EGF = std::make_unique<ElementOpGroup<ElementOpFluid>>
            (elgrp->mGroupName, nRecSteps, sampleIntv,
             elgrp->mTmin - dt / 2., elgrp->mTmax + dt / 2., elgrp->mBufferSize,
             elgrp->mWCS, elgrp->mUserChannels,
             npnts, phisToUse, elgrp->mNaSpace, elementIO);
        } else {
            EGS = std::make_unique<ElementOpGroup<ElementOpSolid>>
            (elgrp->mGroupName, nRecSteps, sampleIntv,
             elgrp->mTmin - dt / 2., elgrp->mTmax + dt / 2., elgrp->mBufferSize,
             elgrp->mWCS, elgrp->mUserChannels,
             npnts, phisToUse, elgrp->mNaSpace, elementIO);
        }
        
        // add elements to group
        for (int ielem = 0; ielem < opQuadsUse.size(); ielem++) {
            const std::vector<int> &ipnts =
            (elgrp->mEdgeDim < 0) ? ipntsElemUniform : ipntsElem[ielem];
            if (elgrp->mFluid) {
                std::unique_ptr<ElementOpFluid> eop =
                std::make_unique<ElementOpFluid>(ipnts);
                eop->setElement(quads[opQuadsUse[ielem]].getFluidElement());
                EGF->addElementOp(eop);
            } else {
                std::unique_ptr<ElementOpSolid> eop =
                std::make_unique<ElementOpSolid>(ipnts);
                eop->setElement(quads[opQuadsUse[ielem]].getSolidElement());
                EGS->addElementOp(eop);
            }
        }
        
        // add group to domain
        if (elgrp->mFluid) {
            domain.addElementOpGroupInFluid(EGF);
        } else {
            domain.addElementOpGroupInSolid(EGS);
        }
        timer::gPreloopTimer.ended("Releasing element group");
        
        
        //////////// verbose ////////////
        int nelemGlobal = mpi::sum((int)opQuadsUse.size());
        if (io::gVerbose != io::VerboseLevel::None && mpi::root()) {
            timer::gPreloopTimer.begin("Verbosing element group");
            ss << elgrp->verbose(dt, nRecSteps, npnts,
                                 (int)phisToUse.size(), nelemGlobal);
            timer::gPreloopTimer.ended("Verbosing element group");
        }
        
        timer::gPreloopTimer.ended("Element group " + toString(gindex));
    }
    
    // verbose
    if (io::gVerbose != io::VerboseLevel::None && mpi::root()) {
        if (groupCount == 0) {
            ss << "* No element groups in this simulation.\n";
        }
        ss << boxBaseline() << "\n\n";
        io::cout << ss.str();
    }
}
