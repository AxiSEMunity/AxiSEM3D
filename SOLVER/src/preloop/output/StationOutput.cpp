//
//  StationOutput.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 5/25/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  station output: generator of StationGroup in core

#include "StationOutput.hpp"
#include "SE_Model.hpp"
#include "inparam.hpp"
#include "mpi.hpp"
#include "timer.hpp"
#include "io.hpp"

// for computeSPZ
#include "Source.hpp"

// in core
#include "StationIO_Ascii.hpp"
#include "StationIO_NetCDF.hpp"
#include "StationIO_ParNetCDF.hpp"
#include "StationGroup.hpp"
#include "Domain.hpp"

using namespace bstring;

// build from inparam
std::shared_ptr<const StationOutput>
StationOutput::buildInparam(int gindex, double dt,
                            double tmin_simu, double tmax_simu) {
    // group name and key in inparam
    const InparamYAML &gm = inparam::gInparamOutput;
    std::string root = "list_of_station_groups";
    const std::string &groupName = gm.get<std::string>
    (root + ":{" + toString(gindex) + "}");
    root += ":[" + toString(gindex) + "]:" + groupName;
    
    // locations
    const std::string &rootl = root + ":locations";
    const std::string &fileName = gm.get<std::string>(rootl + ":station_file");
    bool sourceCentered = gm.getWithLimits<bool>(rootl + ":horizontal_x1_x2", {
        {"LATITUDE_LONGITUDE", false}, {"DISTANCE_AZIMUTH", true},
        {"XY_CARTESIAN", true}});
    bool xy = gm.getWithLimits<bool>(rootl + ":horizontal_x1_x2", {
        {"LATITUDE_LONGITUDE", false}, {"DISTANCE_AZIMUTH", false},
        {"XY_CARTESIAN", true}});
    bool useDepth = gm.getWithLimits<bool>(rootl + ":vertical_x3", {
        {"RADIUS", false}, {"DEPTH", true}});
    bool ellipticity = false;
    if (!sourceCentered) {
        ellipticity = gm.get<bool>(rootl + ":ellipticity");
    }
    bool depthSolid = false;
    if (useDepth) {
        depthSolid = gm.get<bool>(rootl + ":depth_below_solid_surface");
    }
    bool undulated = gm.get<bool>(rootl + ":undulated_geometry");
    
    // fields
    const std::string &rootw = root + ":wavefields";
    channel::WavefieldCS wcs = gm.getWithLimits<channel::WavefieldCS>
    (rootw + ":coordinate_frame", {
        {"spz", channel::WavefieldCS::spz},
        {"RTZ", channel::WavefieldCS::RTZ},
        {"ENZ", channel::WavefieldCS::ENZ},
        {"xyz", channel::WavefieldCS::xyz}});
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
    // format
    Format format = gm.getWithLimits<Format>(rootf + ":format", {
        {"ASCII_STATION", Format::AsciiStation},
        {"ASCII_CHANNEL", Format::AsciiChannel},
        {"NETCDF", Format::NetCDF}});
    // buffer size
    int bufferSize = gm.getWithBounds(rootf + ":buffer_size", 1);
    // flush
    bool flush = gm.get<bool>(rootf + ":flush");
    
    // construct and return
    return std::make_shared
    <const StationOutput>(groupName, fileName,
                          sourceCentered, xy, ellipticity,
                          useDepth, depthSolid, undulated,
                          wcs, fluid, userChannels,
                          samplingPeriod, tmin, tmax,
                          format, bufferSize, flush);
}

// verbose
std::string StationOutput::
verbose(double dt, int numRecordSteps, int numStations) const {
    // title
    std::stringstream ss;
    ss << boxSubTitle(0, mGroupName + " ", '~');
    int width = mUseDepth ? 25 : 24;
    
    //////// locations ////////
    ss << boxSubTitle(2, "Locations");
    ss << boxEquals(4, width, "station location file", mFileName);
    std::string crd = "(";
    if (mXY) {
        crd += "X, Y, ";
    } else {
        crd += (mSourceCentered ? "distance, azimuth, "
                : "latitude, longitude, ");
    }
    crd += (mUseDepth ? "depth)" : "radius)");
    ss << boxEquals(4, width, "coordinates in file", crd);
    if (!mSourceCentered) {
        ss << boxEquals(4, width, "ellipticity correction", mEllipticity);
    }
    if (mUseDepth) {
        ss << boxEquals(4, width, "depth below solid surface", mDepthSolid);
    }
    ss << boxEquals(4, width, "radial geometry", mUndulatedGeometry ?
                    "undulated" : "reference");
    
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
    if (mFormat == Format::AsciiStation) {
        ss << boxEquals(4, width, "output file format", "ascii (per station)");
    } else if (mFormat == Format::AsciiChannel) {
        ss << boxEquals(4, width, "output file format", "ascii (per channel)");
    } else {
#ifdef _USE_PARALLEL_NETCDF
        ss << boxEquals(4, width, "output file format", "NetCDF (parallel)");
#else
        ss << boxEquals(4, width, "output file format", "NetCDF (serial)");
#endif
    }
    ss << boxEquals(4, width, "buffer size", mBufferSize);
    ss << boxEquals(4, width, "flush file", mFlush);
    
    //////// dimensions ////////
    ss << boxSubTitle(2, "Dimensions");
    ss << boxEquals(4, width, "# stations", numStations);
    ss << boxEquals(4, width, "# channels", stdChs.size());
    ss << boxEquals(4, width, "# time steps", numRecordSteps);
    return ss.str();
}

// release stations to domain
void StationOutput::release(const SE_Model &sem, Domain &domain, double dt,
                            double tmin_simu, double tmax_simu,
                            int nTotalSteps) {
    // group count
    int groupCount =
    inparam::gInparamOutput.get<int>("list_of_station_groups:[?]");
    
    // verbose
    std::stringstream ss;
    if (io::gVerbose != io::VerboseLevel::None && mpi::root()) {
        ss << boxTitle("Station groups");
    }
    
    // go over groups
    std::set<std::string> groupNames;
    for (int gindex = 0; gindex < groupCount; gindex++) {
        timer::gPreloopTimer.begin("Station group " + toString(gindex));
        
        //////////// inparam ////////////
        timer::gPreloopTimer.begin("Building station group from inparam");
        std::shared_ptr<const StationOutput> stgrp =
        buildInparam(gindex, dt, tmin_simu, tmax_simu);
        timer::gPreloopTimer.message("group name = " + stgrp->mGroupName);
        if (groupNames.find(stgrp->mGroupName) != groupNames.end()) {
            throw std::runtime_error("StationOutput::release || "
                                     "Station group name must be unique. || "
                                     "Duplicated station group name: " +
                                     stgrp->mGroupName);
        }
        groupNames.insert(stgrp->mGroupName);
        timer::gPreloopTimer.ended("Building station group from inparam");
        
        
        //////////// read ////////////
        timer::gPreloopTimer.begin("Reading and broadcasting station file");
        std::vector<std::string> stKeys;
        eigen::DMatX3 stCrds;
        // open and read
        if (mpi::root()) {
            // use map to check key uniqueness while recording order in map
            std::map<std::string, std::pair<int, eigen::DRow3>> keyCrds;
            const std::vector<std::string> &lines = readLines
            (io::popInputDir(stgrp->mFileName), "StationOutput::release");
            int numStations = 0;
            for (const std::string &line: lines) {
                const std::vector<std::string> &words = split(line, " \t");
                // empty or comment line
                if (words[0].front() == '#' || words[0] == "") {
                    continue;
                }
                // size check
                if (words.size() != 5 && words.size() != 6) {
                    throw std:: runtime_error
                    ("StationOutput::release || "
                     "Input file for stations must contain five or six columns."
                     " || Ascii file: " + io::popInputDir(stgrp->mFileName) +
                     " || Station group name: " + stgrp->mGroupName);
                }
                // key
                const std::string &stKey = words[1] + "." + words[0];
                if (keyCrds.find(stKey) != keyCrds.end()) {
                    throw std::runtime_error
                    ("StationOutput::release || Station key must be unique. || "
                     "Duplicated station key: " + stKey + " || "
                     "Station group name: " + stgrp->mGroupName);
                }
                // coords
                static eigen::DRow3 crd;
                crd(0) = cast<double>(words[2], "StationOutput::release");
                crd(1) = cast<double>(words[3], "StationOutput::release");
                crd(2) = cast<double>(words.back(), "StationOutput::release");
                keyCrds.insert({stKey, {numStations, crd}});
                // increment number
                numStations++;
            }
            // cast to vector
            stKeys.resize(numStations);
            stCrds.resize(numStations, 3);
            for (auto it = keyCrds.begin(); it != keyCrds.end(); ++it) {
                int ist = std::get<0>(it->second);
                stKeys[ist] = it->first;
                stCrds.row(ist) = std::get<1>(it->second);
            }
        }
        // broad cast
        mpi::bcast(stKeys);
        mpi::bcastEigen(stCrds);
        timer::gPreloopTimer.ended("Reading and broadcasting station file");
        
        
        //////////// locate ////////////
        timer::gPreloopTimer.begin("Locating stations in mesh");
        std::map<int, int> stationRank;
        std::map<int, int> stationQuad;
        int numStations = (int)stKeys.size();
        for (int ist = 0; ist < numStations; ist++) {
            // compute spz
            const eigen::DRow3 &spz = Source::
            computeSPZ(sem, stCrds.row(ist),
                       stgrp->mSourceCentered, stgrp->mXY, stgrp->mEllipticity,
                       stgrp->mUseDepth, stgrp->mDepthSolid,
                       stgrp->mUndulatedGeometry,
                       "Station group name: " + stgrp->mGroupName, false);
            
            // locate in mesh
            int quadTag =
            sem.locateInplane(spz({0, 2}).transpose(), stgrp->mFluid);
            if (quadTag != -1) {
                stationRank.insert({ist, mpi::rank()});
                stationQuad.insert({ist, quadTag});
            }
        }
        timer::gPreloopTimer.ended("Locating stations in mesh");
        
        // mpi assemble
        // use minimum rank if station is located on more than one ranks
        timer::gPreloopTimer.begin("Assembling and verifying locations");
        mpi::aggregate(stationRank, MPI_ALL, MPI_MIN);
        
        // check if all stations are located
        if (stationRank.size() != numStations) {
            for (int ist = 0; ist < numStations; ist++) {
                if (stationRank.find(ist) == stationRank.end()) {
                    throw std::runtime_error
                    ("StationOutput::release || Error locating station in mesh."
                     " || Station key: " + stKeys[ist] +
                     " || Station group name: " + stgrp->mGroupName);
                }
            }
        }
        timer::gPreloopTimer.ended("Assembling and verifying locations");
        
        
        //////////// release ////////////
        timer::gPreloopTimer.begin("Releasing station group");
        // smapling interval
        int sampleIntv = (int)(stgrp->mSamplingPeriod / dt);
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
            if (t < stgrp->mTmin - dt / 2. || t > stgrp->mTmax + dt / 2.) {
                continue;
            }
            nRecSteps++;
        }
        if (nRecSteps == 0) {
            throw std::runtime_error
            ("StationOutput::release || "
             "No time step to record in the given time window. ||"
             "Station group name: " + stgrp->mGroupName);
        }
        
        // io
        std::unique_ptr<StationIO> stationIO = nullptr;
        if (stgrp->mFormat == Format::AsciiStation) {
            stationIO = std::make_unique<StationIO_Ascii>(true);
        } else if (stgrp->mFormat == Format::AsciiChannel) {
            stationIO = std::make_unique<StationIO_Ascii>(false);
        } else {
#ifdef _USE_PARALLEL_NETCDF
            stationIO = std::make_unique<StationIO_ParNetCDF>();
#else
            stationIO = std::make_unique<StationIO_NetCDF>();
#endif
        }
        // set flush
        stationIO->setFlush(stgrp->mFlush);
        
        // create both solid and fluid groups
        std::unique_ptr<StationGroup<StationFluid>> SGF = nullptr;
        std::unique_ptr<StationGroup<StationSolid>> SGS = nullptr;
        // but only initialize one (stationIO can be used only once)
        if (stgrp->mFluid) {
            SGF = std::make_unique<StationGroup<StationFluid>>
            (stgrp->mGroupName, nRecSteps, sampleIntv,
             stgrp->mTmin - dt / 2., stgrp->mTmax + dt / 2., stgrp->mBufferSize,
             stgrp->mWCS, stgrp->mUserChannels, stationIO);
        } else {
            SGS = std::make_unique<StationGroup<StationSolid>>
            (stgrp->mGroupName, nRecSteps, sampleIntv,
             stgrp->mTmin - dt / 2., stgrp->mTmax + dt / 2., stgrp->mBufferSize,
             stgrp->mWCS, stgrp->mUserChannels, stationIO);
        }
        
        // get quads
        const std::vector<Quad> &quads = sem.getQuads();
        // add stations to group
        for (int ist = 0; ist < numStations; ist++) {
            if (stationRank.at(ist) != mpi::rank()) {
                // handled by another rank
                continue;
            }
            
            // compute spz
            const eigen::DRow3 &spz = Source::
            computeSPZ(sem, stCrds.row(ist),
                       stgrp->mSourceCentered, stgrp->mXY, stgrp->mEllipticity,
                       stgrp->mUseDepth, stgrp->mDepthSolid,
                       stgrp->mUndulatedGeometry,
                       "Station group name: " + stgrp->mGroupName, false);
            
            // theta and baz
            double theta = geodesy::sz2rtheta(spz, true, 0, 2, 2, 0)(0);
            const eigen::DRow3 &llr = geodesy::spz2llr(spz, true, false);
            double baz = geodesy::backAzimuth(llr, true)(0);
            
            // inplane interpolation
            int quadTag = stationQuad.at(ist);
            const eigen::DRowN &inplaneFactor =
            sem.computeInplaneFactor(spz({0, 2}).transpose(), quadTag);
            
            // station
            if (stgrp->mFluid) {
                std::unique_ptr<StationFluid> st =
                std::make_unique<StationFluid>(stKeys[ist], spz(1), theta, baz);
                st->setElement(quads[stationQuad[ist]].getFluidElement(),
                               inplaneFactor);
                SGF->addStation(st);
            } else {
                std::unique_ptr<StationSolid> st =
                std::make_unique<StationSolid>(stKeys[ist], spz(1), theta, baz);
                st->setElement(quads[stationQuad[ist]].getSolidElement(),
                               inplaneFactor);
                SGS->addStation(st);
            }
        }
        
        // add group to domain
        if (stgrp->mFluid) {
            domain.addStationGroupInFluid(SGF);
        } else {
            domain.addStationGroupInSolid(SGS);
        }
        timer::gPreloopTimer.ended("Releasing station group");
        
        
        //////////// verbose ////////////
        if (io::gVerbose != io::VerboseLevel::None && mpi::root()) {
            timer::gPreloopTimer.begin("Verbosing station group");
            ss << stgrp->verbose(dt, nRecSteps, numStations);
            timer::gPreloopTimer.ended("Verbosing station group");
        }
        
        timer::gPreloopTimer.ended("Station group " + toString(gindex));
    }
    
    // verbose
    if (io::gVerbose != io::VerboseLevel::None && mpi::root()) {
        if (groupCount == 0) {
            ss << "* No station groups in this simulation.\n";
        }
        ss << boxBaseline() << "\n\n";
        io::cout << ss.str();
    }
}
