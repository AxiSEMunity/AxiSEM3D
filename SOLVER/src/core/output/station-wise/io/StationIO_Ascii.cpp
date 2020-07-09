//
//  StationIO_Ascii.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/23/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  Ascii IO for station output
//  only for a small number of stations

#include "StationIO_Ascii.hpp"
#include "io.hpp"
#include "mpi.hpp"

// initialize
void StationIO_Ascii::initialize(const std::string &groupName,
                                 int numRecordSteps,
                                 const std::vector<std::string> &channels,
                                 const std::vector<std::string> &stKeys) {
    // finalize
    finalize();
    
    // base
    StationIO::initialize(groupName, numRecordSteps, channels, stKeys);
    
    // need rank_station.info for channel-centric format
    if (!mStationCentric) {
        createRankStation(groupName, stKeys);
    }
    
    // nothing locally without stations
    if (stKeys.size() == 0) {
        return;
    }
    
    // group dir
    const std::string &gdir = io::gOutputDirectory + "/stations/" + groupName;
    
    // precesion
    int precisionTime = 16;
#ifdef _USE_DOUBLE
    int precisionWave = 16;
#else
    int precisionWave = 8;
#endif
    
    // time_points (shared by ranks)
    if (mpi::rank() == mRankWithMaxNumStations) {
        createFileStream(gdir + "/time_points", precisionTime);
    }
    
    // data
    if (mStationCentric) {
        for (const std::string &stKey: stKeys) {
            createFileStream(gdir + "/" + stKey, precisionWave);
        }
        // channel_order.info
        if (mpi::rank() == mRankWithMaxNumStations) {
            std::ofstream fout(gdir + "/channel_order.info");
            if (!fout) {
                throw std::runtime_error("StationIO_Ascii::initialize || "
                                         "Error creating channel index file: ||"
                                         + gdir + "/channel_order.info");
            }
            for (const std::string &channel: channels) {
                fout << channel << "\n";
            }
            fout.close();
        }
    } else {
        // need to seperate ranks
        const std::string &rdir = gdir + "/dir_rank" + mpi::strRank();
        io::mkdir(rdir);
        for (const std::string &channel: channels) {
            createFileStream(rdir + "/" + channel, precisionWave);
        }
        // station_order.info
        std::ofstream fout(rdir + "/station_order.info");
        if (!fout) {
            throw std::runtime_error("StationIO_Ascii::initialize || "
                                     "Error creating station index file: ||"
                                     + rdir + "/station_order.info");
        }
        for (const std::string &stKey: stKeys) {
            fout << stKey << "\n";
        }
        fout.close();
    }
}

// finalize
void StationIO_Ascii::finalize() {
    // close files
    for (const std::unique_ptr<std::ofstream> &fs: mFileStreams) {
        fs->close();
    }
    mFileStreams.clear();
}

// dump to file
void StationIO_Ascii::dumpToFile(const eigen::DColX &bufferTime,
                                 const eigen::RTensor3 &bufferFields,
                                 int bufferLine) {
    // no station
    int nst = (int)bufferFields.dimensions()[2];
    if (nst == 0) {
        return;
    }
    
    // no line, avoid redundant "\n"
    if (bufferLine == 0) {
        return;
    }
    
    // time on max rank
    int tfile = 0;
    if (mpi::rank() == mRankWithMaxNumStations) {
        Eigen::internal::set_is_malloc_allowed(true);
        (*mFileStreams[0]) << bufferTime.topRows(bufferLine) << "\n";
        Eigen::internal::set_is_malloc_allowed(false);
        tfile = 1;
    }
    
    // wavefields
    int nch = (int)bufferFields.dimensions()[1];
    if (mStationCentric) {
        for (int ist = 0; ist < nst; ist++) {
            eigen::IArray3 loc = {0, 0, ist};
            eigen::IArray3 len = {bufferLine, nch, 1};
            Eigen::internal::set_is_malloc_allowed(true);
            (*mFileStreams[ist + tfile]) << bufferFields.slice(loc, len).
            reshape(eigen::IArray2{bufferLine, nch}) << "\n";
            Eigen::internal::set_is_malloc_allowed(false);
        }
    } else {
        for (int ich = 0; ich < nch; ich++) {
            eigen::IArray3 loc = {0, ich, 0};
            eigen::IArray3 len = {bufferLine, 1, nst};
            Eigen::internal::set_is_malloc_allowed(true);
            (*mFileStreams[ich + tfile]) << bufferFields.slice(loc, len).
            reshape(eigen::IArray2{bufferLine, nst}) << "\n";
            Eigen::internal::set_is_malloc_allowed(false);
        }
    }
}

// create file
void StationIO_Ascii::
createFileStream(const std::string &varKey, int precision) {
    const std::string &fname = varKey + ".ascii";
    std::unique_ptr<std::ofstream> fs = std::make_unique<std::ofstream>(fname);
    if (!(*fs)) {
        throw std::runtime_error("StationIO_Ascii::createFileStream || "
                                 "Error creating ascii output file: || "
                                 + fname + " || "
                                 + "Use NetCDF or channel-centric "
                                 "for a large number of stations.");
    }
    fs->precision(precision);
    mFileStreams.push_back(std::move(fs));
}
