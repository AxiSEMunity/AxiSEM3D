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
    
    // data_time (shared by ranks)
    if (mpi::rank() == mRankWithMaxNumStations) {
        createFileStream(gdir + "/data_time", precisionTime);
    }
    
    // data
    if (mStationCentric) {
        for (const std::string &stKey: stKeys) {
            createFileStream(gdir + "/" + stKey, precisionWave);
        }
        // list_channel.info
        if (mpi::rank() == mRankWithMaxNumStations) {
            std::ofstream fout(gdir + "/list_channel.info");
            if (!fout) {
                throw std::runtime_error("StationIO_Ascii::initialize || "
                                         "Error creating channel index file: ||"
                                         + gdir + "/list_channel.info");
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
        // list_station.info
        std::ofstream fout(rdir + "/list_station.info");
        if (!fout) {
            throw std::runtime_error("StationIO_Ascii::initialize || "
                                     "Error creating station index file: ||"
                                     + rdir + "/list_station.info");
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
    int nst = (int)bufferFields.dimensions()[0];
    if (nst == 0) {
        return;
    }
    
    // time on max rank
    int tfile = 0;
    if (mpi::rank() == mRankWithMaxNumStations) {
        Eigen::internal::set_is_malloc_allowed(true);
        (*mFileStreams[0]) << bufferTime.topRows(bufferLine) << "\n";
        if (mFlush) {
            mFileStreams[0]->flush();
        }
        Eigen::internal::set_is_malloc_allowed(false);
        tfile = 1;
    }
    
    // wavefields
    int nch = (int)bufferFields.dimensions()[1];
    static const eigen::IArray2 shuffle = {1, 0};
    Eigen::internal::set_is_malloc_allowed(true);
    if (mStationCentric) {
        eigen::IArray3 loc = {0, 0, 0};
        eigen::IArray3 len = {1, nch, bufferLine};
        eigen::IArray2 shape = {nch, bufferLine};
        for (int ist = 0; ist < nst; ist++) {
            loc[0] = ist;
            (*mFileStreams[ist + tfile]) << bufferFields.slice(loc, len).
            reshape(shape).shuffle(shuffle) << "\n";
            if (mFlush) {
                mFileStreams[ist + tfile]->flush();
            }
        }
    } else {
        eigen::IArray3 loc = {0, 0, 0};
        eigen::IArray3 len = {nst, 1, bufferLine};
        eigen::IArray2 shape = {nst, bufferLine};
        for (int ich = 0; ich < nch; ich++) {
            loc[1] = ich;
            (*mFileStreams[ich + tfile]) << bufferFields.slice(loc, len).
            reshape(shape).shuffle(shuffle) << "\n";
            if (mFlush) {
                mFileStreams[ich + tfile]->flush();
            }
        }
    }
    Eigen::internal::set_is_malloc_allowed(false);
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
