//
//  io.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/8/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  IO interfaces

#include "io.hpp"

// io::cout
#include <iostream>
#include <fstream>

// mkdir
extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
};
#include "mpi.hpp"

// setup
#include "inparam.hpp"

// verbose
using namespace bstring;

namespace io {
    //////////////////// input/output dirs ////////////////////
    std::string gInputDirectory = "";
    std::string gOutputDirectory = "";
    
    // check dir existence
    bool dirExists(const std::string &path) {
        struct stat info;
        if (stat(path.c_str(), &info) != 0) {
            return false;
        } else if (info.st_mode & S_IFDIR) {
            return true;
        } else {
            return false;
        }
    }
    
    // mkdir
    void mkdir(const std::string &path) {
        if (!dirExists(path)) {
            ::mkdir(path.c_str(), ACCESSPERMS);
        }
    }
    
    // verify input/output dirs under executable dir
    void verifyDirectories(int argc, char *argv[], std::string &warning) {
        // find path of executable
        std::string argv0(argv[0]);
        std::size_t found = argv0.find_last_of("/\\");
        const std::string &execDirectory =
        (found == std::string::npos) ? "." : argv0.substr(0, found);
        
        // input and output
        gInputDirectory = execDirectory + "/input";
        gOutputDirectory = execDirectory + "/output";
        warning = "";
        if (mpi::root()) {
            // check input
            if (!dirExists(gInputDirectory)) {
                throw std::runtime_error("io::verifyDirectories || "
                                         "Missing input directory: || "
                                         + gInputDirectory);
            }
            // check output
            if (dirExists(gOutputDirectory)) {
                // backup the old
                const std::string &backup =
                gOutputDirectory + "__backup@" + bstring::currentDateTime();
                ::rename(gOutputDirectory.c_str(), backup.c_str());
                // warning
                warning = bstring::warning("io::verifyDirectories || "
                                           "Output directory exists; "
                                           "old output renamed to || " +
                                           backup);
            }
            // create output folders
            mkdir(gOutputDirectory);
            mkdir(gOutputDirectory + "/stations");
            mkdir(gOutputDirectory + "/elements");
            mkdir(gOutputDirectory + "/develop");
            mkdir(gOutputDirectory + "/plots");
        }
    }
    
    // pop input dir before filename
    std::string popInputDir(const std::string &fname) {
        if (fname.front() == '/') {
            // absolute path
            return fname;
        }
        return gInputDirectory + "/" + fname;
    }
    
    
    //////////////////// runtime verbose ////////////////////
    // verbose control
    VerboseLevel gVerbose = VerboseLevel::Detailed;
    bool gVerboseWarnings = true;
    
    // setup verbose
    void setupVerbose() {
        // channel
        const std::string &channel =
        inparam::gInparamAdvanced.get<std::string>("verbose:channel");
        if (channel != "STDOUT") {
            // change cout to file
            const std::string &fname = gOutputDirectory + "/" + channel;
            cout.mFileStream = std::make_unique<std::ofstream>(fname);
            if (!*cout.mFileStream) {
                throw std::runtime_error("io::setupVerbose || Error opening or "
                                         "creating stdout file: || " + fname);
            }
            cout.mCoutStream = &(*cout.mFileStream);
        }
        // level
        io::gVerbose = inparam::gInparamAdvanced.
        getWithLimits<io::VerboseLevel>("verbose:level", {
            {"NONE", io::VerboseLevel::None},
            {"ESSENTIAL", io::VerboseLevel::Essential},
            {"DETAILED", io::VerboseLevel::Detailed}});
        // warnings
        io::gVerboseWarnings =
        inparam::gInparamAdvanced.get<bool>("verbose:warnings");
    }
    
    // print welcome messages
    std::string welcome() {
        std::string welc = "\n"
        "{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~}"
        "[                                                            ]"
        "[    A            |i .|'''.|'||''''E'||    ||'  ____'||''|.  ]"
        "[   |||   ... ...... ||..  ' ||  .   |||  |||   ` // ||   || ]"
        "[  |  ||   '|..'  ||  ''|||. ||''|   |'|..'||    //  ||    ||]"
        "[ .''''|.   .x.   ||      '||||      | 'M' ||    \\  ||    ||]"
        "[.|.  .||..|  ||..||.|'...|S.||....|.|. | .||.    3' D|...|' ]"
        "[.............................................   //          ]"
        "[                                               /'     v x.y ]"
        "[                                                            ]"
        "[Copyright (c) 2019 Kuangdai Leng & friends, MIT License     ]"
        "[Source, docs and issues: github.com/kuangdai/AxiSEM-3D      ]"
        "{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~}\n\n";
        
        std::string space = filled((fmt::gBoxWidth - 60) / 2);
        std::string tilde = filled((fmt::gBoxWidth - 60) / 2, '~');
        welc = replace(welc, "{", tilde);
        welc = replace(welc, "}", tilde + "\n");
        welc = replace(welc, "[", space);
        welc = replace(welc, "]", space + "\n");
        welc = replace(welc, "x.y", _VERSION);
        return replace(welc, "\\", "\\\\");
    }
    
    // verbose
    std::string verbose() {
        std::stringstream ss;
        ss << boxTitle("IO");
        ss << boxSubTitle(0, "Directories");
        ss << boxEquals(2, 8, "input", gInputDirectory);
        ss << boxEquals(2, 8, "output", gOutputDirectory);
        ss << boxEquals(2, 8, "source", gProjectDirectory);
        ss << boxSubTitle(0, "Verbose");
        if (gVerbose == VerboseLevel::Essential) {
            ss << boxEquals(2, 8, "level", "essential");
        } else if (gVerbose == VerboseLevel::Detailed) {
            ss << boxEquals(2, 8, "level", "detailed");
        } else {
            ss << boxEquals(2, 8, "level", "none");
        }
        ss << boxEquals(2, 8, "warnings", gVerboseWarnings);
        ss << boxBaseline() << "\n\n";
        return ss.str();
    }
    
    ////////////////////////////// cout on root //////////////////////////////
    mpi_root_cout::mpi_root_cout()
    :mMyWorldRank(-1), mCoutWorldRank(0), mCoutStream(&(std::cout)) {
        // nothing
    }
    mpi_root_cout cout;
}
