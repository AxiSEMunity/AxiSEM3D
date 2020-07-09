//
//  timer.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/9/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  Measuring wall-clock runtime of a multi-level process

#include "timer.hpp"

// simple timer
#include "vector_tools.hpp"

// multiple timer
#include "io.hpp"
#include "mpi.hpp"
#include "bstring.hpp"

/////////////////////// simple timer ///////////////////////
using namespace std::chrono;

// elapsed seconds at each pause
std::vector<double> SimpleTimer::elapsedAll() const {
    std::vector<double> elapsed;
    int npairs = (int)mTimePoints.size() / 2;
    elapsed.reserve(npairs + 1);
    for (int ipair = 0; ipair < npairs; ipair++) {
        double elap = (duration_cast<duration<double>>
                       (mTimePoints[2 * ipair + 1] -
                        mTimePoints[2 * ipair])).count();
        elapsed.push_back(elap);
    }
    if (running()) {
        double elap = (duration_cast<duration<double>>
                       (high_resolution_clock::now() -
                        mTimePoints.back())).count();
        elapsed.push_back(elap);
    }
    return elapsed;
}

// total elapsed seconds
double SimpleTimer::elapsedTotal() const {
    const std::vector<double> &elapsed = elapsedAll();
    return vector_tools::sum(elapsed);
}

// clock resolution
double SimpleTimer::getClockResolution() {
    SimpleTimer cpu;
    high_resolution_clock::time_point start_time, current_time;
    start_time = high_resolution_clock::now();
    current_time = start_time;
    while (current_time == start_time) {
        current_time = high_resolution_clock::now();
    }
    return (duration_cast<duration<double>>
            (current_time - start_time)).count();
}


/////////////////////// multi-level timer ///////////////////////
// constructor
MultilevelTimer::MultilevelTimer(int nLevels):
mTimers(std::vector<SimpleTimer>(nLevels, SimpleTimer())),
mFile(), mCurrentLevel(0) {
    // nothing
}

// activate
void MultilevelTimer::activate(const std::string &fileName) {
    if (mpi::root()) {
        // deactivate existing
        deactivate();
        // file
        mFile.open(fileName);
        if (!mFile) {
            throw std::runtime_error("MultilevelTimer::activate || "
                                     "Error creating log file: || " +
                                     fileName);
        }
        // timers
        for (SimpleTimer &timer: mTimers) {
            timer.clear();
        }
        mCurrentLevel = 0;
    }
}

// deactivate
void MultilevelTimer::deactivate() {
    if (isActivated()) {
        mFile.close();
    }
}

// begin a sub-process
void MultilevelTimer::begin(const std::string &procName, char indentFill) {
    if (isActivated()) {
        mFile << bstring::filled(mCurrentLevel * 4, indentFill);
        mFile << "<BEGIN>  "<< procName << "\n";
        mFile.flush();
        mTimers[mCurrentLevel].start();
        mCurrentLevel++;
    }
}

// messaging during a sub-process
void MultilevelTimer::message(const std::string &what, char indentFill) {
    if (isActivated()) {
        mFile << bstring::filled(mCurrentLevel * 4, indentFill);
        mFile << "<MESSG>  " << what << "\n";
        mFile.flush();
    }
}

// end a sub-process
void MultilevelTimer::ended(const std::string &procName, char indentFill) {
    if (isActivated()) {
        mCurrentLevel--;
        mTimers[mCurrentLevel].pause();
        mFile << bstring::filled(mCurrentLevel * 4, indentFill);
        mFile << "<ENDED>  " << procName << " [ELAPSED = ";
        mFile << mTimers[mCurrentLevel].elapsedTotal();
        mFile << " sec]" << "\n";
        mFile.flush();
    }
}

// pause during a sub-process
void MultilevelTimer::pause() {
    if (isActivated()) {
        for (int ilvl = 0; ilvl < mCurrentLevel; ilvl++) {
            mTimers[ilvl].pause();
        }
    }
}

// resume during a sub-process
void MultilevelTimer::resume() {
    if (isActivated()) {
        for (int ilvl = 0; ilvl < mCurrentLevel; ilvl++) {
            mTimers[ilvl].resume();
        }
    }
}

// current level elapsed
double MultilevelTimer::elapsedCurrentLevel() const {
    if (isActivated() && mCurrentLevel > 1) {
        return mTimers[mCurrentLevel - 1].elapsedTotal();
    } else {
        // never called begin()
        return 0.;
    }
}


/////////////////////// global multi-level timer ///////////////////////
namespace timer {
    MultilevelTimer gPreloopTimer = MultilevelTimer(10);
    
    // setup
    void setupPreloopDiagnosis(bool activate) {
        if (activate) {
            gPreloopTimer.activate(io::gOutputDirectory +
                                   "/develop/preloop_diagnosis.log");
        }
    }
    
    // verbose
    std::string verbose() {
        std::stringstream ss;
        ss << bstring::boxTitle("Preloop Diagnosis");
        if (gPreloopTimer.isActivated()) {
            ss << "* Preloop diagnosis is turned on.\n";
            ss << "log file  =  output/develop/preloop_diagnosis.log\n";
        } else {
            ss << "* Preloop diagnosis is turned off.\n";
        }
        ss << bstring::boxBaseline() << "\n\n";
        return ss.str();
    }
}
