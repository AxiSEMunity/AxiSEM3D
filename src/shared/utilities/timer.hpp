//
//  timer.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/9/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  Measuring wall-clock runtime of a multi-level process

#ifndef timer_hpp
#define timer_hpp

// simple timer
#include <vector>
#include <chrono>

// multiple timer
#include <fstream>

/////////////////////// simple timer ///////////////////////
class SimpleTimer {
public:
    // start
    void start() {
        clear();
        resume();
    }
    
    // pause
    void pause() {
        if (running()) {
            mTimePoints.push_back(std::chrono::high_resolution_clock::now());
        }
    }
    
    // resume
    void resume() {
        if (!running()) {
            mTimePoints.push_back(std::chrono::high_resolution_clock::now());
        }
    }
    
    // clear
    void clear() {
        mTimePoints.clear();
    }
    
    // running status
    bool running() const {
        return mTimePoints.size() % 2 == 1;
    }
    
    // elapsed seconds at each pause
    std::vector<double> elapsedAll() const;
    
    // total elapsed seconds
    double elapsedTotal() const;
    
    // clock resolution
    static double getClockResolution();
    
private:
    // recorded time points
    std::vector<std::chrono::high_resolution_clock::time_point> mTimePoints;
};


/////////////////////// multi-level timer ///////////////////////
class MultilevelTimer {
public:
    // constructor
    MultilevelTimer(int nLevels);
    
    // destructor
    ~MultilevelTimer() {
        deactivate();
    }
    
    // activate
    void activate(const std::string &fileName);
    
    // deactivate
    void deactivate();
    
    // activated
    bool isActivated() const {
        return mFile.is_open();
    }
    
    // begin a sub-process
    void begin(const std::string &procName, char indentFill = '.');
    
    // messaging during a sub-process
    void message(const std::string &what, char indentFill = '.');
    
    // end a sub-process
    void ended(const std::string &procName, char indentFill = '.');
    
    // pause during a sub-process
    void pause();
    
    // resume during a sub-process
    void resume();
    
    // current-level elapsed
    double elapsedCurrentLevel() const;
    
private:
    // by constructor
    std::vector<SimpleTimer> mTimers;
    
    // by activate
    std::ofstream mFile;
    int mCurrentLevel = 0;
};


/////////////////////// global multi-level timer ///////////////////////
namespace timer {
    // global timer
    extern MultilevelTimer gPreloopTimer;
    
    // setup
    void setupPreloopDiagnosis(bool activate);
    
    // verbose
    std::string verbose();
}

#endif /* timer_hpp */
