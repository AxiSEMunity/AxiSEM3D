//
//  SymplecticTimeScheme.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/30/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  symplectic time scheme

#include "SymplecticTimeScheme.hpp"
#include "SolidPoint.hpp"
#include "FluidPoint.hpp"
#include "Domain.hpp"
#include "timer.hpp"

// solve
void SymplecticTimeScheme::solve() const {
    // loop timer
    std::vector<std::string> timerNames = {
        "TOTAL", "SOURCE", "STIFFNESS", "MASS_TERM", "BOUNDARIES",
        "TIME_MARCH", "MPI_COMM", "WAVE_OUTPUT", "SCANNING"};
    std::map<std::string, SimpleTimer> timers;
    for (const std::string &name: timerNames) {
        timers.insert({name, SimpleTimer()});
    }
    
    // verbose
    verboseBegin("SYMPLECTIC");
    
    // start timer
    timers.at("TOTAL").resume();
    
    // disable eigen malloc from now on
    Eigen::internal::set_is_malloc_allowed(false);
    
    // simulation time
    double t = mT0;
    
    ////////////////////////// loop //////////////////////////
    for (int tstep = 0; tstep < mNumTimeSteps; tstep++) {
        // wavefield output: record and dump
        // Need to do this here to have the same time-axis with Newmark
        timers.at("WAVE_OUTPUT").resume();
        mDomain->recordOutput(tstep, t);
        timers.at("WAVE_OUTPUT").pause();
        
        // wavefield scanning
        timers.at("SCANNING").resume();
        mDomain->doScanning(tstep);
        timers.at("SCANNING").pause();
        
        // apply source
        timers.at("SOURCE").resume();
        mDomain->applySources(tstep, t);
        timers.at("SOURCE").pause();
        
        // launch
        timers.at("TIME_MARCH").resume();
        launch(mDomain->getSolidPoints());
        launch(mDomain->getFluidPoints());
        timers.at("TIME_MARCH").pause();
        
        // sub iteration
        for (int isub = 1; isub <= 4; isub++) {
            // compute stiffness
            timers.at("STIFFNESS").resume();
            mDomain->computeStiffness();
            timers.at("STIFFNESS").pause();
            
            // boundary conditions before assembling stiffness
            timers.at("BOUNDARIES").resume();
            mDomain->applyBC_BeforeAssemblingStiff();
            timers.at("BOUNDARIES").pause();
            
            // assemble phase 1: gather + send + recv
            timers.at("MPI_COMM").resume();
            mDomain->mpiGatherSendRecv();
            timers.at("MPI_COMM").pause();
            
            // stability
            if (tstep % mStabilityInterval == 0 && isub == 1) {
                mDomain->checkStability(tstep + 1, t, mDt);
            }
            
            // assemble phase 2: wait + scatter
            timers.at("MPI_COMM").resume();
            mDomain->mpiWaitScatter();
            timers.at("MPI_COMM").pause();
            
            // boundary conditions after assembling stiffness
            timers.at("BOUNDARIES").resume();
            mDomain->applyBC_AfterAssemblingStiff();
            timers.at("BOUNDARIES").pause();
            
            // stiff => accel
            timers.at("MASS_TERM").resume();
            mDomain->computeStiffToAccel();
            timers.at("MASS_TERM").pause();
            
            // boundary conditions after computing acceleration
            timers.at("BOUNDARIES").resume();
            mDomain->applyBC_AfterComputingAccel();
            timers.at("BOUNDARIES").pause();
            
            // update fields
            timers.at("TIME_MARCH").resume();
            update(mDomain->getSolidPoints(), isub);
            update(mDomain->getFluidPoints(), isub);
            timers.at("TIME_MARCH").pause();
        } // end of sub iteration
        
        // verbose
        verboseIter(timers.at("TOTAL").elapsedTotal(), tstep + 1, t);
        
        // march
        t += mDt;
    } // end of time loop
    
    // wavefield output: dump remaining buffers and finalize
    timers.at("WAVE_OUTPUT").resume();
    mDomain->dumpOutput();
    mDomain->finalizeOutput();
    timers.at("WAVE_OUTPUT").pause();
    
    // enable eigen malloc from now on
    Eigen::internal::set_is_malloc_allowed(true);
    
    // wavefield scanning
    timers.at("SCANNING").resume();
    mDomain->reportScanning();
    timers.at("SCANNING").pause();
    
    // total
    timers.at("TOTAL").pause();
    
    // report timers
    reportLoopTimers(timerNames, timers);
    
    // verbose
    verboseEnd("SYMPLECTIC", timers.at("TOTAL").elapsedTotal(), t);
}
