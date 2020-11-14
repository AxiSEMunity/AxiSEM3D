//
//  NewmarkTimeScheme.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/30/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  Newmark time scheme

#include "NewmarkTimeScheme.hpp"
#include "SolidPoint.hpp"
#include "FluidPoint.hpp"
#include "Domain.hpp"
#include "timer.hpp"

// solve
void NewmarkTimeScheme::solve() const {
    // loop timers
    std::vector<std::string> timerNames = {
        "TOTAL", "SOURCE", "STIFFNESS", "MASS_TERM", "BOUNDARIES",
        "TIME_MARCH", "MPI_COMM", "WAVE_OUTPUT", "SCANNING"};
    std::map<std::string, SimpleTimer> timers;
    for (const std::string &name: timerNames) {
        timers.insert({name, SimpleTimer()});
    }
    
    // verbose
    verboseBegin("NEWMARK");
    
    // start timer
    timers.at("TOTAL").resume();
    
    // disable eigen malloc from now on
    Eigen::internal::set_is_malloc_allowed(false);
    
    // simulation time
    double t = mT0;
    
    ////////////////////////// loop //////////////////////////
    for (int tstep = 0; tstep < mNumTimeSteps; tstep++) {
        // apply source
        timers.at("SOURCE").resume();
        mDomain->applySources(tstep, t);
        timers.at("SOURCE").pause();
        
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
        
        // wavefield output: record and dump
        timers.at("WAVE_OUTPUT").resume();
        mDomain->recordOutput(tstep, t);
        timers.at("WAVE_OUTPUT").pause();
        
        // wavefield scanning
        timers.at("SCANNING").resume();
        mDomain->doScanning(tstep);
        timers.at("SCANNING").pause();
        
        // stability
        if (tstep % mStabilityInterval == 0) {
            mDomain->checkStability(tstep + 1, t, mDt);
        }
        
        // verbose
        verboseIter(timers.at("TOTAL").elapsedTotal(), tstep + 1, t);
        
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
        update(mDomain->getSolidPoints());
        update(mDomain->getFluidPoints());
        timers.at("TIME_MARCH").pause();
        
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
    verboseEnd("NEWMARK", timers.at("TOTAL").elapsedTotal(), t);
}
