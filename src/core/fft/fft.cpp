//
//  fft.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/22/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  global FFT solvers for core

#include "fft.hpp"

// to use template member functions implemented in cpp
#include "SolverFFTW.cpp"

// verbose
#include "bstring.hpp"
#include "mpi.hpp"
#include "timer.hpp"
#include "vector_tools.hpp"

namespace fft {
    // global FFT solvers
    SolverFFTW<Real, 1> gFFT_1;
    SolverFFTW<Real, 3> gFFT_3;
    SolverFFTW<Real, nPEM * 3> gFFT_N3;
    SolverFFTW<Real, nPEM * 6> gFFT_N6;
    SolverFFTW<Real, nPEM * 9> gFFT_N9;
    
    // create plans
    void createPlans(double timeLimitForPlanning) {
        // split time limit
        double tf1 = gFFT_1.timeFactorForPlanning();
        double tf3 = gFFT_3.timeFactorForPlanning();
        double tfN3 = gFFT_N3.timeFactorForPlanning();
        double tfN6 = gFFT_N6.timeFactorForPlanning();
        double tfN9 = gFFT_N9.timeFactorForPlanning();
        double tfSum = std::max(tf1 + tf3 + tfN3 + tfN6 + tfN9,
                                numerical::dEpsilon);
        double timeLimitUnit = timeLimitForPlanning / tfSum;
        
        // 1
        timer::gPreloopTimer.begin("FFT with HOWMANY = 1");
        gFFT_1.createPlans(timeLimitUnit * tf1);
        internal::diagnosticInfo(gFFT_1.getPlannedNRs(), 1);
        timer::gPreloopTimer.ended("FFT with HOWMANY = 1");
        
        // 3
        timer::gPreloopTimer.begin("FFT with HOWMANY = 3");
        gFFT_3.createPlans(timeLimitUnit * tf3);
        internal::diagnosticInfo(gFFT_3.getPlannedNRs(), 3);
        timer::gPreloopTimer.ended("FFT with HOWMANY = 3");
        
        // N3
        timer::gPreloopTimer.begin("FFT with HOWMANY = # GLL per Element * 3");
        gFFT_N3.createPlans(timeLimitUnit * tfN3);
        internal::diagnosticInfo(gFFT_N3.getPlannedNRs(), nPEM * 3);
        timer::gPreloopTimer.ended("FFT with HOWMANY = # GLL per Element * 3");
        
        // N6
        timer::gPreloopTimer.begin("FFT with HOWMANY = # GLL per Element * 6");
        gFFT_N6.createPlans(timeLimitUnit * tfN6);
        internal::diagnosticInfo(gFFT_N6.getPlannedNRs(), nPEM * 6);
        timer::gPreloopTimer.ended("FFT with HOWMANY = # GLL per Element * 6");
        
        // N9
        timer::gPreloopTimer.begin("FFT with HOWMANY = # GLL per Element * 9");
        gFFT_N9.createPlans(timeLimitUnit * tfN9);
        internal::diagnosticInfo(gFFT_N9.getPlannedNRs(), nPEM * 9);
        timer::gPreloopTimer.ended("FFT with HOWMANY = # GLL per Element * 9");
    }
    
    // verbose
    std::string verbose(const std::string &stageKey) {
        using namespace bstring;
        std::stringstream ss;
        ss << boxTitle("FFT Solvers for " + stageKey);
        ss << boxEquals(0, 43, "real number precision", typeName<Real>());
        ss << boxEquals(0, 43, "\"howmany\" on GLL points", "{1, 3}");
        ss << boxEquals(0, 43, "\"howmany\" on elements", "{N3, N6, N9}");
        ss << internal::verbose("FFT on GLL points, \"howmany\" = 1",
                                gFFT_1.getPlannedNRs(), 1);
        ss << internal::verbose("FFT on GLL points, \"howmany\" = 3",
                                gFFT_3.getPlannedNRs(), 3);
        ss << internal::verbose("FFT on elements, \"howmany\" = N3",
                                gFFT_N3.getPlannedNRs(), nPEM * 3);
        ss << internal::verbose("FFT on elements, \"howmany\" = N6",
                                gFFT_N6.getPlannedNRs(), nPEM * 6);
        ss << internal::verbose("FFT on elements, \"howmany\" = N9",
                                gFFT_N9.getPlannedNRs(), nPEM * 9);
        ss << boxBaseline() << "\n\n";
        return ss.str();
    }
    
    
    //////////////// internal tools for verbose ////////////////
    namespace internal {
        // statistics for diagnosis and verbose
        void statistics(const std::vector<int> &plannedNRs, int HOWMANY,
                        int &numNR_MaxMPI, int &maxNR_MaxMPI,
                        int &sumNR_MaxMPI, double &memGB_MaxMPI) {
            // local
            int numNR = (int)plannedNRs.size();
            int maxNR = (numNR == 0) ? 0 : plannedNRs.back();
            int sumNR = vector_tools::sum(plannedNRs);
            double memGB = sizeof(Real) * maxNR * HOWMANY * 2. / 1e9;
            
            // mpi global
            numNR_MaxMPI = mpi::max(numNR);
            maxNR_MaxMPI = mpi::max(maxNR);
            sumNR_MaxMPI = mpi::max(sumNR);
            memGB_MaxMPI = mpi::max(memGB);
        }
        
        // diagnostic info after creating plans
        void diagnosticInfo(const std::vector<int> &plannedNRs, int HOWMANY) {
            // statistics
            int numNR_MaxMPI = 0;
            int maxNR_MaxMPI = 0;
            int sumNR_MaxMPI = 0;
            double memGB_MaxMPI = 0.;
            statistics(plannedNRs, HOWMANY,
                       numNR_MaxMPI, maxNR_MaxMPI, sumNR_MaxMPI, memGB_MaxMPI);
            
            // write diagnostic info
            using timer::gPreloopTimer;
            gPreloopTimer.message("HOWMANY = " + bstring::toString(HOWMANY));
            gPreloopTimer.message("# NR's (max over MPI ranks) = " +
                                  bstring::toString(numNR_MaxMPI));
            gPreloopTimer.message("max NR (max over MPI ranks) = " +
                                  bstring::toString(maxNR_MaxMPI));
            gPreloopTimer.message("sum NR (max over MPI ranks) = " +
                                  bstring::toString(sumNR_MaxMPI));
            gPreloopTimer.message("*** non-scalable memory allocation "
                                  "(max over MPI ranks) = " +
                                  bstring::toString(memGB_MaxMPI) + " GB ***");
        }
        
        // verbose
        std::string verbose(const std::string &subtitle,
                            const std::vector<int> &plannedNRs, int HOWMANY) {
            // statistics
            int numNR_MaxMPI = 0;
            int maxNR_MaxMPI = 0;
            int sumNR_MaxMPI = 0;
            double memGB_MaxMPI = 0.;
            statistics(plannedNRs, HOWMANY,
                       numNR_MaxMPI, maxNR_MaxMPI, sumNR_MaxMPI, memGB_MaxMPI);
            
            // verbose
            std::stringstream ss;
            ss << bstring::boxSubTitle(0, subtitle);
            ss << bstring::boxEquals(2, 41, "\"howmany\"", HOWMANY);
            ss << bstring::boxEquals(2, 41, "# logical sizes (NR) "
                                     "(max over MPI ranks)", numNR_MaxMPI);
            ss << bstring::boxEquals(2, 41, "maximum logical size "
                                     "(max over MPI ranks)", maxNR_MaxMPI);
            ss << bstring::boxEquals(2, 41, "sum of logical sizes "
                                     "(max over MPI ranks)", sumNR_MaxMPI);
            return ss.str();
        }
    }
}
