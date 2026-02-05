//
//  main.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/6/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  AxiSEM3D main

#include "main.hpp"

extern void fenv_setup();

int main(int argc, char *argv[]) {
    // setup floating-point environment
    fenv_setup();
    
    // initialize MPI
    mpi::initialize(&argc, &argv);
    
    // all other codes go into try
    try {
        // setup environment
        setupEnvironment(argc, argv);
        
        // start preloop timer
        timer::gPreloopTimer.begin("Pre-timeloop processing", '*');
        
        // build exodus mesh
        timer::gPreloopTimer.begin("Exodus mesh", '*');
        std::unique_ptr<ExodusMesh> exodusMesh = buildExodusMesh();
        timer::gPreloopTimer.ended("Exodus mesh", '*');
        
        // build ABC parameters
        timer::gPreloopTimer.begin("ABC parameters", '*');
        std::unique_ptr<const ABC> abc = buildABC(*exodusMesh);
        timer::gPreloopTimer.ended("ABC parameters", '*');
        
        // compute nr field
        timer::gPreloopTimer.begin("Nr(s,z) and weights", '*');
        eigen::DColX nrWeights = computeNrFieldAndWeights(*exodusMesh);
        timer::gPreloopTimer.ended("Nr(s,z) and weights", '*');
        
        // build nr-weighted local mesh
        timer::gPreloopTimer.begin("Nr-weighted local mesh", '*');
        std::unique_ptr<LocalMesh> localMesh =
        buildLocalMesh(*exodusMesh, nrWeights, "Nr(s,z)", "Stage-I");
        nrWeights.resize(0); // free memory
        timer::gPreloopTimer.ended("Nr-weighted local mesh", '*');
        
        // build 3D models
        timer::gPreloopTimer.begin("3D models", '*');
        std::vector<std::shared_ptr<const Model3D>> models3D;
        buildModels3D(*exodusMesh, *localMesh, models3D, "Stage-I", *abc);
        timer::gPreloopTimer.ended("3D models", '*');
        
        // build SE model
        timer::gPreloopTimer.begin("Spectral element model (ROOT)", '*');
        std::unique_ptr<SE_Model> sem =
        buildSE_Model(*exodusMesh, *abc, *localMesh, models3D, "Stage-I");
        timer::gPreloopTimer.ended("Spectral element model (ROOT)", '*');
        
        // dt
        timer::gPreloopTimer.begin("Time step", '*');
        double dt = computeDt(*sem, *abc);
        timer::gPreloopTimer.ended("Time step", '*');
        
        // attenuation
        timer::gPreloopTimer.begin("Attenuation builder", '*');
        std::unique_ptr<const AttBuilder> attBuilder =
        buildAttenuation(*exodusMesh, dt);
        timer::gPreloopTimer.ended("Attenuation builder", '*');
        
        // time scheme
        timer::gPreloopTimer.begin("Time Scheme", '*');
        std::unique_ptr<TimeScheme> timeScheme = buildTimeScheme();
        timer::gPreloopTimer.ended("Time Scheme", '*');
        
        // domain
        timer::gPreloopTimer.begin("Nr-weighted domain", '*');
        std::shared_ptr<Domain> domain = Domain::createEmpty();
        releaseToDomain(*sem, *abc, *localMesh, attBuilder, *timeScheme,
                        *domain, "Stage-I");
        timer::gPreloopTimer.ended("Nr-weighted domain", '*');
        
        // FFT
        timer::gPreloopTimer.begin("Preloop FFT", '*');
        initalizeFFT("Preloop");
        timer::gPreloopTimer.ended("Preloop FFT", '*');
        
        // Stage-I: nr-weighted domain
        /////////////////////////////////////
        // Stage-II: cost-weighted domain
        
        // measure cost
        timer::gPreloopTimer.begin("Cost measurement", '*');
        eigen::DColX costWeights =
        measureCost(*sem, *exodusMesh, *localMesh, *timeScheme);
        timer::gPreloopTimer.ended("Cost measurement", '*');
        
        // free memory
        timer::gPreloopTimer.begin("Freeing memory Stage-I", '*');
        localMesh.reset();
        sem.reset();
        domain.reset();
        timer::gPreloopTimer.ended("Freeing memory Stage-I", '*');
        
        // build cost-weighted local mesh
        timer::gPreloopTimer.begin("Cost-weighted local mesh", '*');
        localMesh = buildLocalMesh(*exodusMesh, costWeights,
                                   "measured computational cost", "Stage-II");
        costWeights.resize(0); // free memory
        timer::gPreloopTimer.ended("Cost-weighted local mesh", '*');
        
        // rebuild MPI-dependent 3D models
        timer::gPreloopTimer.begin("3D models (rebuild)", '*');
        buildModels3D(*exodusMesh, *localMesh, models3D, "Stage-II", *abc);
        timer::gPreloopTimer.ended("3D models (rebuild)", '*');
        
        // build SE model
        timer::gPreloopTimer.begin("Spectral element model (ROOT)", '*');
        sem =
        buildSE_Model(*exodusMesh, *abc, *localMesh, models3D, "Stage-II");
        timer::gPreloopTimer.ended("Spectral element model (ROOT)", '*');
        
        // domain
        timer::gPreloopTimer.begin("Cost-weighted domain", '*');
        domain = Domain::createEmpty();
        releaseToDomain(*sem, *abc, *localMesh, attBuilder, *timeScheme,
                        *domain, "Stage-II");
        timer::gPreloopTimer.ended("Cost-weighted domain", '*');
        
        // free memory
        timer::gPreloopTimer.begin("Freeing memory Stage-II", '*');
        // needed for wavefield scanning
        double meshPeriod = exodusMesh->getGlobalVariable("minimum_period");
        // needed for element output
        double distTol = exodusMesh->getGlobalVariable("dist_tolerance");
        exodusMesh.reset();
        localMesh.reset();
        abc.reset();
        attBuilder.reset();
        timer::gPreloopTimer.ended("Freeing memory Stage-II", '*');
        
        // Stage-II: cost-weighted domain
        /////////////////////////////////////
        // Stage-III: source and output
        
        // source
        timer::gPreloopTimer.begin("Source", '*');
        releaseSources(*sem, *domain, dt, *timeScheme);
        timer::gPreloopTimer.ended("Source", '*');
        
        // station output
        timer::gPreloopTimer.begin("Station groups", '*');
        StationOutput::release(*sem, *domain, dt,
                               timeScheme->getT0(), timeScheme->getT1(),
                               timeScheme->getNumTimeSteps());
        timer::gPreloopTimer.ended("Station groups", '*');
        
        // element output
        timer::gPreloopTimer.begin("Element groups", '*');
        ElementOutput::release(*sem, *domain, dt,
                               timeScheme->getT0(), timeScheme->getT1(),
                               timeScheme->getNumTimeSteps(), distTol);
        timer::gPreloopTimer.ended("Element groups", '*');
        
        // wavefield scanning
        timer::gPreloopTimer.begin("Wavefield scanning", '*');
        setupWavefieldScanning(dt, meshPeriod, timeScheme->getNumTimeSteps(),
                               *sem, *domain);
        timer::gPreloopTimer.ended("Wavefield scanning", '*');
        
        // free memory
        timer::gPreloopTimer.begin("Freeing memory Stage-III", '*');
        sem.reset();
        timer::gPreloopTimer.ended("Freeing memory Stage-III", '*');
        
        // wavefield output: initialize
        timer::gPreloopTimer.begin("Output files", '*');
        domain->initializeOutput();
        timer::gPreloopTimer.ended("Output files", '*');
        
        // FFT (must be initialized after adding sources and receivers)
        timer::gPreloopTimer.begin("Time loop FFT", '*');
        initalizeFFT("Time Loop");
        timer::gPreloopTimer.ended("Time loop FFT", '*');
        
        // end preloop timer
        timer::gPreloopTimer.ended("Pre-timeloop processing", '*');
        
        // Preloop
        /////////////////////////////////////
        // Time loop
        
        // hand over domain to time scheme
        timeScheme->setDomain(domain);
        
        // go go go
        timeScheme->solve();
        
    } catch (const std::exception &e) {
        std::cerr << bstring::exception(e);
        std::cerr.flush();
        mpi::abort(1);
    }
    
    // finalize mpi
    mpi::finalize();
    return 0;
}


//////////////////// top-level functions ////////////////////
// setup environment
void setupEnvironment(int argc, char *argv[]) {
    // initialize IO dirs from CLI first
    io::parseIOArguments(argc, argv);

    // initialize IO
    std::string warningIO;
    io::verifyDirectories(warningIO);
    
    // inparam setup
    inparam::setup();
    
    // io verbose
    io::setupVerbose();
    
    // mpi group
    mpi::setupGroup(inparam::gInparamAdvanced.get<int>
                    ("mpi:nproc_per_group"));
    
    // welcome
    if (io::gVerbose != io::VerboseLevel::None) {
        io::cout << io::welcome();
    }
    
    // verbose mpi, io, inparam
    if (io::gVerbose == io::VerboseLevel::Detailed) {
        io::cout << mpi::verbose();
        io::cout << io::verbose();
        io::cout << inparam::verbose();
    }
    
    // io warning
    if (io::gVerboseWarnings) {
        io::cout << warningIO;
    }
    
    // preloop timer
    timer::setupPreloopDiagnosis(inparam::gInparamAdvanced.get<bool>
                                 ("develop:diagnose_preloop"));
    if (io::gVerbose == io::VerboseLevel::Detailed) {
        io::cout << timer::verbose();
    }
}

// build exodus mesh
std::unique_ptr<ExodusMesh> buildExodusMesh() {
    // exodus mesh
    timer::gPreloopTimer.begin("Building Exodus mesh");
    std::unique_ptr<ExodusMesh> exodusMesh = std::make_unique<ExodusMesh>
    (inparam::gInparamModel.get<std::string>("model1D:exodus_mesh"));
    if (io::gVerbose != io::VerboseLevel::None) {
        io::cout << exodusMesh->verbose();
    }
    timer::gPreloopTimer.ended("Building Exodus mesh");
    
    // geodesy
    timer::gPreloopTimer.begin("Setting up geodesy");
    geodesy::setup(*exodusMesh);
    if (io::gVerbose != io::VerboseLevel::None) {
        io::cout << geodesy::verbose(exodusMesh->getDiscontinuities());
    }
    timer::gPreloopTimer.ended("Setting up geodesy");
    return exodusMesh;
}

// ABC
std::unique_ptr<const ABC> buildABC(const ExodusMesh &exodusMesh) {
    std::unique_ptr<ABC> abc = ABC::buildInparam(exodusMesh);
    if (io::gVerbose != io::VerboseLevel::None) {
        io::cout << abc->verbose();
    }
    return abc;
}

// compute nr field
eigen::DColX computeNrFieldAndWeights(ExodusMesh &exodusMesh) {
    mpi::enterSuper();
    if (mpi::super()) {
        // create nr field on super (timer inside)
        std::unique_ptr<const NrField> nrField =
        NrField::buildInparam(exodusMesh);
        if (io::gVerbose != io::VerboseLevel::None) {
            io::cout << nrField->verbose();
        }
        
        // compute nr on mesh
        timer::gPreloopTimer.begin("Computing Nr at mesh nodes");
        // additional parameters
        bool boundByInplane =
        inparam::gInparamNr.get<bool>("bound_Nr_by_inplane");
        bool useLuckyNumbers =
        inparam::gInparamAdvanced.get<bool>("develop:fftw_lucky_numbers");
        // compute
        exodusMesh.formNrAtNodes(*nrField, boundByInplane, useLuckyNumbers);
        if (io::gVerbose != io::VerboseLevel::None) {
            io::cout << exodusMesh.verboseNr(boundByInplane, useLuckyNumbers);
        }
        timer::gPreloopTimer.ended("Computing Nr at mesh nodes");
    }
    mpi::enterWorld();
    
    timer::gPreloopTimer.begin("Computing Nr-weights");
    eigen::DColX weights;
    mpi::enterSuper();
    if (mpi::super()) {
        weights = mpi::
        nodalToElemental(exodusMesh.getConnectivity(),
                         exodusMesh.getNrAtNodes().cast<double>().eval(), true);
        // also differentiate solid (x10) and fluid (x1)
        weights.array() *=
        (10 - 9 * exodusMesh.getIsElementFluid().array()).cast<double>();
    }
    mpi::enterWorld();
    mpi::bcastEigen(weights);
    timer::gPreloopTimer.ended("Computing Nr-weights");
    return weights;
}

// build nr-weighted local mesh
std::unique_ptr<LocalMesh>
buildLocalMesh(const ExodusMesh &exodusMesh, const eigen::DColX &weights,
               const std::string &weightsKey, const std::string &stageKey) {
    // build
    timer::gPreloopTimer.begin("Building local mesh");
    std::unique_ptr<LocalMesh> localMesh =
    std::make_unique<LocalMesh>(exodusMesh, weights);
    timer::gPreloopTimer.ended("Building local mesh");
    
    // verbose
    if (io::gVerbose == io::VerboseLevel::Detailed) {
        io::cout << localMesh->verbose("Mesh Partition " + stageKey,
                                       weightsKey, weights);
    }
    
    // plot
    bool plotDD = inparam::gInparamAdvanced.get<bool>
    ("mpi:plot_domain_decomposition");
    if (plotDD) {
        timer::gPreloopTimer.begin("Plotting domain decompositon");
        localMesh->plotDD("domain_decomposition_" + stageKey + ".nc", weights);
        timer::gPreloopTimer.ended("Plotting domain decompositon");
    }
    
    // return
    return localMesh;
}

// build 3D models
void buildModels3D(const ExodusMesh &exodusMesh, const LocalMesh &localMesh,
                   std::vector<std::shared_ptr<const Model3D>> &models3D,
                   const std::string &stageKey, const ABC &abc) {
    // build
    timer::gPreloopTimer.begin("Building 3D models");
    bool rebuilding = !(stageKey == "Stage-I");
    Model3D::buildInparam(exodusMesh, localMesh, models3D, rebuilding);
    timer::gPreloopTimer.ended("Building 3D models");
    
    // verbose
    if (io::gVerbose != io::VerboseLevel::None && !rebuilding) {
        timer::gPreloopTimer.begin("Verbosing 3D models");
        std::stringstream ss;
        ss << bstring::boxTitle("3D Models");
        for (const std::shared_ptr<const Model3D> &m: models3D) {
            ss << m->verbose();
        }
        if (models3D.size() == 0) {
            ss << "* No 3D models in this simulation.\n";
        }
        ss << bstring::boxBaseline() << "\n\n";
        io::cout << ss.str();
        timer::gPreloopTimer.ended("Verbosing 3D models");
    }
    
    // check surface: ABC vs Ocean Load
    if (!io::gVerboseWarnings || rebuilding){
        return;
    }
    bool oceanload = false;
    for (const std::shared_ptr<const Model3D> &m: models3D) {
        if (std::dynamic_pointer_cast<const OceanLoad3D>(m)) {
            oceanload = true;
            break;
        }
    }
    const std::vector<std::string> &keys = abc.getBoundaryKeys();
    bool topABC = (std::find(keys.begin(), keys.end(), "TOP") != keys.end());
    if (oceanload && topABC) {
        io::cout << bstring::warning
        ("main::buildModels3D || "
         "The mesh surface is both specified as an absorbing boundary || "
         "and applied with ocean load. Though allowed by AxiSEM3D, || "
         "such a mixed boundary condition seems strange.");
    }
}

// build SE model
std::unique_ptr<SE_Model>
buildSE_Model(const ExodusMesh &exodusMesh,
              const ABC &abc, LocalMesh &localMesh,
              const std::vector<std::shared_ptr<const Model3D>> &models3D,
              const std::string &stageKey) {
    // build model
    bool useLuckyNumbers =
    inparam::gInparamAdvanced.get<bool>("develop:fftw_lucky_numbers");
    std::unique_ptr<SE_Model> sem =
    std::make_unique<SE_Model>(exodusMesh, abc, localMesh, models3D,
                               useLuckyNumbers);
    // free dummy memory
    localMesh.freeMemorySE_ModelBuilt();
    
    // verbose
    if (io::gVerbose == io::VerboseLevel::Detailed) {
        io::cout << sem->verbose("SE Model (" + stageKey + ")");
    }
    return sem;
}

// compute dt
double computeDt(const SE_Model &sem, const ABC &abc) {
    // get courant
    double courant =
    inparam::gInparamSource.getWithBounds("time_axis:Courant_number", 0.0);
    // automatically determined by mesh
    eigen::DCol2 sz;
    double dtMesh = sem.computeDt(courant, abc, sz);
    // enforced
    double dtEnforce = inparam::gInparamSource.
    getWithOptions<double>("time_axis:enforced_dt", {{"NONE", -1.}});
    // finally used
    bool useMesh = dtEnforce < numerical::dEpsilon;
    double dtUse = useMesh ? dtMesh : dtEnforce;
    
    // verbose
    if (io::gVerbose != io::VerboseLevel::None) {
        std::stringstream ss;
        ss << bstring::boxTitle("Time Step");
        ss << bstring::boxEquals(0, 22, "Δt determined by mesh", dtMesh);
        ss << bstring::boxEquals(0, 21, "   Courant number", courant);
        if (geodesy::isCartesian()) {
            ss << bstring::boxEquals(0, 21, "   location (s,z)",
                                     bstring::range(sz(0), sz(1)));
        } else {
            const auto &rt = geodesy::sz2rtheta(sz, false);
            ss << bstring::boxEquals(0, 22, "   location (r,θ)",
                                     bstring::range(rt(0), rt(1)));
        }
        if (useMesh) {
            ss << bstring::boxEquals(0, 22, "Δt enforced by user", "NONE");
        } else {
            ss << bstring::boxEquals(0, 22, "Δt enforced by user", dtEnforce);
        }
        ss << bstring::boxEquals(0, 22, "Δt to be used", dtUse);
        ss << bstring::boxBaseline() << "\n\n";
        io::cout << ss.str();
    }
    // warning
    if (io::gVerboseWarnings && dtEnforce > dtMesh) {
        io::cout <<
        bstring::warning("main::computeDt || "
                         "Enforcing a Δt greater than mesh-determined || "
                         "can cause numerical instability.");
    }
    // diagnose
    timer::gPreloopTimer.message("mesh Δt = " + bstring::toString(dtMesh));
    timer::gPreloopTimer.message("user Δt = " + bstring::toString(dtEnforce));
    timer::gPreloopTimer.message("used Δt = " + bstring::toString(dtUse));
    return dtUse;
}

// attenuation
std::unique_ptr<const AttBuilder>
buildAttenuation(const ExodusMesh &exodusMesh, double dt) {
    // options
    int cg4 = inparam::gInparamModel.getWithLimits<int>
    ("attenuation", {{"NONE", -1}, {"FULL", 0}, {"CG4", 1}});
    
    // create
    std::unique_ptr<const AttBuilder> attBuilder = nullptr;
    if (cg4 >= 0 && exodusMesh.hasAttenuation()) {
        attBuilder = std::make_unique<AttBuilder>(exodusMesh, cg4, dt);
    }
    
    // verbose
    if (io::gVerbose != io::VerboseLevel::None) {
        io::cout << AttBuilder::verbose(attBuilder);
    }
    return attBuilder;
}

// time scheme
std::unique_ptr<TimeScheme> buildTimeScheme() {
    // key
    const std::string &tsKey = inparam::gInparamSource.
    getWithLimits<std::string>("time_axis:integrator", {
        {"NEWMARK", "NEWMARK"},
        {"SYMPLECTIC", "SYMPLECTIC"}});
    // intervals
    int vint =
    inparam::gInparamAdvanced.getWithBounds("verbose:loop_info_interval", 1);
    int sint =
    inparam::gInparamAdvanced.getWithBounds("verbose:stability_interval", 1);
    
    // create time scheme
    if (tsKey == "NEWMARK") {
        return std::make_unique<NewmarkTimeScheme>(vint, sint);
    } else if (tsKey == "SYMPLECTIC") {
        return std::make_unique<SymplecticTimeScheme>(vint, sint);
    }
    // never reach here
    return nullptr;
}

// release to domain
void releaseToDomain(SE_Model &sem, const ABC &abc, LocalMesh &localMesh,
                     const std::unique_ptr<const AttBuilder> &attBuilder,
                     const TimeScheme &timeScheme, Domain &domain,
                     const std::string &stageKey) {
    // domain
    sem.release(abc, localMesh, attBuilder, timeScheme, domain);
    // free dummy memory
    localMesh.freeMemorySE_ModelReleased();
    
    // verbose domain
    bool stageI = (stageKey == "Stage-I");
    if (io::gVerbose == io::VerboseLevel::Detailed) {
        io::cout << domain.verbose("Computational Domain (" + stageKey + ")");
    } else if (io::gVerbose == io::VerboseLevel::Essential && !stageI) {
        io::cout << domain.verbose("Computational Domain");
    }
    
    // initialize gradient, just do once
    if (stageI) {
        GradientQuadrature<numerical::Real>::
        setGMat(spectrals::gGMatrixGLL, spectrals::gGMatrixGLJ);
    }
}

// initalize FFT
void initalizeFFT(const std::string &stageKey) {
    // initialize FFT
    double timeLimit = inparam::gInparamAdvanced.
    getWithBounds("develop:time_limit_for_fftw_planning", 0.);
    fft::createPlans(timeLimit);
    
    // verbose FFT
    if (io::gVerbose == io::VerboseLevel::Detailed) {
        io::cout << fft::verbose(stageKey);
    }
}

// measure cost
eigen::DColX
measureCost(const SE_Model &sem, const ExodusMesh &exodusMesh,
            const LocalMesh &localMesh, const TimeScheme &timeScheme) {
    bool measurePoint = inparam::gInparamAdvanced.getWithLimits<bool>
    ("mpi::weight_for_load_balancing", {
        {"ELEMENT", false}, {"ELEMENT_POINT", true}});
    return sem.measureCost(exodusMesh, localMesh, timeScheme, measurePoint);
}

// release sources
void releaseSources(SE_Model &sem, Domain &domain, double dt,
                    TimeScheme &timeScheme) {
    // form R-Tree
    timer::gPreloopTimer.begin("Forming R-tree on mesh");
    sem.formInplaneRTree();
    timer::gPreloopTimer.ended("Forming R-tree on mesh");
    
    // release sources
    timer::gPreloopTimer.begin("Releasing sources");
    double t0 = 0.;
    Source::release(sem, domain, dt, t0);
    timer::gPreloopTimer.ended("Releasing sources");
    
    // set simulation time
    timer::gPreloopTimer.begin("Setting simulation time");
    // read length and max steps
    double t1 =
    inparam::gInparamSource.getWithBounds("time_axis:record_length", 0.);
    int nStepMax =
    inparam::gInparamAdvanced.get<int>("develop:max_num_time_steps");
    if (nStepMax <= 0) {
        nStepMax = std::numeric_limits<int>::max();
    }
    
    // set up time scheme
    int step0 = (int)floor(t0 / dt);
    int step1 = (int)ceil(t1 / dt);
    int nStep = step1 - step0 + 1;
    timeScheme.setTime(step0 * dt, dt, std::min(nStep, nStepMax));
    
    // verbose time
    using namespace bstring;
    if (io::gVerbose != io::VerboseLevel::None) {
        std::stringstream ss;
        ss << boxTitle("Simulation Time");
        ss << boxEquals(0, 25, "Δt of simulation", dt);
        ss << boxEquals(0, 24, "t0 determined by sources", t0);
        ss << boxEquals(0, 24, "t1 specified by user", t1);
        ss << boxEquals(0, 24, "[t0, t1] used by solver",
                        range(step0 * dt, step1 * dt));
        ss << boxEquals(0, 24, "# time steps", nStep);
        ss << boxBaseline() << "\n\n";
        io::cout << ss.str();
    }
    timer::gPreloopTimer.message("Δt = " + toString(dt));
    timer::gPreloopTimer.message("t0 = " + toString(step0 * dt));
    timer::gPreloopTimer.message("t1 = " + toString(step1 * dt));
    timer::gPreloopTimer.message("nstep = " + toString(nStep));
    timer::gPreloopTimer.ended("Setting simulation time");
}

// setup wavefield scanning
void setupWavefieldScanning(double dt, double period, int numTotalSteps,
                            const SE_Model &sem, Domain &domain) {
    bool enableScanning =
    inparam::gInparamNr.get<bool>("wavefield_scanning:enable_scanning");
    if (enableScanning) {
        WavefieldScanning::setup(dt, period, numTotalSteps, sem, domain);
    }
}
