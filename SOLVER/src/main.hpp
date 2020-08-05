//
//  main.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/13/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  AxiSEM3D main

#ifndef main_hpp
#define main_hpp

//////////////////// includes ////////////////////
// environment
#include <iostream>
#include "mpi.hpp"
#include "io.hpp"
#include "inparam.hpp"
#include "timer.hpp"
#include "eigen_tools.hpp"

// mesh stage
#include "ExodusMesh.hpp"
#include "geodesy.hpp"
#include "ABC.hpp"
#include "NrField.hpp"
#include "LocalMesh.hpp"

// model stage
#include "Model3D.hpp"
#include "OceanLoad3D.hpp"
#include "SE_Model.hpp"

// domain stage
#include "AttBuilder.hpp"
#include "NewmarkTimeScheme.hpp"
#include "SymplecticTimeScheme.hpp"
#include "Domain.hpp"
#include "fft.hpp"

// source-receiver stage
#include "Source.hpp"
#include "StationOutput.hpp"
#include "ElementOutput.hpp"

// wavefield scanning
#include "WavefieldScanning.hpp"


//////////////////// top-level functions ////////////////////
// setup environment
void setupEnvironment(int argc, char *argv[]);

// build exodus mesh
std::unique_ptr<ExodusMesh> buildExodusMesh();

// ABC
std::unique_ptr<const ABC> buildABC(const ExodusMesh &exodusMesh);

// compute nr field and weights
eigen::DColX computeNrFieldAndWeights(ExodusMesh &exodusMesh);

// build nr-weighted local mesh
std::unique_ptr<LocalMesh>
buildLocalMesh(const ExodusMesh &exodusMesh, const eigen::DColX &weights,
               const std::string &weightsKey, const std::string &stageKey);

// build 3D models
void buildModels3D(const ExodusMesh &exodusMesh, const LocalMesh &localMesh,
                   std::vector<std::shared_ptr<const Model3D>> &models3D,
                   const std::string &stageKey, const ABC &abc);

// build SE model
std::unique_ptr<SE_Model>
buildSE_Model(const ExodusMesh &exodusMesh,
              const ABC &abc, LocalMesh &localMesh,
              const std::vector<std::shared_ptr<const Model3D>> &models3D,
              const std::string &stageKey);

// compute dt
double computeDt(const SE_Model &sem, const ABC &abc);

// attenuation
std::unique_ptr<const AttBuilder>
buildAttenuation(const ExodusMesh &exodusMesh, double dt);

// time scheme
std::unique_ptr<TimeScheme> buildTimeScheme();

// release to domain
void releaseToDomain(SE_Model &sem, const ABC &abc, LocalMesh &localMesh,
                     const std::unique_ptr<const AttBuilder> &attBuilder,
                     const TimeScheme &timeScheme, Domain &domain,
                     const std::string &stageKey);

// initalize FFT
void initalizeFFT(const std::string &stageKey);

// measure cost
eigen::DColX
measureCost(const SE_Model &sem, const ExodusMesh &exodusMesh,
            const LocalMesh &localMesh, const TimeScheme &timeScheme);

// release sources
void releaseSources(SE_Model &sem, Domain &domain, double dt,
                    TimeScheme &timeScheme);

// setup wavefield scanning
void setupWavefieldScanning(double dt, double period, int numTotalSteps,
                            const SE_Model &sem, Domain &domain);

#endif /* main_hpp */
