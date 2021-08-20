//
//  LocalMesh.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/17/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  local mesh: skeleton + element-GLL vicinity + mpi communication

#include "LocalMesh.hpp"
#include "ExodusMesh.hpp"
#include "metis.hpp"

// tools
#include "mpi.hpp"
#include "vector_tools.hpp"
#include "vicinity.hpp"

// plot and verbose
#include "NetCDF_Writer.hpp"
#include "io.hpp"
#include "bstring.hpp"
#include "timer.hpp"
#include "eigen_tools.hpp"

/////////////////////////// build //////////////////////////
// constructor
LocalMesh::
LocalMesh(const ExodusMesh &exodusMesh, const eigen::DColX &weights) {
    // domain decomposition
    eigen::IColX elemRank;
    mpi::enterSuper();
    if (mpi::super()) {
        timer::gPreloopTimer.begin("Mesh partitioning by Metis");
        decomposeExodusMesh(exodusMesh, weights, elemRank);
        timer::gPreloopTimer.ended("Mesh partitioning by Metis");
    }
    mpi::enterWorld();
    
    // build local skeleton
    mpi::enterInfer();
    timer::gPreloopTimer.begin("Building local skeleton");
    buildLocalSkeleton(exodusMesh, elemRank);
    timer::gPreloopTimer.ended("Building local skeleton");
    
    // build element-GLL vicinity and mpi communication
    timer::gPreloopTimer.begin("Building element-GLL vicinity "
                               "and mpi communication");
    buildElementGLL_CommMPI(exodusMesh, elemRank);
    timer::gPreloopTimer.ended("Building element-GLL vicinity "
                               "and mpi communication");
    mpi::enterWorld();
}

// domain decomposition
void LocalMesh::decomposeExodusMesh(const ExodusMesh &exodusMesh,
                                    const eigen::DColX &weights,
                                    eigen::IColX &elemRank) {
    // metis
    double obj = metis::decompose(exodusMesh.getConnectivity(),
                                  exodusMesh.getIsElementFluid(),
                                  weights, mpi::nprocWorld(),
                                  mpi::rankWorld(), elemRank);
    // choose highest quality
    std::vector<double> objAll;
    mpi::gather(obj, objAll, MPI_ALL);
    int minProc = (int)(std::min_element(objAll.begin(), objAll.end())
                        - objAll.begin());
    mpi::bcastEigen(elemRank, minProc);
    timer::gPreloopTimer.message
    (eigen_tools::memoryInfo(elemRank, "element mpi-rank (super-only)"));
}

// build local skeleton
void LocalMesh::buildLocalSkeleton(const ExodusMesh &exodusMesh,
                                   const eigen::IColX &elemRank) {
    // gather world rank
    timer::gPreloopTimer.begin("Gathering world rank in mpi-groups");
    std::vector<int> vecWorldRank;
    mpi::gather(mpi::rankWorld(), vecWorldRank, 0);
    timer::gPreloopTimer.ended("Gathering world rank in mpi-groups");
    
    // compute local skeleton on root (leader of inferior group)
    timer::gPreloopTimer.begin("Building local skeleton on mpi-group leaders");
    std::vector<eigen::IColX> vmL2GE;
    std::vector<eigen::IMatX4_RM> vmConn;
    std::vector<eigen::DMatX2_RM> vmCrds;
    std::vector<eigen::IColX> vmGeom;
    std::vector<eigen::IColX> vmIsFl;
    std::vector<eigen::IColX> vmNdNr;
    if (mpi::root()) {
        // reserve
        vmL2GE.reserve(mpi::nproc());
        vmConn.reserve(mpi::nproc());
        vmCrds.reserve(mpi::nproc());
        vmGeom.reserve(mpi::nproc());
        vmIsFl.reserve(mpi::nproc());
        vmNdNr.reserve(mpi::nproc());
        
        // loop over inferior ranks
        for (int rank = 0; rank < mpi::nproc(); rank++) {
            // element-rank to rank-element
            std::vector<int> vL2GE;
            for (int ieg = 0; ieg < elemRank.rows(); ieg++) {
                if (elemRank(ieg) == vecWorldRank[rank]) {
                    vL2GE.push_back(ieg);
                }
            }
            
            // check empty processor
            if (vL2GE.size() == 0) {
                throw
                std::runtime_error("LocalMesh::buildLocalSkeleton || "
                                   "One processor has no element. || "
                                   "Try decreasing the number of processors.");
            }
            
            // global element tag
            vmL2GE.push_back(Eigen::Map<eigen::IColX>(vL2GE.data(),
                                                      vL2GE.size()));
            
            // elemental fields simply sliced by vL2GE
            vmGeom.push_back(exodusMesh.getGeometryType()(vL2GE));
            vmIsFl.push_back(exodusMesh.getIsElementFluid()(vL2GE));
            
            // nodes: local to global
            const eigen::IMatX4_RM &matConnG =
            exodusMesh.getConnectivity()(vL2GE, Eigen::all);
            std::vector<int> nL2G(matConnG.data(),
                                  matConnG.data() + matConnG.size());
            vector_tools::sortUnique(nL2G);
            
            // nodal fields simply sliced by nL2G
            vmCrds.push_back(exodusMesh.getNodalCoords()(nL2G, Eigen::all));
            vmNdNr.push_back(exodusMesh.getNrAtNodes()(nL2G));
            
            // sliced connectivity with local nodes
            // form an inverse map from global to local for fast search
            std::map<int, int> nG2L;
            for (int inl = 0; inl < nL2G.size(); inl++) {
                nG2L.insert({nL2G[inl], inl});
            }
            // node tag from global to local
            eigen::IMatX4_RM mConn(vL2GE.size(), 4);
            for (int iel = 0; iel < vL2GE.size(); iel++) {
                for (int ivt = 0; ivt < 4; ivt++) {
                    mConn(iel, ivt) = nG2L.at(matConnG(iel, ivt));
                }
            }
            vmConn.push_back(mConn);
        }
    }
    timer::gPreloopTimer.ended("Building local skeleton on mpi-group leaders");
    
    // scatter
    timer::gPreloopTimer.begin("Scattering local skeleton in mpi-groups");
    mpi::scatterEigen(vmL2GE, mL2G_Element, 0);
    mpi::scatterEigen(vmConn, mConnectivity, 0);
    mpi::scatterEigen(vmCrds, mNodalCoords, 0);
    mpi::scatterEigen(vmGeom, mGeometryType, 0);
    mpi::scatterEigen(vmIsFl, mIsElementFluid, 0);
    mpi::scatterEigen(vmNdNr, mNodalNr, 0);
    timer::gPreloopTimer.ended("Scattering local skeleton in mpi-groups");
}

// build element-GLL vicinity
void LocalMesh::buildElementGLL_CommMPI(const ExodusMesh &exodusMesh,
                                        const eigen::IColX &elemRank) {
    // local GLL
    timer::gPreloopTimer.begin("Constructing local element-GLL vicinity");
    std::vector<eigen::IColX> localNeighbours; // useless
    int nLocalGLL = vicinity::
    connectivityToElementGLL(mConnectivity, mElementGLL, localNeighbours);
    timer::gPreloopTimer.ended("Constructing local element-GLL vicinity");
    
    // gather
    timer::gPreloopTimer.begin("Gathering input data in mpi-groups");
    std::vector<int> vecNLocalGLL;
    std::vector<eigen::IMatXN_RM> vecElementGLL;
    std::vector<eigen::IColX> vecL2G_Element;
    std::vector<int> vecWorldRank;
    mpi::gather(nLocalGLL, vecNLocalGLL, 0);
    mpi::gatherEigen(mElementGLL, vecElementGLL, 0);
    mpi::gatherEigen(mL2G_Element, vecL2G_Element, 0);
    mpi::gather(mpi::rankWorld(), vecWorldRank, 0);
    timer::gPreloopTimer.ended("Gathering input data in mpi-groups");
    
    // root
    timer::gPreloopTimer.begin("Building L2G-GLL and mpi communication "
                               "on mpi-group leaders");
    std::vector<eigen::IColX> vecL2G_GLL;
    std::vector<std::vector<int>> vecVecCommProcGLL;
    if (mpi::root()) {
        // global element-GLL vicinity
        timer::gPreloopTimer.begin("Constructing global element-GLL vicinity");
        eigen::IMatXN_RM globalElementGLL;
        std::vector<eigen::IColX> globalNeighbours;
        vicinity::connectivityToElementGLL(exodusMesh.getConnectivity(),
                                           globalElementGLL, globalNeighbours);
        // message large variables (here is a memory peak)
        timer::gPreloopTimer.message
        (eigen_tools::memoryInfo(globalElementGLL,
                                 "global element-GLL vicinity (super-only)"));
        int gNeighbSize = (int)vector_tools::totalSize(globalNeighbours);
        timer::gPreloopTimer.message
        (eigen_tools::memoryInfo(eigen::IColX(gNeighbSize),
                                 "global neighbourhood (super-only)"));
        timer::gPreloopTimer.ended("Constructing global element-GLL vicinity");
        
        // L2G GLL matching by (element, nPol^2)
        // L2G_Element is the bridge between local and global GLL tags
        timer::gPreloopTimer.begin("Matching L2G-GLL by (element, nPol^2)");
        vecL2G_GLL.reserve(mpi::nproc());
        // loop over inferior ranks
        for (int rank = 0; rank < mpi::nproc(); rank++) {
            const eigen::IColX &myL2G_GLL =
            vicinity::buildL2G_GLL(vecNLocalGLL[rank],
                                   vecL2G_Element[rank],
                                   vecElementGLL[rank],
                                   globalElementGLL);
            vecL2G_GLL.push_back(myL2G_GLL);
        }
        timer::gPreloopTimer.ended("Matching L2G-GLL by (element, nPol^2)");
        
        // mpi communication
        timer::gPreloopTimer.begin("Building mpi communication");
        vecVecCommProcGLL.reserve(mpi::nproc());
        // loop over inferior ranks
        for (int rank = 0; rank < mpi::nproc(); rank++) {
            const std::vector<int> &vecCommProcGLL =
            vicinity::buildCommMPI(vecL2G_Element[rank],
                                   globalNeighbours,
                                   elemRank, vecWorldRank[rank],
                                   exodusMesh.getConnectivity(),
                                   vecElementGLL[rank],
                                   vecL2G_GLL[rank]);
            vecVecCommProcGLL.push_back(vecCommProcGLL);
        }
        timer::gPreloopTimer.ended("Building mpi communication");
    }
    timer::gPreloopTimer.ended("Building L2G-GLL and mpi communication "
                               "on mpi-group leaders");
    
    // scatter
    timer::gPreloopTimer.begin("Scattering output data in mpi-groups");
    mpi::scatterEigen(vecL2G_GLL, mL2G_GLL, 0);
    std::vector<int> vecCommProcGLL;
    mpi::scatter(vecVecCommProcGLL, vecCommProcGLL, 0);
    
    // cast mpi communication
    for (auto it = vecCommProcGLL.begin(); it != vecCommProcGLL.end(); ++it) {
        mCommProc.push_back(*it);
        mCommMyGLL.push_back(std::vector<int>());
        for (++it; it != vecCommProcGLL.end(); ++it) {
            if ((*it) == -1) {
                // find separator, go to next rank
                break;
            }
            mCommMyGLL.back().push_back(*it);
        }
    }
    timer::gPreloopTimer.ended("Scattering output data in mpi-groups");
}


/////////////////////////// info //////////////////////////
// verbose
std::string LocalMesh::
verbose(const std::string &meshTitle, const std::string &weightsKey,
        const eigen::DColX &weights) const {
    // element distribution
    std::vector<int> evec;
    mpi::gather((int)mL2G_Element.rows(), evec, 0);
    
    // weights distribution
    std::vector<double> wvec;
    mpi::gather(weights(mL2G_Element).sum(), wvec, 0);
    
    // GLL points
    std::vector<int> gvec;
    mpi::gather((int)mL2G_GLL.size(), gvec, 0);
    
    // communication
    std::vector<int> crvec;
    std::vector<int> cgvec;
    mpi::gather((int)mCommProc.size(), crvec, 0);
    mpi::gather((int)vector_tools::totalSize(mCommMyGLL), cgvec, 0);
    
    // verbose
    std::stringstream ss;
    if (mpi::root()) {
        using namespace bstring;
        ss << boxTitle(meshTitle);
        ss << boxEquals(0, 24, "weights used", weightsKey);
        ss << boxEquals(0, 24, "# processors", mpi::nproc());
        
        // weights
        eigen::DColX wcol = Eigen::Map<eigen::DColX>(wvec.data(), wvec.size());
        wcol /= wcol.sum();
        double mean = wcol.mean();
        double sd = sqrt((wcol.array() - mean).square().sum() / wcol.rows());
        ss << boxSubTitle(0, "Weights distribution");
        ss << boxEquals(2, 22, "normalized total", 1);
        ss << boxEquals(2, 22, "min local distribution", wcol.minCoeff());
        ss << boxEquals(2, 22, "max local distribution", wcol.maxCoeff());
        ss << boxEquals(2, 24, "unbalance factor (σ/μ)", sd / mean);
        
        // element
        const eigen::DColX &ecol =
        Eigen::Map<const eigen::IColX>(evec.data(), evec.size()).cast<double>();
        mean = ecol.mean();
        sd = sqrt((ecol.array() - mean).square().sum() / ecol.rows());
        ss << boxSubTitle(0, "Element distribution");
        ss << boxEquals(2, 22, "total # elements", ecol.sum());
        ss << boxEquals(2, 22, "min local distribution", ecol.minCoeff());
        ss << boxEquals(2, 22, "max local distribution", ecol.maxCoeff());
        ss << boxEquals(2, 24, "unbalance factor (σ/μ)", sd / mean);
        
        // GLL
        const eigen::DColX &gcol =
        Eigen::Map<const eigen::IColX>(gvec.data(), gvec.size()).cast<double>();
        mean = gcol.mean();
        sd = sqrt((gcol.array() - mean).square().sum() / gcol.rows());
        ss << boxSubTitle(0, "GLL-point distribution");
        ss << boxEquals(2, 22, "total # GLL points", gcol.sum());
        ss << boxEquals(2, 22, "min local distribution", gcol.minCoeff());
        ss << boxEquals(2, 22, "max local distribution", gcol.maxCoeff());
        ss << boxEquals(2, 24, "unbalance factor (σ/μ)", sd / mean);
        
        // communication
        ss << boxSubTitle(0, "mpi communication");
        const eigen::IColX &crcol =
        Eigen::Map<const eigen::IColX>(crvec.data(), crvec.size());
        const eigen::IColX &cgcol =
        Eigen::Map<const eigen::IColX>(cgvec.data(), cgvec.size());
        ss << boxEquals(2, 22, "# rank-to-rank pairs", crcol.sum() / 2);
        ss << boxEquals(2, 22, "# GLL-to-GLL paris", cgcol.sum() / 2);
        ss << boxBaseline() << "\n\n";
    }
    return ss.str();
}

// plot domain decomposition
void LocalMesh::plotDD(const std::string &fname,
                       const eigen::DColX &weights) const {
    // compute coords at element center
    int nelem = (int)weights.rows();
    eigen::DMatX2_RM center = eigen::DMatX2_RM::Zero(nelem, 2);
    center(mL2G_Element, Eigen::all) =
    mpi::nodalToElemental(mConnectivity, mNodalCoords, false);
    mpi::sumEigen(center);
    
    // gather rank info
    eigen::IColX elemRank = eigen::IColX::Zero(nelem);
    elemRank(mL2G_Element).array() = mpi::rank();
    mpi::sumEigen(elemRank);
    
    // write to file
    if (mpi::root()) {
        // open
        NetCDF_Writer nw(io::gOutputDirectory + "/plots/" + fname, true);
        // coords
        nw.defineVariable("coords", {
            {"dim_elem", nelem}, {"dim_sz", 2}}, numerical::dErr);
        nw.writeWholeVariable("coords", center);
        // weights
        nw.defineVariable("weights", {
            {"dim_elem", nelem}}, numerical::dErr);
        nw.writeWholeVariable("weights", weights);
        // rank
        nw.defineVariable("mpi_rank", {
            {"dim_elem", nelem}}, (int)-1);
        nw.writeWholeVariable("mpi_rank", elemRank);
    }
}
