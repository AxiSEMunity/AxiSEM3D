//
//  SE_Model.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/26/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  spectral element model
//  Exodus mesh ---
//                | -> mesh --------
//  local mesh ----                | -> SE model (this) -> domain
//                     3D models ---


#include "SE_Model.hpp"
// constructor
#include "ExodusMesh.hpp"
#include "ABC.hpp"
#include "LocalMesh.hpp"
#include "Model3D.hpp"
#include "inparam.hpp"
// components
#include "GLLPoint.hpp"
#include "Quad.hpp"
#include "mpi.hpp"
// utils
#include "timer.hpp"
// release
#include "Messaging.hpp"
#include "Domain.hpp"
// measure
#include "Element.hpp"
#include "SolidPoint.hpp"
#include "FluidPoint.hpp"
#include "point_time.hpp"
#include "io.hpp"
#include <fstream>
// source-receiver
#include "Geometric3D.hpp"
// wavefield scanning
#include "vicinity.hpp"

// step 1: constructor
SE_Model::
SE_Model(const ExodusMesh &exodusMesh,
         const ABC &abc, const LocalMesh &localMesh,
         const std::vector<std::shared_ptr<const Model3D>> &models3D,
         bool useLuckyNumbers):
mDistToleranceMesh(exodusMesh.getGlobalVariable("dist_tolerance")) {
    ///////////////// generate empty GLL points /////////////////
    timer::gPreloopTimer.begin("Generating empty points");
    int ngll = (int)(localMesh.mL2G_GLL.rows());
    mGLLPoints.reserve(ngll);
    for (int igll = 0; igll < ngll; igll++) {
        mGLLPoints.push_back(GLLPoint());
    }
    timer::gPreloopTimer.ended("Generating empty points");
    
    ///////////////// generate Quads /////////////////
    timer::gPreloopTimer.begin("Generating Quads");
    int nquad = (int)(localMesh.mL2G_Element.rows());
    mQuads.reserve(nquad);
    for (int iquad = 0; iquad < nquad; iquad++) {
        mQuads.push_back(Quad(exodusMesh, localMesh, iquad, useLuckyNumbers));
    }
    timer::gPreloopTimer.ended("Generating Quads");
    
    ///////////////// 3D models /////////////////
    timer::gPreloopTimer.begin("Applying 3D models to Quads");
    for (const std::shared_ptr<const Model3D> &m: models3D) {
        timer::gPreloopTimer.begin("Applying 3D model:" + m->getModelName());
        m->applyTo(mQuads);
        timer::gPreloopTimer.ended("Applying 3D model:" + m->getModelName());
    }
    timer::gPreloopTimer.ended("Applying 3D models to Quads");
    
    // finished 3D properties
    timer::gPreloopTimer.begin("Finishing applying 3D models to Quads");
    for (Quad &quad: mQuads) {
        quad.finishing3D();
    }
    Undulation::finished3D();
    for (Quad &quad: mQuads) {
        quad.finished3D();
    }
    timer::gPreloopTimer.ended("Finishing applying 3D models to Quads");
    
    ///////////////// setup GLL points by Quads /////////////////
    timer::gPreloopTimer.begin("Setting up GLL points");
    // setup GLL by quads
    for (const Quad &quad: mQuads) {
        quad.setupGLL(abc, localMesh, mGLLPoints);
    }
    timer::gPreloopTimer.ended("Setting up GLL points");
    
    ///////////////// assemble mpi /////////////////
    timer::gPreloopTimer.begin("Assembling mass and normal");
    // init mpi buffer
    timer::gPreloopTimer.begin("Creating MPI buffers");
    std::vector<eigen::DColX> bufSend, bufRecv;
    int nCommProc = (int)localMesh.mCommProc.size();
    for (int iproc = 0; iproc < nCommProc; iproc++) {
        int commSize = 0;
        for (int igll: localMesh.mCommMyGLL[iproc]) {
            commSize += mGLLPoints[igll].getCommSize();
        }
        bufSend.push_back(eigen::DColX::Zero(commSize));
        bufRecv.push_back(eigen::DColX::Zero(commSize));
    }
    timer::gPreloopTimer.ended("Creating MPI buffers");
    
    // feed
    timer::gPreloopTimer.begin("Feeding MPI buffers");
    for (int iproc = 0; iproc < nCommProc; iproc++) {
        int row = 0;
        for (int igll: localMesh.mCommMyGLL[iproc]) {
            mGLLPoints[igll].feedComm(bufSend[iproc], row);
        }
    }
    timer::gPreloopTimer.ended("Feeding MPI buffers");
    
    // send and recv
    timer::gPreloopTimer.begin("Exchanging MPI buffers");
    std::vector<MPI_Request> reqSend(nCommProc, MPI_REQUEST_NULL);
    std::vector<MPI_Request> reqRecv(nCommProc, MPI_REQUEST_NULL);
    for (int iproc = 0; iproc < localMesh.mCommProc.size(); iproc++) {
        int commProc = localMesh.mCommProc[iproc];
        mpi::isend(commProc, bufSend[iproc], reqSend[iproc]);
        mpi::irecv(commProc, bufRecv[iproc], reqRecv[iproc]);
    }
    mpi::waitAll(reqRecv);
    timer::gPreloopTimer.ended("Exchanging MPI buffers");
    
    // extract
    timer::gPreloopTimer.begin("Extracting MPI buffers");
    for (int iproc = 0; iproc < nCommProc; iproc++) {
        int row = 0;
        for (int igll: localMesh.mCommMyGLL[iproc]) {
            mGLLPoints[igll].extractComm(bufRecv[iproc], row);
        }
    }
    mpi::waitAll(reqSend);
    timer::gPreloopTimer.ended("Extracting MPI buffers");
    timer::gPreloopTimer.ended("Assembling mass and normal");
    
    ///////////////// source-receiver /////////////////
    // geometric model
    timer::gPreloopTimer.begin("Casting geometric models");
    for (const std::shared_ptr<const Model3D> &m: models3D) {
        std::shared_ptr<const Geometric3D> g3D =
        std::dynamic_pointer_cast<const Geometric3D>(m);
        if (g3D) {
            timer::gPreloopTimer.message("Adding geometric3D model: " +
                                         g3D->getModelName());
            mModelsG3D.push_back(g3D);
        }
    }
    timer::gPreloopTimer.ended("Casting geometric models");
    
    // range
    timer::gPreloopTimer.begin("Computing local mesh range");
    const eigen::DMatX2_RM &sz = localMesh.mNodalCoords;
    mRangeS(0) = sz.col(0).minCoeff() - mDistToleranceMesh;
    mRangeS(1) = sz.col(0).maxCoeff() + mDistToleranceMesh;
    mRangeZ(0) = sz.col(1).minCoeff() - mDistToleranceMesh;
    mRangeZ(1) = sz.col(1).maxCoeff() + mDistToleranceMesh;
    const eigen::DMatX2_RM &rt = geodesy::sz2rtheta(sz, true);
    mRangeR(0) = rt.col(0).minCoeff() - mDistToleranceMesh;
    mRangeR(1) = rt.col(0).maxCoeff() + mDistToleranceMesh;
    mRangeT(0) = rt.col(1).minCoeff();
    mRangeT(1) = rt.col(1).maxCoeff();
    timer::gPreloopTimer.ended("Computing local mesh range");
}

// step 2: get dt for attenuation
double SE_Model::
computeDt(double courant, const ABC &abc, eigen::DCol2 &sz) const {
    // minimum dt over quads
    double dtMin = std::numeric_limits<double>::max();
    for (const Quad &quad: mQuads) {
        double dt = quad.computeDt(courant, abc);
        if (dtMin > dt) {
            dtMin = dt;
            sz = quad.getNodalSZ().rowwise().mean();
        }
    }
    
    // minimum dt over ranks
    std::vector<double> dtV;
    std::vector<eigen::DCol2> szV;
    mpi::gather(dtMin, dtV, MPI_ALL);
    mpi::gatherEigen(sz, szV, MPI_ALL);
    
    // rank with minimum dt
    Eigen::Index rankDt = -1;
    dtMin = Eigen::Map<eigen::DColX>(dtV.data(), dtV.size()).minCoeff(&rankDt);
    sz = szV[rankDt];
    
    // round to 4 significant figures to avoid dumping error
    std::stringstream ss;
    ss << std::scientific << std::setprecision(3) << dtMin;
    return bstring::cast<double>(ss.str(), "SE_Model::getDeltaT");
}

// step 3: release to domain
void SE_Model::release(const ABC &abc, const LocalMesh &localMesh,
                       const std::unique_ptr<const AttBuilder> &attBuilder,
                       const TimeScheme &timeScheme, Domain &domain) {
    // points
    timer::gPreloopTimer.begin("Releasing GLL points");
    for (GLLPoint &point: mGLLPoints) {
        point.release(abc, timeScheme, domain);
    }
    timer::gPreloopTimer.ended("Releasing GLL points");
    
    // elements
    timer::gPreloopTimer.begin("Releasing elements");
    for (Quad &quad: mQuads) {
        quad.release(localMesh, mGLLPoints, attBuilder, domain);
    }
    timer::gPreloopTimer.ended("Releasing elements");
    
    // mpi
    timer::gPreloopTimer.begin("Releasing messaging");
    std::unique_ptr<Messaging> msg = std::make_unique<Messaging>();
    for (int icom = 0; icom < localMesh.mCommProc.size(); icom++) {
        // create and add message on a single rank
        int rankOther = localMesh.mCommProc[icom];
        std::vector<MessageRank::MeshPoint> meshPoints;
        for (int igll: localMesh.mCommMyGLL[icom]) {
            meshPoints.push_back({
                mGLLPoints[igll].getGlobalTag(),
                mGLLPoints[igll].getSolidPoint(),
                mGLLPoints[igll].getFluidPoint()});
        }
        std::unique_ptr<MessageRank> msgRank =
        std::make_unique<MessageRank>(rankOther, meshPoints);
        msg->addRankComm(msgRank);
    }
    domain.setMessaging(msg);
    timer::gPreloopTimer.ended("Releasing messaging");
}

// step 4: measure
eigen::DColX SE_Model::
measureCost(const ExodusMesh &exodusMesh, const LocalMesh &localMesh,
            const TimeScheme &timeScheme, bool measurePoint) const {
    // options for measurement
    const double minTime = 1e4 * SimpleTimer::getClockResolution();
    const int minCount = 5;
    const int multiplyCount = 5;
    const int numSameKind = 4;
    
    // lambdas
    // element signature
    auto elemSign = [this](int iquad) -> std::string {
        std::shared_ptr<const Element> elem = this->mQuads[iquad].getElement();
        return elem->costSignature();
    };
    
    // element cost
    auto elemCost = [this](int iquad, int count) -> double {
        std::shared_ptr<const Element> elem = this->mQuads[iquad].getElement();
        return elem->measure(count);
    };
    
    // point signature
    auto pointSign = [this](int igll) -> std::string {
        std::shared_ptr<SolidPoint> sp = mGLLPoints[igll].getSolidPoint();
        std::shared_ptr<FluidPoint> fp = mGLLPoints[igll].getFluidPoint();
        std::string psign;
        if (sp && fp) {
            // on solid-fluid boundary
            psign = sp->costSignature() + "+" + fp->costSignature();
        } else {
            psign = sp ? sp->costSignature() : fp->costSignature();
        }
        return psign;
    };
    
    // point cost
    auto pointCost = [this, &timeScheme](int igll, int count) -> double {
        std::shared_ptr<SolidPoint> sp = mGLLPoints[igll].getSolidPoint();
        std::shared_ptr<FluidPoint> fp = mGLLPoints[igll].getFluidPoint();
        return ((sp ? point_time::measure(*sp, count, timeScheme) : 0.) +
                (fp ? point_time::measure(*fp, count, timeScheme) : 0.));
    };
    
    
    ////////// measure elements //////////
    timer::gPreloopTimer.begin("Measuring element cost library");
    std::map<std::string, double> elemCostLibrary;
    for (int iquad = 0; iquad < mQuads.size(); iquad++) {
        // get cost signature
        const std::string &esign = elemSign(iquad);
        // create in library
        elemCostLibrary.insert({esign, -1.});
        // perform measurement only if it is new (measured = -1.)
        if (elemCostLibrary.at(esign) > 0.) {
            continue;
        }
        // increase measured count until measured time > minimum time
        int count = minCount;
        double measured = elemCost(iquad, count);
        while (measured < minTime) {
            count *= multiplyCount;
            measured = elemCost(iquad, count);
        }
        // store in library
        elemCostLibrary.at(esign) = measured / count;
        // measure elements with the same signature
        int sameKindFound = 0;
        for (int jquad = iquad + 1; jquad < mQuads.size(); jquad++) {
            // check cost signature
            if (elemSign(jquad) == esign) {
                // measure
                double jmeasured = elemCost(jquad, count) / count;
                // use minimum
                elemCostLibrary.at(esign) =
                std::min(jmeasured, elemCostLibrary.at(esign));
                if (++sameKindFound == numSameKind) {
                    break;
                }
            }
        }
    }
    timer::gPreloopTimer.ended("Measuring element cost library");
    
    ////////// measure points //////////
    timer::gPreloopTimer.begin("Measuring point cost library");
    std::map<std::string, double> pointCostLibrary;
    for (int igll = 0; igll < mGLLPoints.size(); igll++) {
        // get cost signature
        const std::string &psign = pointSign(igll);
        // create in library
        pointCostLibrary.insert({psign, -1.});
        // perform measurement only if it is new (measured = -1.)
        if (pointCostLibrary.at(psign) > 0.) {
            continue;
        }
        // increase measured count until measured time > minimum time
        int count = minCount;
        double measured = pointCost(igll, count);
        while (measured < minTime) {
            count *= multiplyCount;
            measured = pointCost(igll, count);
        }
        // store in library
        pointCostLibrary.at(psign) = measured / count;
        // measure points with the same signature
        int sameKindFound = 0;
        for (int jgll = igll + 1; jgll < mGLLPoints.size(); jgll++) {
            // check cost signature
            if (pointSign(jgll) == psign) {
                // measure
                double jmeasured = pointCost(jgll, count) / count;
                // use minimum
                pointCostLibrary.at(psign) =
                std::min(jmeasured, pointCostLibrary.at(psign));
                if (++sameKindFound == numSameKind) {
                    break;
                }
            }
        }
    }
    timer::gPreloopTimer.ended("Measuring point cost library");
    
    // gather libraries using minimum measurements
    timer::gPreloopTimer.begin("Gathering cost libraries");
    mpi::aggregate(elemCostLibrary, 0, MPI_MIN);
    mpi::aggregate(pointCostLibrary, 0, MPI_MIN);
    // write libraries to file
    if (mpi::root() && timer::gPreloopTimer.isActivated()) {
        const std::string &fname = (io::gOutputDirectory +
                                    "/develop/measured_element_costs.log");
        std::ofstream ofs(fname);
        if (!ofs) {
            throw std::runtime_error("SE_Model::measureCost || Error opening "
                                     "or creating log file: || " + fname);
        }
        ofs << "Elements:\n";
        for (auto it = elemCostLibrary.begin();
             it != elemCostLibrary.end(); ++it) {
            ofs << it->first << " " << it->second << "\n";
        }
        ofs << "\nGLL points:\n";
        for (auto it = pointCostLibrary.begin();
             it != pointCostLibrary.end(); ++it) {
            ofs << it->first << " " << it->second << "\n";
        }
        if (!measurePoint) {
            ofs << "* Point costs are not accounted in domain partitioning.\n";
        }
        timer::gPreloopTimer.message("Cost libraries written to " + fname);
    }
    timer::gPreloopTimer.ended("Gathering cost libraries");
    
    // read libraries locally
    timer::gPreloopTimer.begin("Computing weights from libraries");
    eigen::DColX weightsLocal(mQuads.size());
    for (int iquad = 0; iquad < mQuads.size(); iquad++) {
        // element
        weightsLocal(iquad) = elemCostLibrary.at(elemSign(iquad));
        // points on this elements
        if (measurePoint) {
            for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
                int igll = localMesh.mElementGLL(iquad, ipnt);
                weightsLocal(iquad) += (pointCostLibrary.at(pointSign(igll))
                                        / mGLLPoints[igll].getElementCount());
            }
        }
    }
    timer::gPreloopTimer.ended("Computing weights from libraries");
    
    // gather weights
    timer::gPreloopTimer.begin("Gathering and broadcasting weights");
    // gather
    std::vector<eigen::DColX> weightsLocalVec;
    std::vector<eigen::IColX> L2G_ElementVec;
    mpi::gatherEigen(weightsLocal, weightsLocalVec, 0);
    mpi::gatherEigen(localMesh.mL2G_Element, L2G_ElementVec, 0);
    // map global on root
    eigen::DColX weightsGlobal;
    if (mpi::root()) {
        // allocate and fill with a crazy number
        weightsGlobal = eigen::DColX::Constant(exodusMesh.getNumQuads(),
                                               numerical::dErr);
        for (int iproc = 0; iproc < mpi::nproc(); iproc++) {
            for (int ielem = 0; ielem < L2G_ElementVec[iproc].size(); ielem++) {
                int iglob = L2G_ElementVec[iproc][ielem];
                double weight = weightsLocalVec[iproc][ielem];
                weightsGlobal(iglob) = weight;
            }
        }
    }
    mpi::bcastEigen(weightsGlobal);
    timer::gPreloopTimer.ended("Gathering and broadcasting weights");
    return weightsGlobal;
}

// verbose
std::string SE_Model::verbose(const std::string &titleSEM) const {
    // verbose on rank with maximum number of quads
    int nQuad = (int)mQuads.size();
    int nGLL = (int)mGLLPoints.size();
    std::vector<int> nQuadV;
    std::vector<int> nGLLV;
    mpi::gather(nQuad, nQuadV, MPI_ALL);
    mpi::gather(nGLL, nGLLV, MPI_ALL);
    // find rank with max
    Eigen::Index rankMax;
    Eigen::Map<eigen::IColX>(nQuadV.data(), nQuadV.size()).maxCoeff(&rankMax);
    
    // verbose
    using namespace bstring;
    std::stringstream ss;
    ss << boxTitle(titleSEM);
    ss << boxEquals(0, 24, "Rank with max # elements", rankMax);
    ss << boxEquals(0, 24, "# elements on the rank", nQuadV[rankMax]);
    ss << boxEquals(0, 24, "# GLL points on the rank", nGLLV[rankMax]);
    ss << "* Not much info here but this was a time-consuming step.\n";
    ss << boxBaseline() << "\n\n";
    return ss.str();
}


////////////////// source-receiver //////////////////
// get total undulation
double SE_Model::getTotalUndulation(const eigen::DRow3 &spzRef) const {
    double totalUndulation = 0.;
    static eigen::DColX undulation(1);
    for (const std::shared_ptr<const Geometric3D> &m: mModelsG3D) {
        m->getUndulation(spzRef, undulation);
        totalUndulation += undulation(0);
    }
    return totalUndulation;
}

// undulated to reference
eigen::DRow3 SE_Model::undulatedToReference(const eigen::DRow3 &spzUnd) const {
    // no geometric model
    if (mModelsG3D.size() == 0) {
        return spzUnd;
    }
    
    // lambda to compute R from spz
    static auto toR = geodesy::isCartesian() ?
    [](const eigen::DRow3 &spz) -> double {return spz(2);} :
    [](const eigen::DRow3 &spz) -> double {
        return sqrt(spz(0) * spz(0) + spz(2) * spz(2));};
    
    // lambda to compute spz from R
    static auto toSPZ = geodesy::isCartesian() ?
    [](double R, eigen::DRow3 &spz, double sint, double cost) {spz(2) = R;} :
    [](double R, eigen::DRow3 &spz, double sint, double cost) {
        spz(0) = R * sint; spz(2) = R * cost;};
    
    // compute theta
    double sint = 0., cost = 0.;
    if (!geodesy::isCartesian()) {
        double theta = geodesy::sz2rtheta(spzUnd, true, 0, 2, 2, 0)(0);
        sint = sin(theta);
        cost = cos(theta);
    }
    
    /////////////////// solve the inverse problem ///////////////////
    // target
    double undR = toR(spzUnd);
    
    // intitial guess
    eigen::DRow3 spzRef = spzUnd;
    if (undR < geodesy::getInnerRadius()) {
        toSPZ(geodesy::getInnerRadius(), spzRef, sint, cost);
    }
    if (undR > geodesy::getOuterRadius()) {
        toSPZ(geodesy::getOuterRadius(), spzRef, sint, cost);
    }
    
    // iteration
    double current = toR(spzRef);
    double upper = geodesy::getOuterRadius();
    double lower = geodesy::getInnerRadius();
    int maxIter = 100;
    int iter = 0;
    while (iter++ <= maxIter) {
        double diff = getTotalUndulation(spzRef) + current - undR;
        if (std::abs(diff) < mDistToleranceMesh) {
            return spzRef;
        }
        // divide and search
        if (diff > 0.) {
            upper = current;
        } else {
            lower = current;
        }
        current = .5 * (lower + upper);
        toSPZ(current, spzRef, sint, cost);
    }
    throw std::runtime_error("SE_Model::undulatedToReference || "
                             "Failed to find reference radius.");
}

// form inplane RTree
void SE_Model::formInplaneRTree() {
    mRTreeFluid = std::make_unique<RTreeND<2, 1, int>>();
    mRTreeSolid = std::make_unique<RTreeND<2, 1, int>>();
    for (int iquad = 0; iquad < mQuads.size(); iquad++) {
        const eigen::DCol2 &sz = mQuads[iquad].getNodalSZ().rowwise().mean();
        if (mQuads[iquad].fluid()) {
            mRTreeFluid->addLeaf(sz, iquad);
        } else {
            mRTreeSolid->addLeaf(sz, iquad);
        }
    }
}

// locate inplane
int SE_Model::locateInplane(const eigen::DCol2 &sz, bool inFluid) const {
    // mesh range
    if ((sz(0) < mRangeS(0) || sz(0) > mRangeS(1)) ||
        (sz(1) < mRangeZ(0) || sz(1) > mRangeZ(1))) {
        return -1;
    }
    const eigen::DCol2 &rt = geodesy::sz2rtheta(sz, false);
    if ((rt(0) < mRangeR(0) || rt(0) > mRangeR(1)) ||
        (rt(1) * rt(0) < mRangeT(0) * rt(0) - mDistToleranceMesh ||
         rt(1) * rt(0) > mRangeT(1) * rt(0) + mDistToleranceMesh)) {
        return -1;
    }
    
    // find the nearest
    std::vector<double> dists;
    std::vector<Eigen::Matrix<int, 1, 1>> vals;
    // one point is shared by at most 6 elements
    // However, it is possible that:
    // (a) the element the target point falls in has a relatively large size
    // (b) the target point is located near the edge of that host element
    // In the case of (a) && (b), the centroid of the host element will be
    // further away from the target point than the other non-hosting elements.
    // Therefore, we have to use a larger candidate pool
    // -- search the host element among 30 nearby candidates
    const int nNear = 30;
    if (inFluid) {
        mRTreeFluid->query(sz, nNear, dists, vals);
    } else {
        mRTreeSolid->query(sz, nNear, dists, vals);
    }
    
    // locate inplane
    static eigen::DCol2 xieta;
    for (int iq = 0; iq < vals.size(); iq++) {
        int iquad = vals[iq](0);
        if (mQuads[iquad].inverseMapping(sz, xieta)) {
            return iquad;
        }
    }
    // not found
    return -1;
}

// compute inplane factor
eigen::DRowN SE_Model::
computeInplaneFactor(const eigen::DCol2 &sz, int iquad) const {
    static eigen::DCol2 xieta;
    if (mQuads[iquad].inverseMapping(sz, xieta)) {
        const eigen::DColP &weightsXi =
        interpLagrange(xieta(0), mQuads[iquad].axial() ?
                       spectrals::gPositionsGLJ : spectrals::gPositionsGLL);
        const eigen::DColP &weightsEta =
        interpLagrange(xieta(1), spectrals::gPositionsGLL);
        return Eigen::Map<const eigen::DRowN>
        (eigen::DMatPP_RM(weightsXi * weightsEta.transpose()).data());
    } else {
        throw std::runtime_error("SE_Model::computeInterpFactor || "
                                 "Location is not in the given quad.");
    }
}


////////////////// wavefield scanning //////////////////
// initialize scanning on points
void SE_Model::initScanningOnPoints(bool vertexOnly) const {
    if (vertexOnly) {
        // enable all vertex points
        for (const Quad &quad: mQuads) {
            for (int ivt = 0; ivt < 4; ivt++) {
                int ipnt = vicinity::constants::gNodeIPnt[ivt];
                quad.getElement()->getPoint(ipnt).enableScanning();
            }
        }
        // disable fluid points on solid-fluid boundary
        for (const GLLPoint &point: mGLLPoints) {
            if (point.getSolidPoint() && point.getFluidPoint()) {
                point.getFluidPoint()->disableScanning();
            }
        }
    } else {
        // enable all vertex points
        for (const GLLPoint &point: mGLLPoints) {
            if (point.getSolidPoint()) {
                point.getSolidPoint()->enableScanning();
            } else {
                // only enable fluid if not on solid-fluid boundary
                point.getFluidPoint()->enableScanning();
            }
        }
    }
}
