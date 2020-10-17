//
//  Domain.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 8/1/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  computational domain

#include "Domain.hpp"
// mesh
#include "SolidPoint.hpp"
#include "FluidPoint.hpp"
#include "SolidElement.hpp"
#include "FluidElement.hpp"
// boundary
#include "SolidFluidBoundary.hpp"
#include "AbsorbingBoundary.hpp"
#include "AxialBoundary.hpp"
#include "FluidSurfaceBoundary.hpp"
// mpi
#include "Messaging.hpp"
// source
#include "WavefieldInjection.hpp"
#include "ElementSource.hpp"
// output
#include "StationGroup.hpp"
#include "ElementOpGroup.hpp"
// wavefield scanning
#include "WavefieldScanning.hpp"
#include "NetCDF_Writer.hpp"
// verbose
#include "io.hpp"
#include "mpi.hpp"
#include "bstring.hpp"
#include "vector_tools.hpp"
#include "geodesy.hpp"

// constructor
Domain::Domain() {
    // set all single unique_ptr's to nullptr
    // cannot do this in .hpp as it requires #include
    mSolidFluidBoundary = std::make_unique<SolidFluidBoundary>();
    mAbsorbingBoundary = std::make_unique<AbsorbingBoundary>();
    mAxialBoundary = std::make_unique<AxialBoundary>();
    mFluidSurfaceBoundary = std::make_unique<FluidSurfaceBoundary>();
    mMessaging = nullptr;
    mWavefieldInjection = nullptr;
    mWavefieldScanning = nullptr;
}

////////////////////// domain construction //////////////////////
// add a solid point
void Domain::addSolidPoint(const std::shared_ptr<SolidPoint> &point) {
    point->setDomainTag((int)mSolidPoints.size());
    mSolidPoints.push_back(point);
}

// add a fluid point
void Domain::addFluidPoint(const std::shared_ptr<FluidPoint> &point) {
    point->setDomainTag((int)mFluidPoints.size());
    mFluidPoints.push_back(point);
}

// add a solid element
void Domain::addSolidElement(const std::shared_ptr<SolidElement> &element) {
    element->setDomainTag((int)mSolidElements.size());
    mSolidElements.push_back(element);
}

// add a fluid element
void Domain::addFluidElement(const std::shared_ptr<FluidElement> &element) {
    element->setDomainTag((int)mFluidElements.size());
    mFluidElements.push_back(element);
}

// replace a solid element
void Domain::replaceSolidElement(const std::shared_ptr<SolidElement> &element) {
    // check domain tag
    int domainTag = element->getDomainTag();
    if (domainTag == -1) {
        throw std::runtime_error("Domain::replaceSolidElement || "
                                 "Domain tag must be set before replacement.");
    }
    // check quad tag
    if (element->getQuadTag() != mSolidElements[domainTag]->getQuadTag()) {
        throw std::runtime_error("Domain::replaceSolidElement || "
                                 "The new and the old elements must have "
                                 "the same Quad tag.");
    }
    //  check old reference
    if (mSolidElements[domainTag].use_count() > 1) {
        throw std::runtime_error("Domain::replaceSolidElement || "
                                 "The old element has been referred outside "
                                 "this domain and is thus irreplaceable.");
    }
    // after replacement, the old will be deleted
    mSolidElements[element->getDomainTag()] = element;
}

// replace a fluid element
void Domain::replaceFluidElement(const std::shared_ptr<FluidElement> &element) {
    // check domain tag
    int domainTag = element->getDomainTag();
    if (domainTag == -1) {
        throw std::runtime_error("Domain::replaceFluidElement || "
                                 "Domain tag must be set before replacement.");
    }
    // check quad tag
    if (element->getQuadTag() != mFluidElements[domainTag]->getQuadTag()) {
        throw std::runtime_error("Domain::replaceFluidElement || "
                                 "The new and the old elements must have "
                                 "the same Quad tag.");
    }
    //  check old reference
    if (mFluidElements[domainTag].use_count() > 1) {
        throw std::runtime_error("Domain::replaceFluidElement || "
                                 "The old element has been referred outside "
                                 "this domain and is thus irreplaceable.");
    }
    // after replacement, the old will be deleted
    mFluidElements[element->getDomainTag()] = element;
}

// set mpi messaging
void Domain::setMessaging(std::unique_ptr<Messaging> &msg) {
    mMessaging = std::move(msg);
}

// verbose
std::string Domain::verbose(const std::string &title) const {
    //////////////////////////// count info ////////////////////////////
    // point
    std::map<std::string, int> typeCountPoint;
    for (const auto &point: mSolidPoints) {
        if (!mMessaging->pointInSmallerRank(point)) {
            vector_tools::aggregate(typeCountPoint, point->typeInfo(), 1);
        }
    }
    for (const auto &point: mFluidPoints) {
        if (!mMessaging->pointInSmallerRank(point)) {
            vector_tools::aggregate(typeCountPoint, point->typeInfo(), 1);
        }
    }
    mpi::aggregate(typeCountPoint, 0, MPI_SUM);
    
    // element
    std::map<std::string, int> typeCountElement;
    for (const auto &elem: mSolidElements) {
        vector_tools::aggregate(typeCountElement, elem->typeInfo(), 1);
    }
    for (const auto &elem: mFluidElements) {
        vector_tools::aggregate(typeCountElement, elem->typeInfo(), 1);
    }
    mpi::aggregate(typeCountElement, 0, MPI_SUM);
    
    // solid-fluid boundary
    std::map<std::string, int> typeCountSFB =
    mSolidFluidBoundary->countInfo(*mMessaging);
    mpi::aggregate(typeCountSFB, 0, MPI_SUM);
    
    // absorbing boundary
    std::map<std::string, int> typeCountABB =
    mAbsorbingBoundary->countInfo(*mMessaging);
    mpi::aggregate(typeCountABB, 0, MPI_SUM);
    
    // axial boundary
    int axSolid = 0, axFluid = 0;
    mAxialBoundary->countInfo(*mMessaging, axSolid, axFluid);
    axSolid = mpi::sum(axSolid);
    axFluid = mpi::sum(axFluid);
    
    // fluid surface boundary
    int fluidSurf = mFluidSurfaceBoundary->countInfo(*mMessaging);
    fluidSurf = mpi::sum(fluidSurf);
    
    // mpi
    int nComm = mMessaging->getNumRankComm();
    nComm = mpi::sum(nComm);
    nComm /= 2; // one comm counted by two ranks
    
    
    //////////////////////////// max key size ////////////////////////////
    int mkl = (int)std::string("# rank-to-rank communications").size();
    mkl = std::max(mkl, vector_tools::maxKeyLength(typeCountPoint));
    mkl = std::max(mkl, vector_tools::maxKeyLength(typeCountElement));
    mkl = std::max(mkl, vector_tools::maxKeyLength(typeCountSFB));
    mkl = std::max(mkl, vector_tools::maxKeyLength(typeCountABB));
    
    
    //////////////////////////// write box ////////////////////////////
    std::stringstream ss;
    ss << bstring::boxTitle(title);
    // point
    int total = vector_tools::sumValues(typeCountPoint);
    ss << bstring::boxSubTitle(0, "GLL points");
    ss << bstring::boxEquals(2, mkl, typeCountPoint);
    ss << bstring::boxEquals(2, mkl + 1, "Σ", total);
    // element
    total = vector_tools::sumValues(typeCountElement);
    ss << bstring::boxSubTitle(0, "Spectral elements");
    ss << bstring::boxEquals(2, mkl, typeCountElement);
    ss << bstring::boxEquals(2, mkl + 1, "Σ", total);
    // axial boundary
    ss << bstring::boxSubTitle(0, "Axial boundary");
    ss << bstring::boxEquals(2, mkl, "Solid", axSolid);
    ss << bstring::boxEquals(2, mkl, "Fluid", axFluid);
    ss << bstring::boxEquals(2, mkl + 1, "Σ ", axSolid + axFluid);
    // absorbing boundary
    total = vector_tools::sumValues(typeCountABB);
    ss << bstring::boxSubTitle(0, "Absorbing boundary");
    if (total > 0) {
        ss << bstring::boxEquals(2, mkl, typeCountABB);
    }
    ss << bstring::boxEquals(2, mkl + 1, "Σ", total);
    // fluid surface boundary
    ss << bstring::boxSubTitle(0, "Fluid surface boundary");
    ss << bstring::boxEquals(2, mkl + 1, "Σ", fluidSurf);
    // solid-fluid boundary
    total = vector_tools::sumValues(typeCountSFB);
    ss << bstring::boxSubTitle(0, "Solid-fluid boundary");
    if (total > 0) {
        ss << bstring::boxEquals(2, mkl, typeCountSFB);
    }
    ss << bstring::boxEquals(2, mkl + 1, "Σ", total);
    // mpi
    ss << bstring::boxSubTitle(0, "Domain decomposition");
    ss << bstring::boxEquals(2, mkl, "# sub-domains (nproc)", mpi::nproc());
    ss << bstring::boxEquals(2, mkl, "# rank-to-rank communications", nComm);
    // end
    ss << bstring::boxBaseline() << "\n\n";
    return ss.str();
}

// MESHING STAGE //
///////////////////
// SOLVING STAGE //

// set wavefield injection
void Domain::
setWavefieldInjection(std::unique_ptr<WavefieldInjection> &wj) {
    wj->setInDomain(*this);
    mWavefieldInjection = std::move(wj);
}

void Domain::
addElementSource(std::unique_ptr<const ElementSource> &esrc) {
    mElementSources.push_back(std::move(esrc));
}

// add station group in solid
void Domain::
addStationGroupInSolid(std::unique_ptr<StationGroupInSolid> &stgrp) {
    mStationGroupInSolids.push_back(std::move(stgrp));
}

// add station group in fluid
void Domain::
addStationGroupInFluid(std::unique_ptr<StationGroupInFluid> &stgrp) {
    mStationGroupInFluids.push_back(std::move(stgrp));
}

// add element group in solid
void Domain::
addElementOpGroupInSolid(std::unique_ptr<ElementOpGroupInSolid> &elgrp) {
    mElementOpGroupInSolids.push_back(std::move(elgrp));
}

// add element group in fluid
void Domain::
addElementOpGroupInFluid(std::unique_ptr<ElementOpGroupInFluid> &elgrp) {
    mElementOpGroupInFluids.push_back(std::move(elgrp));
}

// set wavefield scanning
void Domain::
setWavefieldScanning(std::unique_ptr<WavefieldScanning> &ws) {
    mWavefieldScanning = std::move(ws);
}


////////////////////// timeloop //////////////////////
// apply sources
void Domain::applySources(int tstep, double time) const {
    // wavefield injection
    if (mWavefieldInjection) {
        mWavefieldInjection->apply(tstep, time);
    }
    
    // general sources
    for (const std::unique_ptr<const ElementSource> &src: mElementSources) {
        src->apply(time);
    }
}

// compute stiffness
void Domain::computeStiffness() const {
    // solid
    for (const std::shared_ptr<const SolidElement> &element: mSolidElements) {
        element->computeStiff();
    }
    
    // fluid
    for (const std::shared_ptr<const FluidElement> &element: mFluidElements) {
        element->computeStiff();
    }
}

// boundary conditions before assembling stiffness
// * Clayton ABC
// * solid-fluild
void Domain::applyBC_BeforeAssemblingStiff() const {
    // the following order matters!
    mAbsorbingBoundary->applyClayton();
    mSolidFluidBoundary->apply();
}

// boundary conditions after assembling stiffness
// * axial
// * fluid surface
void Domain::applyBC_AfterAssemblingStiff() const {
    // the following order matters!
    mAxialBoundary->apply();
    mFluidSurfaceBoundary->apply();
}

// boundary conditions after computing acceleration
// * sponge ABC
void Domain::applyBC_AfterComputingAccel() const {
    mAbsorbingBoundary->applySponge();
}

// stiff to accel on points
void Domain::computeStiffToAccel() const {
    // solid
    for (const std::shared_ptr<SolidPoint> &point: mSolidPoints) {
        point->computeStiffToAccel();
    }
    
    // fluid
    for (const std::shared_ptr<FluidPoint> &point: mFluidPoints) {
        point->computeStiffToAccel();
    }
}

// mpi phase 1: commGatherSendRecv
void Domain::mpiGatherSendRecv() const {
    mMessaging->commGatherSendRecv();
}

// mpi phase 2: commWaitScatter
void Domain::mpiWaitScatter() const {
    mMessaging->commWaitScatter();
}

// initialize output
void Domain::initializeOutput() const {
    // solid stations
    for (const std::unique_ptr<StationGroupInSolid> &stgrp:
         mStationGroupInSolids) {
        stgrp->initialize();
    }
    
    // fluid stations
    for (const std::unique_ptr<StationGroupInFluid> &stgrp:
         mStationGroupInFluids) {
        stgrp->initialize();
    }
    
    // solid elements
    for (const std::unique_ptr<ElementOpGroupInSolid> &elgrp:
         mElementOpGroupInSolids) {
        elgrp->initialize();
    }
    
    // fluid elements
    for (const std::unique_ptr<ElementOpGroupInFluid> &elgrp:
         mElementOpGroupInFluids) {
        elgrp->initialize();
    }
}

// record output
void Domain::recordOutput(int tstep, double time) const {
    // solid stations
    for (const std::unique_ptr<StationGroupInSolid> &stgrp:
         mStationGroupInSolids) {
        stgrp->record(tstep, time);
    }
    
    // fluid stations
    for (const std::unique_ptr<StationGroupInFluid> &stgrp:
         mStationGroupInFluids) {
        stgrp->record(tstep, time);
    }
    
    // solid elements
    for (const std::unique_ptr<ElementOpGroupInSolid> &elgrp:
         mElementOpGroupInSolids) {
        elgrp->record(tstep, time);
    }
    
    // fluid elements
    for (const std::unique_ptr<ElementOpGroupInFluid> &elgrp:
         mElementOpGroupInFluids) {
        elgrp->record(tstep, time);
    }
}

// dump output
void Domain::dumpOutput() const {
    // solid stations
    for (const std::unique_ptr<StationGroupInSolid> &stgrp:
         mStationGroupInSolids) {
        stgrp->dumpToFile();
    }
    
    // fluid stations
    for (const std::unique_ptr<StationGroupInFluid> &stgrp:
         mStationGroupInFluids) {
        stgrp->dumpToFile();
    }
    
    // solid elements
    for (const std::unique_ptr<ElementOpGroupInSolid> &elgrp:
         mElementOpGroupInSolids) {
        elgrp->dumpToFile();
    }
    
    // fluid elements
    for (const std::unique_ptr<ElementOpGroupInFluid> &elgrp:
         mElementOpGroupInFluids) {
        elgrp->dumpToFile();
    }
}

// finalize output
void Domain::finalizeOutput() const {
    // solid stations
    for (const std::unique_ptr<StationGroupInSolid> &stgrp:
         mStationGroupInSolids) {
        stgrp->finalize();
    }
    
    // fluid stations
    for (const std::unique_ptr<StationGroupInFluid> &stgrp:
         mStationGroupInFluids) {
        stgrp->finalize();
    }
    
    // solid elements
    for (const std::unique_ptr<ElementOpGroupInSolid> &elgrp:
         mElementOpGroupInSolids) {
        elgrp->finalize();
    }
    
    // fluid elements
    for (const std::unique_ptr<ElementOpGroupInFluid> &elgrp:
         mElementOpGroupInFluids) {
        elgrp->finalize();
    }
}

// do scanning
void Domain::doScanning(int tstep) const {
    if (!mWavefieldScanning) {
        return;
    }
    
    // check step
    if (tstep % mWavefieldScanning->mScanningInterval != 0) {
        return;
    }
    
    // solid
    for (const std::shared_ptr<SolidPoint> &point: mSolidPoints) {
        point->doScanning(mWavefieldScanning->mTolFourierH2,
                          mWavefieldScanning->mRelTolH2,
                          mWavefieldScanning->mAbsTolH2,
                          mWavefieldScanning->mMaxNumPeaks);
    }
    
    // fluid
    for (const std::shared_ptr<FluidPoint> &point: mFluidPoints) {
        point->doScanning(mWavefieldScanning->mTolFourierH2,
                          mWavefieldScanning->mRelTolH2,
                          mWavefieldScanning->mAbsTolH2,
                          mWavefieldScanning->mMaxNumPeaks);
    }
}

// report scanning
void Domain::reportScanning() const {
    if (!mWavefieldScanning) {
        return;
    }
    
    // gather data from points
    std::vector<double> szBuffer;
    std::vector<int> nrBuffer;
    
    // solid
    for (const std::shared_ptr<SolidPoint> &point: mSolidPoints) {
        if (!mMessaging->pointInSmallerRank(point)) {
            int scanNr = point->reportScanningNr();
            if (scanNr != -1) {
                szBuffer.push_back(point->getCoords()(0));
                szBuffer.push_back(point->getCoords()(1));
                nrBuffer.push_back(point->getNr());
                nrBuffer.push_back(scanNr);
            }
        }
    }
    
    // fluid
    for (const std::shared_ptr<FluidPoint> &point: mFluidPoints) {
        if (!mMessaging->pointInSmallerRank(point)) {
            int scanNr = point->reportScanningNr();
            if (scanNr != -1) {
                szBuffer.push_back(point->getCoords()(0));
                szBuffer.push_back(point->getCoords()(1));
                nrBuffer.push_back(point->getNr());
                nrBuffer.push_back(scanNr);
            }
        }
    }
    
    // gather mpi
    std::vector<std::vector<double>> szBufferAll;
    std::vector<std::vector<int>> nrBufferAll;
    int nPoints = (int)szBuffer.size() / 2;
    mpi::gather(szBuffer, szBufferAll, 0);
    mpi::gather(nrBuffer, nrBufferAll, 0);
    nPoints = mpi::sum(nPoints);
    
    // write to file
    if (mpi::root()) {
        // flatten buffer
        eigen::DMatXX_RM sz(nPoints, 2);
        eigen::IColX nrOrig(nPoints);
        eigen::IColX nrScan(nPoints);
        nPoints = 0;
        for (int iproc = 0; iproc < mpi::nproc(); iproc++) {
            for (int ipnt = 0; ipnt < szBufferAll[iproc].size() / 2; ipnt++) {
                sz(nPoints, 0) = szBufferAll[iproc][ipnt * 2];
                sz(nPoints, 1) = szBufferAll[iproc][ipnt * 2 + 1];
                nrOrig(nPoints) = nrBufferAll[iproc][ipnt * 2];
                nrScan(nPoints) = nrBufferAll[iproc][ipnt * 2 + 1];
                nPoints++;
            }
        }
        
        // write
        NetCDF_Writer writer;
        writer.open(io::gOutputDirectory + "/" +
                    mWavefieldScanning->mFileName, true);
        writer.defModeOn();
        writer.defineVariable("pointwise_sz", {
            {"dim_point", nPoints}, {"dim_sz", 2}}, numerical::dErr);
        writer.defineVariable("pointwise_Nr", {
            {"dim_point", nPoints}}, (int)-1);
        writer.defineVariable("starting_Nr_for_scanning", {
            {"dim_point", nPoints}}, (int)-1);
        writer.defModeOff();
        writer.writeWholeVariable("pointwise_sz", sz);
        writer.writeWholeVariable("pointwise_Nr", nrScan);
        writer.writeWholeVariable("starting_Nr_for_scanning", nrOrig);
        writer.close();
    }
}

// check stability
void Domain::checkStability(int tstep, double t, double dt) const {
    std::shared_ptr<Point> unstablePoint = nullptr;
    // solid
    for (const std::shared_ptr<SolidPoint> &point: mSolidPoints) {
        if (!point->stable()) {
            unstablePoint = point;
            break;
        }
    }
    
    // fluid
    if (!unstablePoint) {
        for (const std::shared_ptr<FluidPoint> &point: mFluidPoints) {
            if (!point->stable()) {
                unstablePoint = point;
                break;
            }
        }
    }
    
    // abort if unstable
    if (unstablePoint) {
        using namespace bstring;
        std::stringstream ss;
        ss << "Simulation has blown up, with ΔT = " << dt << " s. || ";
        ss << "Where the instability occurred: || ";
        const eigen::DRow2 &sz = unstablePoint->getCoords();
        if (geodesy::isCartesian()) {
            ss << "* (s, z)  =  " << range(sz(0), sz(1), '(', ')') << " || ";
        } else {
            const eigen::DRow2 &rt = geodesy::sz2rtheta(sz, true);
            ss << "* (r, θ)  =  " << range(rt(0), rt(1), '(', ')') << " || ";
        }
        ss << "When the instability occurred: || ";
        ss << "* current time  =  " << t << " || ";
        ss << "* current step  =  " << tstep;
        throw std::runtime_error("Domain::checkStability || " + ss.str());
    }
}


////////////////////////////////////////
//////////////// static ////////////////
////////////////////////////////////////

// empty creator
// constructor not available from outside becuase of forward declarations
std::shared_ptr<Domain> Domain::createEmpty() {
    return std::make_shared<Domain>();
}
