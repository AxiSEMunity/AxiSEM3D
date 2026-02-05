//
//  WavefieldInjection.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 8/20/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  wavefield injection

#include "WavefieldInjection.hpp"
#include "Domain.hpp"
#include "Point.hpp"
#include "vector_tools.hpp"
#include "vicinity.hpp"
#include "SourceTimeFunctionNetCDF.cpp"

// destructor
WavefieldInjection::~WavefieldInjection() {
    if (mReaderSolid) {
        mReaderSolid->close();
        mReaderFluid->close();
    }
}

// initialize elements by in-quad and ex-quad tags
void WavefieldInjection::
initializeElements(std::vector<int> &interiorQuadTags,
                   std::vector<int> &exteriorQuadTags,
                   const Domain &domain) {
    if (mElementsInitialized) {
        throw std::runtime_error("WavefieldInjection::initializeElements || "
                                 "Re-initialization is forbidden.");
    }
    
    // sort and unique quad tags
    vector_tools::sortUnique(interiorQuadTags);
    vector_tools::sortUnique(exteriorQuadTags);
    
    // point indices on an element boundary
    const std::vector<int> &sidePoints = vicinity::constants::gEdgeIPntAll;
    
    // a pool of exterior point tags for search
    std::vector<int> exPointMeshTags;
    // two pools of interior elements for search
    std::vector<std::shared_ptr<const SolidElement>> inSolidElems;
    std::vector<std::shared_ptr<const FluidElement>> inFluidElems;
    // loop over solid elements
    for (const auto &elem: domain.getSolidElements()) {
        // interior
        if (vector_tools::findSortedUnique(interiorQuadTags,
                                           elem->getQuadTag())) {
            inSolidElems.push_back(elem);
            continue; // cannot be both exterior and interior
        }
        // exterior
        if (vector_tools::findSortedUnique(exteriorQuadTags,
                                           elem->getQuadTag())) {
            for (int pnt: sidePoints) {
                exPointMeshTags.push_back(elem->getPoint(pnt).getMeshTag());
            }
        }
    }
    // loop over fluid elements
    for (const auto &elem: domain.getFluidElements()) {
        // interior
        if (vector_tools::findSortedUnique(interiorQuadTags,
                                           elem->getQuadTag())) {
            inFluidElems.push_back(elem);
            continue; // cannot be both exterior and interior
        }
        // exterior
        if (vector_tools::findSortedUnique(exteriorQuadTags,
                                           elem->getQuadTag())) {
            for (int pnt: sidePoints) {
                exPointMeshTags.push_back(elem->getPoint(pnt).getMeshTag());
            }
        }
    }
    // sort exterior point tags
    vector_tools::sortUnique(exPointMeshTags);
    
    // interior solid elements on the injection boundary
    for (const auto &inElem: inSolidElems) {
        // search boundary points
        std::vector<int> pointsOnWJB =
        inElem->findBoundaryPointsByTag(exPointMeshTags);
        // found an interior element
        if (pointsOnWJB.size() != 0) {
            std::shared_ptr<SolidElementInteriorWJ> inElemWJ =
            std::make_shared<SolidElementInteriorWJ>(*inElem, pointsOnWJB);
            mInteriorSolidElements.push_back(inElemWJ);
        }
    }
    
    // interior fluid elements on the injection boundary
    for (const auto &inElem: inFluidElems) {
        // search boundary points
        std::vector<int> pointsOnWJB =
        inElem->findBoundaryPointsByTag(exPointMeshTags);
        // found an interior element
        if (pointsOnWJB.size() != 0) {
            std::shared_ptr<FluidElementInteriorWJ> inElemWJ =
            std::make_shared<FluidElementInteriorWJ>(*inElem, pointsOnWJB);
            mInteriorFluidElements.push_back(inElemWJ);
        }
    }
    
    // initialized
    mElementsInitialized = true;
}

// initialize elements by in-quad tags and boundary coordinates
void WavefieldInjection::
initializeElements(std::vector<int> &interiorQuadTags,
                   const std::vector<double> &boundaryCrdsRorZ,
                   const std::vector<double> &boundaryCrdsTorS,
                   double distTol, const Domain &domain) {
    if (mElementsInitialized) {
        throw std::runtime_error("WavefieldInjection::initializeElements || "
                                 "Re-initialization is forbidden.");
    }
    
    // sort and unique quad tags
    vector_tools::sortUnique(interiorQuadTags);
    
    // two pools of interior elements for search
    std::vector<std::shared_ptr<const SolidElement>> inSolidElems;
    std::vector<std::shared_ptr<const FluidElement>> inFluidElems;
    // loop over solid elements
    for (const auto &elem: domain.getSolidElements()) {
        // interior
        if (vector_tools::findSortedUnique(interiorQuadTags,
                                           elem->getQuadTag())) {
            inSolidElems.push_back(elem);
        }
    }
    // loop over fluid elements
    for (const auto &elem: domain.getFluidElements()) {
        // interior
        if (vector_tools::findSortedUnique(interiorQuadTags,
                                           elem->getQuadTag())) {
            inFluidElems.push_back(elem);
        }
    }
    
    // interior solid elements on the injection boundary
    for (const auto &inElem: inSolidElems) {
        // search boundary points
        std::vector<int> pointsOnWJB =
        inElem->findBoundaryPointsByCrds(boundaryCrdsRorZ, boundaryCrdsTorS,
                                         distTol);
        // found an interior element
        if (pointsOnWJB.size() != 0) {
            std::shared_ptr<SolidElementInteriorWJ> inElemWJ =
            std::make_shared<SolidElementInteriorWJ>(*inElem, pointsOnWJB);
            mInteriorSolidElements.push_back(inElemWJ);
        }
    }
    
    // interior fluid elements on the injection boundary
    for (const auto &inElem: inFluidElems) {
        // search boundary points
        std::vector<int> pointsOnWJB =
        inElem->findBoundaryPointsByCrds(boundaryCrdsRorZ, boundaryCrdsTorS,
                                         distTol);
        // found an interior element
        if (pointsOnWJB.size() != 0) {
            std::shared_ptr<FluidElementInteriorWJ> inElemWJ =
            std::make_shared<FluidElementInteriorWJ>(*inElem, pointsOnWJB);
            mInteriorFluidElements.push_back(inElemWJ);
        }
    }
    
    // initialized
    mElementsInitialized = true;
}

// initialize source-time functions
void WavefieldInjection::
initializeSTFs(const std::string &ncFileNameSolid,
               const std::string &ncFileNameFluid,
               bool aligned, int bufferSize) {
    if (mSTFsInitialized) {
        throw std::runtime_error("WavefieldInjection::initializeSTFs || "
                                 "Re-initialization is forbidden.");
    }
    
    // no injection
    if (mInteriorSolidElements.size() +
        mInteriorFluidElements.size() == 0) {
        mSTFsInitialized = true;
        return;
    }
    
    // file
    mReaderSolid = std::make_shared<NetCDF_Reader>();
    mReaderFluid = std::make_shared<NetCDF_Reader>();
#ifdef _USE_PARALLEL_NETCDF
    mReaderSolid->openParallel(ncFileNameSolid);
    mReaderFluid->openParallel(ncFileNameFluid);
#else
    mReaderSolid->open(ncFileNameSolid);
    mReaderFluid->open(ncFileNameFluid);
#endif
    
    // solid
    for (const auto &elem: mInteriorSolidElements) {
        elem->setSTF(aligned, bufferSize, mReaderSolid);
    }
    
    // fluid
    for (const auto &elem: mInteriorFluidElements) {
        elem->setSTF(aligned, bufferSize, mReaderFluid);
    }
    
    // initialized
    mSTFsInitialized = true;
}

// set in domain
void WavefieldInjection::setInDomain(Domain &domain) const {
    if (!mElementsInitialized || !mSTFsInitialized) {
        throw std::runtime_error("WavefieldInjection::setInDomain || "
                                 "Initialization is not done.");
    }
    
    // solid
    for (const auto &elem: mInteriorSolidElements) {
        domain.replaceSolidElement(elem);
    }
    
    // fluid
    for (const auto &elem: mInteriorFluidElements) {
        domain.replaceFluidElement(elem);
    }
}

// apply
void WavefieldInjection::
apply(int timeStep, double time) const {
    // solid
    for (const auto &elem: mInteriorSolidElements) {
        elem->applyInjection(timeStep, time);
    }
    
    // fluid
    for (const auto &elem: mInteriorFluidElements) {
        elem->applyInjection(timeStep, time);
    }
}

// info
void WavefieldInjection::
infoSolid(std::vector<std::string> &quadKeys,
          std::vector<int> &quadNrs,
          std::vector<std::vector<int>> &quadBoundaryPnts,
          std::vector<std::vector<std::string>> &recKeys,
          std::vector<eigen::DMatXX> &rec_spz) const {
    // size
    quadKeys.resize(mInteriorSolidElements.size());
    quadNrs.resize(mInteriorSolidElements.size());
    quadBoundaryPnts.resize(mInteriorSolidElements.size());
    recKeys.resize(mInteriorSolidElements.size());
    rec_spz.resize(mInteriorSolidElements.size());
    // key
    for (int ielem = 0; ielem < mInteriorSolidElements.size(); ielem++) {
        quadKeys[ielem] = mInteriorSolidElements[ielem]->quadKey();
        quadNrs[ielem] = mInteriorSolidElements[ielem]->getNr();
        quadBoundaryPnts[ielem] = mInteriorSolidElements[ielem]->pointsOnWJB();
        mInteriorSolidElements[ielem]->receiverInfo(recKeys[ielem],
                                                    rec_spz[ielem]);
    }
}

// info
void WavefieldInjection::
infoFluid(std::vector<std::string> &quadKeys,
          std::vector<int> &quadNrs,
          std::vector<std::vector<int>> &quadBoundaryPnts,
          std::vector<std::vector<std::string>> &recKeys,
          std::vector<eigen::DMatXX> &rec_spz) const {
    // size
    quadKeys.resize(mInteriorFluidElements.size());
    quadNrs.resize(mInteriorFluidElements.size());
    quadBoundaryPnts.resize(mInteriorFluidElements.size());
    recKeys.resize(mInteriorFluidElements.size());
    rec_spz.resize(mInteriorFluidElements.size());
    // key
    for (int ielem = 0; ielem < mInteriorFluidElements.size(); ielem++) {
        quadKeys[ielem] = mInteriorFluidElements[ielem]->quadKey();
        quadNrs[ielem] = mInteriorFluidElements[ielem]->getNr();
        quadBoundaryPnts[ielem] = mInteriorFluidElements[ielem]->pointsOnWJB();
        mInteriorFluidElements[ielem]->receiverInfo(recKeys[ielem],
                                                    rec_spz[ielem]);
    }
}
