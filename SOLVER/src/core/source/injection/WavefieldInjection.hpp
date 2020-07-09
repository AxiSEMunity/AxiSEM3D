//
//  WavefieldInjection.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 8/20/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  wavefield injection

#ifndef WavefieldInjection_hpp
#define WavefieldInjection_hpp

#include "ElementInteriorWJ.hpp"
class Domain;

class WavefieldInjection {
    // element types
    typedef ElementInteriorWJ<SolidElement, 3> SolidElementInteriorWJ;
    typedef ElementInteriorWJ<FluidElement, 1> FluidElementInteriorWJ;
    
public:
    // destructor
    ~WavefieldInjection();
    
    // initialize elements by in-quad and ex-quad tags
    void initializeElements(std::vector<int> &interiorQuadTags,
                            std::vector<int> &exteriorQuadTags,
                            const Domain &domain);
    
    // initialize elements by in-quad tags and boundary coordinates
    void initializeElements(std::vector<int> &interiorQuadTags,
                            const std::vector<double> &boundaryCrdsRorZ,
                            const std::vector<double> &boundaryCrdsTorS,
                            double distTol, const Domain &domain);
    
    // initialize source-time functions
    void initializeSTFs(const std::string &ncFileNameSolid,
                        const std::string &ncFileNameFluid,
                        bool aligned, int bufferSize);
    
    // set in domain
    void setInDomain(Domain &domain) const;
    
    // apply
    void apply(int timeStep, double time) const;
    
    // info in solid
    void infoSolid(std::vector<std::string> &quadKeys,
                   std::vector<int> &quadNrs,
                   std::vector<std::vector<int>> &quadBoundaryPoints,
                   std::vector<std::vector<std::string>> &recKeys,
                   std::vector<eigen::DMatXX> &rec_spz) const;
    
    // info in fluid
    void infoFluid(std::vector<std::string> &quadKeys,
                   std::vector<int> &quadNrs,
                   std::vector<std::vector<int>> &quadBoundaryPnts,
                   std::vector<std::vector<std::string>> &recKeys,
                   std::vector<eigen::DMatXX> &rec_spz) const;
    
private:
    // elements
    std::vector<std::shared_ptr<SolidElementInteriorWJ>> mInteriorSolidElements;
    std::vector<std::shared_ptr<FluidElementInteriorWJ>> mInteriorFluidElements;
    
    // file
    std::shared_ptr<NetCDF_Reader> mReaderSolid = nullptr;
    std::shared_ptr<NetCDF_Reader> mReaderFluid = nullptr;
    
    // flags
    bool mElementsInitialized = false;
    bool mSTFsInitialized = false;
};

#endif /* WavefieldInjection_hpp */
