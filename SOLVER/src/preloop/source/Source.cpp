//
//  Source.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 5/17/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  source: generator of ElementSource in core

#include "Source.hpp"

#include "SE_Model.hpp"
#include "STF.hpp"
#include "Mechanism.hpp"

#include "inparam.hpp"
#include "mpi.hpp"
#include "timer.hpp"
#include "io.hpp"

// build from inparam
std::shared_ptr<const Source> Source::buildInparam(int sindex) {
    // source name and key in inparam
    const InparamYAML &gm = inparam::gInparamSource;
    std::string root = "list_of_sources";
    const std::string &sourceName = gm.get<std::string>
    (root + ":{" + bstring::toString(sindex) + "}");
    root += ":[" + bstring::toString(sindex) + "]:" + sourceName;
    
    // horizontal
    root += ":location";
    bool axial = false, sourceCentered = false, ellipticity = true;
    eigen::DRow3 crd = eigen::DRow3::Zero();
    if (gm.contains(root + ":latitude_longitude")) {
        const std::string &rooth = root + ":latitude_longitude";
        if (gm.get<std::string>(rooth) == "ON_AXIS") {
            axial = true;
        } else {
            axial = false;
            sourceCentered = false;
            crd(0) = gm.getWithBounds(rooth + ":[0]", -90., 90.);
            crd(1) = gm.getWithBounds(rooth + ":[1]", -180., 180.);
            ellipticity = gm.get<bool>(root + ":ellipticity");
        }
    } else if (gm.contains(root + ":distance_azimuth")) {
        const std::string &rooth = root + ":distance_azimuth";
        if (gm.get<std::string>(rooth) == "ON_AXIS") {
            axial = true;
        } else {
            axial = false;
            sourceCentered = true;
            if (geodesy::isCartesian()) {
                crd(0) = gm.getWithBounds(rooth + ":[0]", 0.);
            } else {
                crd(0) = gm.getWithBounds(rooth + ":[0]", 0., numerical::dPi);
            }
            crd(1) = gm.getWithBounds(rooth + ":[1]", 0., 2. * numerical::dPi);
        }
    } else {
        throw std::runtime_error("Source::buildInparam || "
                                 "Either location:latitude_longitude or || "
                                 "location:distance_azimuth must exist. || "
                                 "Source name: " + sourceName);
    }
    
    // first source is special
    if (sindex == 0 && inparam::gInparamModel.
        get<std::string>("geodesy:lat_lon_north_pole_mesh") == "SOURCE") {
        axial = true;
    }
    
    // vertical
    bool useDepth = false, depthSolid = false;
    if (gm.contains(root + ":depth")) {
        useDepth = true;
        crd(2) = gm.get<double>(root + ":depth");
        depthSolid = gm.get<bool>(root + ":depth_below_solid_surface");
    } else if (gm.contains(root + ":radius")) {
        useDepth = false;
        crd(2) = gm.get<double>(root + ":radius");
    } else {
        throw std::runtime_error("Source::buildInparam || "
                                 "Either location:radius or "
                                 "location:depth must exist. || "
                                 "Source name: " + sourceName);
    }
    bool undulated = gm.get<bool>(root + ":undulated_geometry");
    
    // construct and return
    return std::make_shared<const Source>(sourceName, axial,
                                          sourceCentered, ellipticity,
                                          useDepth, depthSolid, undulated, crd);
}

// compute spz
eigen::DRow3 Source::
computeSPZ(const SE_Model &sem, const eigen::DRow3 &crdIn,
           bool sourceCentered, bool xy, bool ellipticity,
           bool useDepth, bool depthSolid, bool undulatedGeometry,
           const std::string &errInfo, bool enforceOnAxis) {
    // result
    static eigen::DRow3 crd;
    crd = crdIn;
    
    // vertical
    // NOTE: here vertical is approximate when using undulated geometry,
    //       but it is still needed to determine the flattening at depth
    if (useDepth) {
        double R = (depthSolid ? geodesy::getOuterSolidRadius() :
                    geodesy::getOuterRadius());
        crd(2) = R - crd(2);
    }
    
    // directly determine s and phi for an axial source
    if (enforceOnAxis) {
        crd(0) = 0.;
        crd(1) = 0.;
    } else {
        // horizontal
        if (sourceCentered) {
            if (xy) {
                // xyz -> spz, for both spherical and Cartesian
                double x = crd(0);
                double y = crd(1);
                crd(0) = sqrt(x * x + y * y);
                crd(1) = (crd(0) < numerical::dEpsilon) ? 0. : atan2(y, x);
            } else {
                if (!geodesy::isCartesian()) {
                    // (theta, phi, r) -> spz
                    double r = crd(2);
                    double theta = crd(0);
                    crd(0) = r * sin(theta);
                    crd(2) = r * cos(theta);
                }
            }
        } else {
            // geographic to source-centered
            crd = geodesy::llr2spz(crd, ellipticity);
            
            // correct spz for Cartesian
            // NOTE: Geographic and Cartesian contradict each other.
            //       We correct spz using s as arc-length and z as radius.
            //       Without such "bending", sqrt(s*s+z*z) on the surface
            //       will exceed the outer radius.
            if (geodesy::isCartesian()) {
                static eigen::DRow3 spzUnbend;
                spzUnbend = geodesy::sz2rtheta(crd, true, 0, 2, 2, 0);
                spzUnbend(0) *= spzUnbend(2);
                spzUnbend(1) = crd(1);
                // z must be significantly larger than s
                if (spzUnbend(2) < spzUnbend(0) * 10.) {
                    throw std::runtime_error
                    ("Source::computeSPZ || "
                     "Invalid geographic location in Cartesian mesh. || "
                     "Set PLANET_RADIUS in Salvus mesher data file (*.bm). || "
                     + errInfo);
                }
                crd = spzUnbend;
            }
        }
    }
    
    // determine accurate vertical location in undulated geometry
    if (undulatedGeometry) {
        // use the undulated surface for depth
        if (useDepth) {
            // reference surface
            double R = (depthSolid ? geodesy::getOuterSolidRadius() :
                        geodesy::getOuterRadius());
            if (geodesy::isCartesian()) {
                // change vertical to surface
                crd(2) = R;
                // add undulation on the surface
                crd(2) = R + sem.getTotalUndulation(crd) - crdIn(2);
            } else {
                double theta = geodesy::sz2rtheta(crd, true, 0, 2, 2, 0)(0);
                crd(0) = R * sin(theta);
                crd(2) = R * cos(theta);
                // add undulation on the surface
                double r = R + sem.getTotalUndulation(crd) - crdIn(2);
                crd(0) = r * sin(theta);
                crd(2) = r * cos(theta);
            }
        }
        // undulated to reference
        crd = sem.undulatedToReference(crd);
    }
    
    // return
    return crd;
}

// compute rotation matrix Q from input to (z, s, phi)
const eigen::DMat33 &Source::computeQzsp(const eigen::DRow3 &spz,
                                         bool ellipticity) const {
    static eigen::DMat33 Qzsp;
    Qzsp.setIdentity();
    
    // axial source needs no rotation
    if (mAxial) {
        return Qzsp;
    }
    
    // first rotation: around vertical axis by back azimuth
    if (!mSourceCentered) {
        // compute back azimuth
        eigen::DRow3 llr = mCrdIn;
        if (mUseDepth) {
            double R = (mDepthSolid ? geodesy::getOuterSolidRadius() :
                        geodesy::getOuterRadius());
            llr(2) = R - llr(2);
        }
        double baz = geodesy::backAzimuth(llr, ellipticity)(0);
        // 1 0        0
        // 0 cos(baz) -sin(baz)
        // 0 sin(baz) cos(baz)
        Qzsp(1, 1) = cos(baz);
        Qzsp(2, 2) = cos(baz);
        Qzsp(1, 2) = -sin(baz);
        Qzsp(2, 1) = sin(baz);
    }
    
    // second rotation: around transpose axis by theta
    if (!geodesy::isCartesian()) {
        // cos(theta) -sin(theta) 0
        // sin(theta) cos(theta)  0
        // 0          0           1
        static eigen::DMat33 Qtheta;
        Qtheta.setIdentity();
        double theta = geodesy::sz2rtheta(spz, true, 0, 2, 2, 0)(0);
        Qtheta(0, 0) = cos(theta);
        Qtheta(1, 1) = cos(theta);
        Qtheta(0, 1) = -sin(theta);
        Qtheta(1, 0) = sin(theta);
        Qzsp.applyOnTheRight(Qtheta);
    }
    return Qzsp;
}

// verbose
std::string Source::verbose(int sindex, const STF &stf,
                            const Mechanism &mechanism) const {
    // title
    using namespace bstring;
    std::stringstream ss;
    ss << boxSubTitle(0, mName + " ", '~');
    
    // central
    bool centralSource =
    sindex == 0 && inparam::gInparamModel.get<std::string>
    ("geodesy:lat_lon_north_pole_mesh") == "SOURCE";
    if (centralSource) {
        ss << "  * This source is used to determine the reference\n";
        ss << "    geographic location on the mesh axis.\n";
    }
    
    // width
    int width = stf.verboseKeyWidth();
    
    // location
    ss << boxSubTitle(2, "Location");
    // horizontal
    if (mAxial && !centralSource) {
        ss << boxEquals(4, width, "horizontal", "on axis");
    } else {
        if (!mSourceCentered) {
            ss << boxEquals(4, width, "latitude", mCrdIn(0));
            ss << boxEquals(4, width, "longitude", mCrdIn(1));
        } else {
            ss << boxEquals(4, width, "distance", mCrdIn(0));
            ss << boxEquals(4, width, "azimuth", mCrdIn(1));
        }
    }
    // vertical
    if (mUseDepth) {
        if (mDepthSolid) {
            ss << boxEquals(4, width, "depth",
                            toString(mCrdIn(2)) + " (below solid surface)");
        } else {
            ss << boxEquals(4, width, "depth",
                            toString(mCrdIn(2)) + " (below outer surface)");
        }
    } else {
        ss << boxEquals(4, width, "radius", mCrdIn(2));
    }
    if (!mAxial && !mSourceCentered) {
        ss << boxEquals(4, width, "ellipticity correction", mEllipticity);
    }
    ss << boxEquals(4, width, "radial geometry", mUndulatedGeometry ?
                    "undulated" : "reference");
    
    // mechanism
    ss << mechanism.verbose(width);
    
    // STF
    ss << stf.verbose();
    return ss.str();
}

// release source to domain
void Source::release(const SE_Model &sem, Domain &domain, double dt,
                     double &minT0) {
    // source count
    int sourceCount = inparam::gInparamSource.get<int>("list_of_sources:[?]");
    std::vector<std::shared_ptr<const Source>> sources;
    sources.reserve(sourceCount);
    
    // first locate source in mesh
    timer::gPreloopTimer.begin("Locating sources in mesh");
    // source index -> quad tag (store because locating is expensive)
    std::map<int, int> sourceIndexQuad;
    std::map<int, int> sourceIndexRank;
    for (int sindex = 0; sindex < sourceCount; sindex++) {
        // create source
        std::shared_ptr<const Source> source = buildInparam(sindex);
        sources.push_back(source);
        
        // compute spz
        const eigen::DRow3 &spz =
        computeSPZ(sem, source->mCrdIn,
                   source->mSourceCentered, false, source->mEllipticity,
                   source->mUseDepth, source->mDepthSolid,
                   source->mUndulatedGeometry,
                   "Source name: " + source->mName, source->mAxial);
        
        // get solid/fluid from mechanism
        std::shared_ptr<const Mechanism> mechanism =
        Mechanism::buildInparam(sindex, source->mName);
        bool inFluid = mechanism->inFluid();
        
        // locate in mesh
        int quadTag = sem.locateInplane(spz({0, 2}).transpose(), inFluid);
        if (quadTag != -1) {
            sourceIndexRank.insert({sindex, mpi::rank()});
            sourceIndexQuad.insert({sindex, quadTag});
        }
    }
    timer::gPreloopTimer.ended("Locating sources in mesh");
    
    // mpi assemble
    // use minimum rank if source is located on more than one ranks
    timer::gPreloopTimer.begin("Assembling and verifying source locations");
    mpi::aggregate(sourceIndexRank, MPI_ALL, MPI_MIN);
    
    // check if all sources are located
    if (sourceIndexRank.size() != sourceCount && mpi::root()) {
        for (int sindex = 0; sindex < sourceCount; sindex++) {
            if (sourceIndexRank.find(sindex) == sourceIndexRank.end()) {
                throw std::runtime_error
                ("Source::release || Error locating source in mesh. || "
                 "Source index = " + bstring::toString(sindex) + " || "
                 "Source name  = " + sources[sindex]->mName);
            }
        }
    }
    timer::gPreloopTimer.ended("Assembling and verifying source locations");
    
    // release
    timer::gPreloopTimer.begin("Computing source patterns and releasing");
    // minimum t0
    minT0 = std::numeric_limits<double>::max();
    // verbose on rank
    std::map<int, std::string> verboseRank;
    // get quads
    const std::vector<Quad> &quads = sem.getQuads();
    // release
    for (int sindex = 0; sindex < sourceCount; sindex++) {
        if (sourceIndexRank.at(sindex) != mpi::rank()) {
            // handled by another rank
            continue;
        }
        
        // compute spz
        const std::shared_ptr<const Source> &source = sources[sindex];
        const eigen::DRow3 &spz =
        computeSPZ(sem, source->mCrdIn,
                   source->mSourceCentered, false, source->mEllipticity,
                   source->mUseDepth, source->mDepthSolid,
                   source->mUndulatedGeometry,
                   "Source name: " + source->mName, source->mAxial);
        
        // inplane interpolation
        int quadTag = sourceIndexQuad.at(sindex);
        const eigen::DRowN &inplaneFactor =
        sem.computeInplaneFactor(spz({0, 2}).transpose(), quadTag);
        
        // STF
        std::unique_ptr<STF> stf = STF::buildInparam(sindex, source->mName, dt);
        // update minimum t0
        minT0 = std::min(minT0, stf->getStartTime());
        
        // mechanism
        std::shared_ptr<const Mechanism> mechanism =
        Mechanism::buildInparam(sindex, source->mName);
        
        // verbose before release
        if (io::gVerbose != io::VerboseLevel::None) {
            verboseRank.insert({sindex, source->
                verbose(sindex, *stf, *mechanism)});
        }
        
        // release
        mechanism->release(source->computeQzsp(spz, source->mEllipticity),
                           source->mAxial, inplaneFactor, spz(1),
                           quads[quadTag], stf, domain);
    }
    // minimum t0 over ranks
    minT0 = mpi::min(minT0);
    if (sourceCount == 0) {
        // no source, use 0.
        minT0 = 0.;
    }
    timer::gPreloopTimer.ended("Computing source patterns and releasing");
    
    // verbose
    if (io::gVerbose == io::VerboseLevel::None) {
        return;
    }
    
    timer::gPreloopTimer.begin("Verbosing sources");
    // gather verbose
    std::vector<std::map<int, std::string>> verboseAll;
    mpi::gather(verboseRank, verboseAll, 0);
    
    // verbose on root
    std::stringstream ss;
    if (mpi::root()) {
        // flatten station-rank map
        std::map<int, std::string> flat;
        for (int iproc = 0; iproc < mpi::nproc(); iproc++) {
            flat.insert(verboseAll[iproc].begin(), verboseAll[iproc].end());
        }
        
        // verbose
        ss << bstring::boxTitle("Sources");
        for (int sindex = 0; sindex < sourceCount; sindex++) {
            ss << flat.at(sindex);
        }
        if (sourceCount == 0) {
            ss << "* No sources in this simulation.\n";
        }
        ss << bstring::boxBaseline() << "\n\n";
    }
    io::cout << ss.str();
    timer::gPreloopTimer.ended("Verbosing sources");
}
