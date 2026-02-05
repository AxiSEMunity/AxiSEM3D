//
//  Quad.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/29/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  Quadrilateral
//  generator of Element in core

#include "Quad.hpp"

// mesh
#include "ExodusMesh.hpp"
#include "ABC.hpp"
#include "LocalMesh.hpp"
#include "NrField.hpp"

// mapping
#include "MappingLinear.hpp"
#include "MappingSpherical.hpp"
#include "MappingSemiSpherical.hpp"

// GLL
#include "GLLPoint.hpp"
#include "vicinity.hpp"

// dt
#include "geodesy.hpp"

// release
#include "Domain.hpp"
#include "SolidElement.hpp"
#include "FluidElement.hpp"

// constructor
Quad::Quad(const ExodusMesh &exodusMesh, const LocalMesh &localMesh,
           int localTag, bool useLuckyNumbers):
mLocalTag(localTag), mGlobalTag(localMesh.mL2G_Element[localTag]),
mFluid(localMesh.mIsElementFluid(localTag)) {
    // model boundary
    mEdgesOnBoundary.insert({"LEFT", exodusMesh.getLeftSide(mGlobalTag)});
    mEdgesOnBoundary.insert({"RIGHT", exodusMesh.getRightSide(mGlobalTag)});
    mEdgesOnBoundary.insert({"BOTTOM", exodusMesh.getBottomSide(mGlobalTag)});
    mEdgesOnBoundary.insert({"TOP", exodusMesh.getTopSide(mGlobalTag)});
    mEdgesOnBoundary.insert({"SOLID_FLUID",
        exodusMesh.getSide("solid_fluid_boundary", mGlobalTag)});
    
    // Nr field
    static eigen::DRow4 nodalNr;
    for (int ivt = 0; ivt < 4; ivt++) {
        int inode = localMesh.mConnectivity(mLocalTag, ivt);
        nodalNr(ivt) = localMesh.mNodalNr(inode);
    }
    mPointNr = std::make_unique<eigen::IRowN>
    (spectrals::interpolateGLL(nodalNr, axial()).array().round().cast<int>());
    // round to lucky numbers
    if (useLuckyNumbers) {
        (*mPointNr) = (*mPointNr).unaryExpr([](int nr) {
            return NrField::nextLuckyNumber(nr);
        });
    }
    
    ////////////// components //////////////
    // 1) mapping
    // form nodal (s,z)
    static eigen::DMat24 nodalSZ;
    for (int ivt = 0; ivt < 4; ivt++) {
        int inode = localMesh.mConnectivity(mLocalTag, ivt);
        nodalSZ.col(ivt) = localMesh.mNodalCoords.row(inode).transpose();
    }
    // shape type
    int gtype = localMesh.mGeometryType(mLocalTag);
    if (gtype < .5) {
        // gtype = 0.0, spherical
        mMapping = std::make_unique<MappingSpherical>(nodalSZ);
    } else if (gtype < 1.5) {
        // gtype = 1.0, linear
        mMapping = std::make_unique<MappingLinear>(nodalSZ);
    } else {
        // gtype = 2.0, semi-spherical
        mMapping = std::make_unique<MappingSemiSpherical>(nodalSZ);
    }
    
    // 2) material
    mMaterial = std::make_unique<Material>(exodusMesh, getNodalSZ(), axial());
    
    // 3) undulation
    mUndulation = std::make_unique<Undulation>();
    
    // 4) ocean load
    mOceanLoad = std::make_unique<OceanLoad>();
}

// setup GLL
void Quad::setupGLL(const ABC &abc, const LocalMesh &localMesh,
                    std::vector<GLLPoint> &GLLPoints) const {
    // compute mass
    static eigen::DMat2N sz;
    static std::array<eigen::DMat22, spectral::nPEM> J;
    const eigen::DRowN &ifact = computeIntegralFactor(sz, J);
    const eigen::arN_DColX &J_PRT = mUndulation->getMassJacobian(sz);
    const eigen::arN_DColX &mass = mMaterial->getMass(ifact, J_PRT, mFluid);
    
    // setup tags, nr, coords and mass
    for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
        // target GLL
        int igll = localMesh.mElementGLL(mLocalTag, ipnt);
        // setup
        GLLPoints[igll].setup(localMesh.mL2G_GLL(igll), // global tag
                              (*mPointNr)(ipnt), // nr
                              sz.col(ipnt), // sz
                              mMapping->getMinEdgeLength() / 1000.); // tol
        // add mass
        GLLPoints[igll].addMass(mass[ipnt], mFluid);
    }
    
    ////////////////////////////////////
    //////////// boundaries ////////////
    ////////////////////////////////////
    
    // solid-fluid
    if (mEdgesOnBoundary.at("SOLID_FLUID") != -1) {
        // compute normals
        static std::vector<int> ipnts;
        static eigen::arP_DMatX3 nSF;
        computeNormal(mEdgesOnBoundary.at("SOLID_FLUID"), sz, J, ipnts, nSF);
        // add normals to points
        for (int ip = 0; ip < spectral::nPED; ip++) {
            int igll = localMesh.mElementGLL(mLocalTag, ipnts[ip]);
            // divid by 2 because fluid and solid elements will add normal twice
            if (mFluid) {
                GLLPoints[igll].addNormalSF(nSF[ip] / 2.);
            } else {
                // negative because normal must point from fluid to solid
                GLLPoints[igll].addNormalSF(-nSF[ip] / 2.);
            }
        }
    }
    
    // clayton ABC
    if (abc.clayton()) {
        for (const std::string &key: abc.getBoundaryKeys()) {
            if (mEdgesOnBoundary.at(key) == -1) {
                continue;
            }
            // compute normals
            static std::vector<int> ipnts;
            static eigen::arP_DMatX3 nABC;
            computeNormal(mEdgesOnBoundary.at(key), sz, J, ipnts, nABC);
            // get properties from material
            eigen::arN_DColX rho, vp, vs;
            mMaterial->getPointwiseRhoVpVs(rho, vp, vs);
            // add ABCs to points
            for (int ip = 0; ip < spectral::nPED; ip++) {
                int ipnt = ipnts[ip];
                int igll = localMesh.mElementGLL(mLocalTag, ipnt);
                GLLPoints[igll].addClaytonABC(key, mFluid, nABC[ip],
                                              rho[ipnt], vp[ipnt], vs[ipnt]);
            }
        }
    }
    
    // sponge ABC
    if (abc.sponge()) {
        // get material
        eigen::arN_DColX rho, vp, vs;
        mMaterial->getPointwiseRhoVpVs(rho, vp, vs);
        
        // loop over gll
        for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
            // regularize 1D/3D
            op1D_3D::regularize1D<eigen::DColX>({
                std::ref(rho[ipnt]),
                std::ref(vp[ipnt]),
                std::ref(vs[ipnt])});
            // compute maximum gamma among all boundaries
            eigen::DColX gammaMax = eigen::DColX::Zero((*mPointNr)(ipnt), 1);
            double distToOuter_min = 1.;
            for (const std::string &key: abc.getBoundaryKeys()) {
                /////////// pattern ///////////
                // geometry
                const auto &outerSpan = abc.getSpongeOuterSpan(key);
                double outer = std::get<0>(outerSpan);
                double span = std::get<1>(outerSpan);
                // coord of me
                double coord = 0.;
                double r = 0.;
                if (geodesy::isCartesian()) {
                    coord = (key == "RIGHT" ? sz(0, ipnt) : sz(1, ipnt));
                    r = sz(1, ipnt);
                } else {
                    const eigen::DCol2 &rt =
                    geodesy::sz2rtheta(sz.col(ipnt).eval(), false);
                    coord = (key == "RIGHT" ? rt(1) : rt(0));
                    r = rt(0);
                }
                // gamma
                double distToOuter = 1. / span * (outer - coord);
                if (distToOuter > distToOuter_min) {
                    // point is inside the inner boundary, skip;
                    // point is closer to other boundary, skip;
                    // there is no need to check distToOuter < 0.
                    continue;
                }
                distToOuter_min = distToOuter;
                
                static const double piHalf = numerical::dPi / 2.;
                double pattern = pow(cos(piHalf * distToOuter), 2.);
                
                /////////// gamma ///////////
                // for theta, change span to arc-length
                if (!geodesy::isCartesian() && key == "RIGHT") {
                    span *= r;
                }
                if (!mFluid) {
                    gammaMax = abc.getU0Solid(std::abs(span),
                                              vp[ipnt], vs[ipnt], rho[ipnt]);
                } else {
                    gammaMax = abc.getU0Fluid(std::abs(span),
                                              vp[ipnt], rho[ipnt]);
                }
                gammaMax *= pattern;
            }
            // release
            if (gammaMax.norm() > numerical::dEpsilon) {
                int igll = localMesh.mElementGLL(mLocalTag, ipnt);
                GLLPoints[igll].addGamma(gammaMax);
            }
        }
    }
    
    // ocean load
    // get data from OceanLoad
    if (mEdgesOnBoundary.at("TOP") != -1 && (*mOceanLoad)) {
        // compute normals
        static std::vector<int> ipnts;
        static eigen::arP_DMatX3 nTop;
        computeNormal(mEdgesOnBoundary.at("TOP"), sz, J, ipnts, nTop);
        // add ocean loads to points
        const eigen::arP_DColX &sumRhoDepth = mOceanLoad->getPointwise();
        for (int ip = 0; ip < spectral::nPED; ip++) {
            int igll = localMesh.mElementGLL(mLocalTag, ipnts[ip]);
            GLLPoints[igll].addOceanLoad(nTop[ip], sumRhoDepth[ip]);
        }
    }
    
    // axial
    if (axial()) {
        const std::vector<int> &ipnts =
        vicinity::constants::gEdgeIPnt[mEdgesOnBoundary.at("LEFT")];
        for (int ipnt: ipnts) {
            int igll = localMesh.mElementGLL(mLocalTag, ipnt);
            GLLPoints[igll].setAxial();
        }
    }
    
    // surface
    if (mEdgesOnBoundary.at("TOP") != -1) {
        const std::vector<int> &ipnts =
        vicinity::constants::gEdgeIPnt[mEdgesOnBoundary.at("TOP")];
        for (int ipnt: ipnts) {
            int igll = localMesh.mElementGLL(mLocalTag, ipnt);
            GLLPoints[igll].setSurface();
        }
    }
}

// compute dt
double Quad::computeDt(double courant, const ABC &abc) const {
    // get vp
    const eigen::DColX &vmaxNr =
    mMaterial->getMaxVelocity().rowwise().maxCoeff();
    
    // 1D coords in reference configuration
    const eigen::DMat2N &szRef = getPointSZ();
    
    // compute coords on slices
    std::vector<eigen::DMat2N> szNr;
    const eigen::DMatXN &dZ = mUndulation->getElemental();
    // both 1D and 3D undulation
    if (geodesy::isCartesian()) {
        for (int nr = 0; nr < dZ.rows(); nr++) {
            eigen::DMat2N sz = szRef;
            sz.row(1) += dZ.row(nr);
            szNr.push_back(sz);
        }
    } else {
        const eigen::DMat2N &rt = geodesy::sz2rtheta(szRef, false);
        const eigen::DRowN &sint = rt.row(1).array().sin();
        const eigen::DRowN &cost = rt.row(1).array().cos();
        for (int nr = 0; nr < dZ.rows(); nr++) {
            eigen::DMat2N sz = szRef;
            sz.row(0) += dZ.row(nr).cwiseProduct(sint);
            sz.row(1) += dZ.row(nr).cwiseProduct(cost);
            szNr.push_back(sz);
        }
    }
    
    // compute hmin and dt slice-wise
    double dt = std::numeric_limits<double>::max();
    int nrMax = std::max((int)vmaxNr.rows(), (int)szNr.size());
    for (int nr = 0; nr < nrMax; nr++) {
        // vmax
        double vmax = (vmaxNr.rows() == 1) ? vmaxNr(0) : vmaxNr(nr);
        
        // hmin
        const eigen::DMat2N &sz = (szNr.size() == 1) ? szNr[0] : szNr[nr];
        double hmin = std::numeric_limits<double>::max();
        for (int ipnt0 = 0; ipnt0 < spectral::nPEM - 1; ipnt0++) {
            int ipol0 = ipnt0 / spectral::nPED;
            int jpol0 = ipnt0 % spectral::nPED;
            for (int ipnt1 = ipnt0 + 1; ipnt1 < spectral::nPEM; ipnt1++) {
                int ipol1 = ipnt1 / spectral::nPED;
                int jpol1 = ipnt1 % spectral::nPED;
                // only consider neighbouring points
                if (std::abs(ipol1 - ipol0) <= 1 &&
                    std::abs(jpol1 - jpol0) <= 1) {
                    hmin = std::min(hmin, (sz.col(ipnt1) -
                                           sz.col(ipnt0)).norm());
                }
            }
        }
        // dt
        dt = std::min(dt, hmin / vmax);
    }
    
    // solid-fluid and clayton BCs are numerically sensitive
    double factorDtForBC = 1.;
    // solid-fluid
    if (mEdgesOnBoundary.at("SOLID_FLUID") != -1) {
        factorDtForBC = .9;
    }
    // clayton ABC
    if (abc.clayton()) {
        for (const std::string &key: abc.getBoundaryKeys()) {
            if (mEdgesOnBoundary.at(key) != -1) {
                factorDtForBC = .5;
                break;
            }
        }
    }
    // decrease DT
    dt *= factorDtForBC;
    
    // courant
    return dt * courant;
}

// release to domain
void Quad::release(const LocalMesh &localMesh,
                   const std::vector<GLLPoint> &GLLPoints,
                   const std::unique_ptr<const AttBuilder> &attBuilder,
                   Domain &domain) {
    // gradient-quadrature operator
    static eigen::DMat2N sz;
    static eigen::DMatPP_RM ifPP;
    std::unique_ptr<const GradientQuadrature<numerical::Real>> grad =
    createGradient<numerical::Real>(sz, ifPP);
    // particle relabelling transformation
    std::unique_ptr<const PRT> prt = mUndulation->createPRT(sz);
    // material and points
    if (mFluid) {
        // acoustic
        std::unique_ptr<const Acoustic> acoustic = mMaterial->createAcoustic();
        // point array
        std::array<std::shared_ptr<FluidPoint>, spectral::nPEM> points;
        for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
            int igll = localMesh.mElementGLL(mLocalTag, ipnt);
            points[ipnt] = GLLPoints[igll].getFluidPoint();
        }
        // element
        mFluidElement = std::make_shared<FluidElement>
        (mGlobalTag, grad, prt, acoustic, points);
        domain.addFluidElement(mFluidElement);
    } else {
        // elastic
        std::unique_ptr<const Elastic> elastic =
        mMaterial->createElastic(attBuilder, computeWeightsCG4(ifPP));
        // point array
        std::array<std::shared_ptr<SolidPoint>, spectral::nPEM> points;
        for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
            int igll = localMesh.mElementGLL(mLocalTag, ipnt);
            points[ipnt] = GLLPoints[igll].getSolidPoint();
        }
        // element
        mSolidElement = std::make_shared<SolidElement>
        (mGlobalTag, grad, prt, elastic, points);
        domain.addSolidElement(mSolidElement);
    }
    
    // free dummy memory
    mPointNr.reset();
    mMaterial.reset();
    mUndulation.reset();
    mOceanLoad.reset();
}


//////////////////////// interal ////////////////////////
// compute integral factor
eigen::DRowN Quad::
computeIntegralFactor(eigen::DMat2N &sz, std::array<eigen::DMat22,
                      spectral::nPEM> &J) const {
    // xieta and weights
    const eigen::DMat2N &xieta = spectrals::getXiEtaElement(axial());
    const eigen::DRowN &weights = spectrals::getWeightsElement(axial());
    
    // compute integral factor
    eigen::DRowN ifact;
    for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
        // element Jacobian
        J[ipnt] = mMapping->jacobian(xieta.col(ipnt));
        // integral weights and Jacobian
        ifact(ipnt) = weights(ipnt) * J[ipnt].determinant();
        // s in integral
        sz.col(ipnt) = mMapping->mapping(xieta.col(ipnt));
        double s = sz(0, ipnt);
        if (axial()) {
            // axial element
            if (ipnt / spectral::nPED == 0) { // ipol == 0
                // axial point, L'Hospital's Rule
                ifact(ipnt) *= J[ipnt](0, 0);
            } else {
                // non-axial point
                ifact(ipnt) *= s / (1. + xieta(0, ipnt));
            }
        } else {
            // non-axial element
            ifact(ipnt) *= s;
        }
    }
    return ifact;
}

// get normal
void Quad::computeNormal(int edge, const eigen::DMat2N &sz,
                         const std::array<eigen::DMat22, spectral::nPEM> &J,
                         std::vector<int> &ipnts,
                         eigen::arP_DMatX3 &normal) const {
    // xieta
    const eigen::DMat2N &xieta = spectrals::getXiEtaElement(axial());
    
    // points on edge
    ipnts = vicinity::constants::gEdgeIPnt[edge];
    
    // normal
    for (int ip = 0; ip < spectral::nPED; ip++) {
        int ipnt = ipnts[ip];
        // 1D normal
        eigen::DCol2 n1D = mMapping->normal(edge, J[ipnt]);
        // integral weights
        int ipol = ipnt / spectral::nPED;
        int jpol = ipnt % spectral::nPED;
        if (axial()) {
            // edge must be even
            if (edge % 2 != 0) {
                throw std::runtime_error("Quad::computeNormal || Impossible.");
            }
            n1D *= spectrals::gWeightsGLJ(ipol);
        } else {
            n1D *= spectrals::gWeightsGLL(edge % 2 == 0 ? ipol : jpol);
        }
        // s in integral
        double s = sz(0, ipnt);
        if (axial()) {
            // axial element
            if (ipol == 0) {
                // axial point, L'Hospital's Rule
                n1D *= J[ipnt](0, 0);
            } else {
                // non-axial point
                n1D *= s / (1. + xieta(0, ipnt));
            }
        } else {
            // non-axial element
            n1D *= s;
        }
        // 3D normal with undulation
        // must use ip instead of varpol
        normal[ip] = mUndulation->computeNormal3D(n1D, sz, ipnt);
    }
}

// weights for CG4 attenuation
eigen::DRow4 Quad::computeWeightsCG4(const eigen::DMatPP_RM &ifPP) const {
    if (spectral::nPol != 4) {
        return eigen::DRow4::Zero();
    }
    // weights on CG4 points (marked O)
    // x x x x x
    // x O x O x
    // x x x x x
    // x O x O x
    // x x x x x
    eigen::DRow4 wcg4;
    wcg4(0) = (ifPP(0, 0) + ifPP(0, 1) + ifPP(1, 0) + ifPP(1, 1) +
               0.50 * (ifPP(0, 2) + ifPP(1, 2) + ifPP(2, 0) + ifPP(2, 1)) +
               0.25 * ifPP(2, 2)) / ifPP(1, 1);
    wcg4(1) = (ifPP(0, 3) + ifPP(0, 4) + ifPP(1, 3) + ifPP(1, 4) +
               0.50 * (ifPP(0, 2) + ifPP(1, 2) + ifPP(2, 3) + ifPP(2, 4)) +
               0.25 * ifPP(2, 2)) / ifPP(1, 3);
    wcg4(2) = (ifPP(3, 0) + ifPP(3, 1) + ifPP(4, 0) + ifPP(4, 1) +
               0.50 * (ifPP(2, 0) + ifPP(2, 1) + ifPP(3, 2) + ifPP(4, 2)) +
               0.25 * ifPP(2, 2)) / ifPP(3, 1);
    wcg4(3) = (ifPP(3, 3) + ifPP(3, 4) + ifPP(4, 3) + ifPP(4, 4) +
               0.50 * (ifPP(2, 3) + ifPP(2, 4) + ifPP(3, 2) + ifPP(4, 2)) +
               0.25 * ifPP(2, 2)) / ifPP(3, 3);
    return wcg4;
}

// get element
std::shared_ptr<const Element> Quad::getElement() const {
    if (mSolidElement && !mFluidElement) {
        return mSolidElement;
    }
    if (mFluidElement && !mSolidElement) {
        return mFluidElement;
    }
    throw std::runtime_error("Quad::getElement || "
                             "Element has not been created and released.");
}

// get solid element
const std::shared_ptr<SolidElement> &Quad::getSolidElement() const {
    if (mSolidElement) {
        return mSolidElement;
    }
    throw std::runtime_error("Quad::getSolidElement || "
                             "Element is not in solid.");
}

// get fluid element
const std::shared_ptr<FluidElement> &Quad::getFluidElement() const {
    if (mFluidElement) {
        return mFluidElement;
    }
    throw std::runtime_error("Quad::getFluidElement || "
                             "Element is not in fluid.");
}
