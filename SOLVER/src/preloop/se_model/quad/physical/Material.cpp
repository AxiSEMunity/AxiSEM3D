//
//  Material.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/29/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  material properties
//  generator of Acoustic and Elastic in core

#include "Material.hpp"

// constructor
#include "ExodusMesh.hpp"
#include "geodesy.hpp"
#include "vector_tools.hpp"

// release
#include "Isotropic.hpp"
#include "TransverselyIsotropic.hpp"
#include "Anisotropic.hpp"
#include "AttBuilder.hpp"

////////////////////////// input //////////////////////////
// constructor
Material::Material(const ExodusMesh &exodusMesh, const eigen::DMat24 &nodalSZ,
                   bool axial) {
    // get radial coords and variables from mesh
    const auto &radialCoords = exodusMesh.getRadialCoords();
    const auto &radialVariables = exodusMesh.getRadialVariables();
    
    // radial coords of me
    const eigen::DRow4 &r = geodesy::isCartesian() ?
    nodalSZ.row(1).eval() : nodalSZ.colwise().norm().eval();
    
    // locate center in radial profile
    int index0 = -1, index1 = -1;
    double factor0 = 0., factor1 = 0.;
    vector_tools::linearInterpSorted(radialCoords, r.mean(),
                                     index0, index1, factor0, factor1);
    double distTol = exodusMesh.getGlobalVariable("dist_tolerance");
    double r0 = radialCoords[index0] - distTol;
    double r1 = radialCoords[index1] + distTol;
    
    // loop over all variables
    for (auto it = radialVariables.begin(); it != radialVariables.end(); ++it) {
        // end points
        const eigen::DColX &rvals = it->second;
        double v0 = rvals(index0);
        double v1 = rvals(index1);
        // interpolation
        const eigen::DRow4 &nvals =
        (v1 - v0) / (r1 - r0) * (r.array() - r0) + v0;
        // add to map
        mProperties.insert({it->first, NodalPhysicalProperty(nvals, axial)});
    }
    
    // reduce TISO to ISO if possible
    if (currentAnisotropy() == AnisotropyType::TISO) {
        const eigen::DRow4 &vpv = mProperties.at("VPV").getNodalData();
        const eigen::DRow4 &vph = mProperties.at("VPH").getNodalData();
        const eigen::DRow4 &vsv = mProperties.at("VSV").getNodalData();
        const eigen::DRow4 &vsh = mProperties.at("VSH").getNodalData();
        const eigen::DRow4 &eta = mProperties.at("ETA").getNodalData();
        // TISO to ISO
        if ((vpv - vph).norm() < vpv.norm() * numerical::dEpsilon &&
            (vsv - vsh).norm() < vsv.norm() * numerical::dEpsilon &&
            (eta - eigen::DRow4::Ones()).norm() < 2. * numerical::dEpsilon) {
            // add ISO
            mProperties.insert({"VP", NodalPhysicalProperty(vpv, axial)});
            mProperties.insert({"VS", NodalPhysicalProperty(vsv, axial)});
            // remove TISO
            mProperties.erase("VPV");
            mProperties.erase("VPH");
            mProperties.erase("VSV");
            mProperties.erase("VSH");
            mProperties.erase("ETA");
        }
    }
}

// add a 3D property
void Material::addProperty3D(const std::string &propKey,
                             const Volumetric3D::ReferenceKind &refKind,
                             const eigen::arN_IColX &inScope,
                             const eigen::arN_DColX &propValue) {
    // evolve anisotropy
    // ISO -> TISO
    if (currentAnisotropy() == AnisotropyType::ISO &&
        (propKey == "VPV" || propKey == "VPH" ||
         propKey == "VSV" || propKey == "VSH" || propKey == "ETA")) {
        evolveISO_TISO();
    }
    // ISO -> TISO -> ANISO
    if (propKey.front() == 'C') {
        if (currentAnisotropy() == AnisotropyType::ISO) {
            evolveISO_TISO();
        }
        if (currentAnisotropy() == AnisotropyType::TISO) {
            evolveTISO_ANISO();
        }
    }
    
    // backward compatibility
    // TISO
    if (currentAnisotropy() == AnisotropyType::TISO &&
        (propKey == "VP" || propKey == "VS")) {
        // recursive call
        addProperty3D(propKey + "V", refKind, inScope, propValue);
        addProperty3D(propKey + "H", refKind, inScope, propValue);
    }
    // ANISO
    if (currentAnisotropy() == AnisotropyType::ANISO &&
        propKey.front() == 'V') {
        // give a clearer instruction
        throw std::runtime_error("Material::getProperty || "
                                 "Setting velocity is prohibited "
                                 "in full anisotropy.");
    }
    
    // compute 3D value
    eigen::arN_DColX valNew;
    if (refKind == Volumetric3D::ReferenceKind::ABS) {
        // absolute
        valNew = propValue;
    } else {
        eigen::arN_DColX val0;
        if (refKind == Volumetric3D::ReferenceKind::REF1D) {
            // 1D as reference
            val0 = getProperty(propKey).getPointwiseNodal();
        } else if (refKind == Volumetric3D::ReferenceKind::REF3D) {
            // 3D as reference
            val0 = getProperty(propKey).getPointwise();
        } else {
            // (3D - 1D) as reference
            val0 = getProperty(propKey).getPointwise();
            const eigen::arN_DColX &val1D =
            getProperty(propKey).getPointwiseNodal();
            for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
                op1D_3D::addTo(-val1D[ipnt], val0[ipnt]);
            }
        }
        // set 3D
        for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
            op1D_3D::times((propValue[ipnt].array() + 1.).matrix(),
                           val0[ipnt], valNew[ipnt]);
        }
    }
    
    // out-of-scope mask
    const eigen::arN_DColX &oriVal = getProperty(propKey).getPointwise();
    for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
        if (oriVal[ipnt].rows() == 1) {
            valNew[ipnt] = inScope[ipnt].select(valNew[ipnt], oriVal[ipnt](0));
        } else {
            valNew[ipnt] = inScope[ipnt].select(valNew[ipnt], oriVal[ipnt]);
        }
    }
    
    // final set
    getProperty(propKey).setGLL(valNew);
}

// finished 3D properties
void Material::finished3D() {
    // create velocity for ANISO using Voigt average (required by ABC)
    if (currentAnisotropy() == AnisotropyType::ANISO) {
        // get density and Cijkl
        const NodalPhysicalProperty &rho = mProperties.at("RHO");
        const NodalPhysicalProperty &C11 = mProperties.at("C11");
        const NodalPhysicalProperty &C12 = mProperties.at("C12");
        const NodalPhysicalProperty &C13 = mProperties.at("C13");
        const NodalPhysicalProperty &C22 = mProperties.at("C22");
        const NodalPhysicalProperty &C23 = mProperties.at("C23");
        const NodalPhysicalProperty &C33 = mProperties.at("C33");
        const NodalPhysicalProperty &C44 = mProperties.at("C44");
        const NodalPhysicalProperty &C55 = mProperties.at("C55");
        const NodalPhysicalProperty &C66 = mProperties.at("C66");
        // kappa, mu in Voigt average
        const NodalPhysicalProperty &kp =
        (C11 + C22 + C33 + (C12 + C23 + C13) * 2.) / 9.;
        const NodalPhysicalProperty &mu =
        (C11 + C22 + C33 - (C12 + C23 + C13) + (C44 + C55 + C66) * 3.) / 15.;
        // vp, vs
        mProperties.insert({"VP", (kp / rho).pow(.5)});
        mProperties.insert({"VS", (mu / rho).pow(.5)});
    }
    
    // try to degenerate TISO to ISO
    if (currentAnisotropy() == AnisotropyType::TISO) {
        eigen::DMatXN vpv = getElemental("VPV");
        eigen::DMatXN vph = getElemental("VPH");
        eigen::DMatXN vsv = getElemental("VSV");
        eigen::DMatXN vsh = getElemental("VSH");
        eigen::DMatXN eta = getElemental("ETA");
        op1D_3D::regularize1D<eigen::DMatXN>({std::ref(vpv), std::ref(vph)});
        op1D_3D::regularize1D<eigen::DMatXN>({std::ref(vsv), std::ref(vsh)});
        // TISO to ISO
        if ((vpv - vph).norm() < vpv.norm() * numerical::dEpsilon &&
            (vsv - vsh).norm() < vsv.norm() * numerical::dEpsilon &&
            (eta.array() - 1.).matrix().norm() <
            eta.norm() * numerical::dEpsilon) {
            // add ISO
            mProperties.insert({"VP", mProperties.at("VPV")});
            mProperties.insert({"VS", mProperties.at("VSV")});
            // remove TISO
            mProperties.erase("VPV");
            mProperties.erase("VPH");
            mProperties.erase("VSV");
            mProperties.erase("VSH");
            mProperties.erase("ETA");
        }
    }
}


////////////////////////// output //////////////////////////
// get maximum velocity for dt
eigen::DMatXN Material::getMaxVelocity() const {
    // three types of anisotropy
    if (currentAnisotropy() == AnisotropyType::ANISO) {
        // get density and Cijkl
        eigen::DMatXN rho = getElemental("RHO");
        eigen::DMatXN C11 = getElemental("C11");
        eigen::DMatXN C12 = getElemental("C12");
        eigen::DMatXN C13 = getElemental("C13");
        eigen::DMatXN C22 = getElemental("C22");
        eigen::DMatXN C23 = getElemental("C23");
        eigen::DMatXN C33 = getElemental("C33");
        op1D_3D::regularize1D<eigen::DMatXN>({std::ref(rho),
            std::ref(C11), std::ref(C12), std::ref(C13),
            std::ref(C22), std::ref(C23), std::ref(C33)});
        // max of C11, C22, C33
        const eigen::DMatXN &major = C11.cwiseMax(C22).cwiseMax(C33);
        // max of C12, C13, C23
        const eigen::DMatXN &minor = C12.cwiseMax(C13).cwiseMax(C23);
        // max kappa
        const eigen::DMatXN &kappa = (major + 2. * minor) / 3.;
        // vp = sqrt(kappa / rho)
        return kappa.cwiseQuotient(rho).cwiseSqrt();
    } else if (currentAnisotropy() == AnisotropyType::TISO) {
        eigen::DMatXN vpv = getElemental("VPV");
        eigen::DMatXN vph = getElemental("VPH");
        op1D_3D::regularize1D<eigen::DMatXN>({std::ref(vpv), std::ref(vph)});
        return vpv.cwiseMax(vph);
    } else {
        return getElemental("VP");
    }
}

// get mass for GLL-point setup
eigen::arN_DColX Material::getMass(const eigen::DRowN &integralFactor,
                                   const eigen::arN_DColX &jacobianPRT,
                                   bool fluid) const {
    // form density
    eigen::arN_DColX mass;
    if (fluid) {
        // get rho and vp
        const eigen::arN_DColX &rho = getPointwise("RHO");
        const eigen::arN_DColX &vp = getPointwise("VP");
        for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
            // kappa = rho * vp ** 2
            op1D_3D::times(rho[ipnt], vp[ipnt].array().square().matrix(),
                           mass[ipnt]);
            // density = 1 / kappa
            mass[ipnt] = mass[ipnt].cwiseInverse();
        }
    } else {
        mass = getPointwise("RHO");
    }
    
    // integrate density to mass
    for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
        mass[ipnt] *= integralFactor(ipnt);
    }
    
    // jacobian of PRT
    if (jacobianPRT[0].rows() > 0) {
        for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
            op1D_3D::times(mass[ipnt], jacobianPRT[ipnt], mass[ipnt]);
        }
    }
    return mass;
}

// create Acoustic
std::unique_ptr<Acoustic> Material::createAcoustic() const {
    // K = 1 / rho
    const eigen::DMatXN &K = getElemental("RHO").cwiseInverse();
    // 1D or 3D
    if (K.rows() == 1) {
        return std::make_unique<Acoustic>(op1D_3D::toPP(K));
    } else {
        return std::make_unique<Acoustic>(K);
    }
}

// create Elastic
std::unique_ptr<Elastic> Material::
createElastic(const std::unique_ptr<const AttBuilder> &attBuilder,
              const eigen::DRow4 &weightsCG4) const {
    // three types of anisotropy
    if (currentAnisotropy() == AnisotropyType::ANISO) {
        return createAnisotropic(attBuilder, weightsCG4);
    } else if (currentAnisotropy() == AnisotropyType::TISO) {
        return createTISO(attBuilder, weightsCG4);
    } else {
        return createIsotropic(attBuilder, weightsCG4);
    }
}

// anisotropic
std::unique_ptr<Elastic> Material::
createAnisotropic(const std::unique_ptr<const AttBuilder> &attBuilder,
                  const eigen::DRow4 &weightsCG4) const {
    using op1D_3D::toPP;
    using std::ref;
    
    // Cijkl
    // C11
    eigen::DMatXN C11 = getElemental("C11");
    eigen::DMatXN C12 = getElemental("C12");
    eigen::DMatXN C13 = getElemental("C13");
    eigen::DMatXN C14 = getElemental("C14");
    eigen::DMatXN C15 = getElemental("C15");
    eigen::DMatXN C16 = getElemental("C16");
    // C22
    eigen::DMatXN C22 = getElemental("C22");
    eigen::DMatXN C23 = getElemental("C23");
    eigen::DMatXN C24 = getElemental("C24");
    eigen::DMatXN C25 = getElemental("C25");
    eigen::DMatXN C26 = getElemental("C26");
    // C33
    eigen::DMatXN C33 = getElemental("C33");
    eigen::DMatXN C34 = getElemental("C34");
    eigen::DMatXN C35 = getElemental("C35");
    eigen::DMatXN C36 = getElemental("C36");
    // C44
    eigen::DMatXN C44 = getElemental("C44");
    eigen::DMatXN C45 = getElemental("C45");
    eigen::DMatXN C46 = getElemental("C46");
    // C55
    eigen::DMatXN C55 = getElemental("C55");
    eigen::DMatXN C56 = getElemental("C56");
    // C66
    eigen::DMatXN C66 = getElemental("C66");
    
    // attenuation
    std::unique_ptr<Attenuation> attenuation = nullptr;
    bool elastic1D = false;
    if (attBuilder) {
        // 1D or 3D
        eigen::DMatXN QKp = getElemental("QKAPPA");
        eigen::DMatXN QMu = getElemental("QMU");
        elastic1D = op1D_3D::regularize1D<eigen::DMatXN>
        ({
            ref(C11), ref(C12), ref(C13), ref(C14), ref(C15), ref(C16),
            ref(C22), ref(C23), ref(C24), ref(C25), ref(C26),
            ref(C33), ref(C34), ref(C35), ref(C36),
            ref(C44), ref(C45), ref(C46),
            ref(C55), ref(C56),
            ref(C66),
            ref(QKp), ref(QMu)});
        
        // build attenuation
        eigen::DMatXN kp =
        (C11 + C22 + C33 + (C12 + C23 + C13) * 2.) / 9.;
        eigen::DMatXN mu =
        (C11 + C22 + C33 - (C12 + C23 + C13) + (C44 + C55 + C66) * 3.) / 15.;
        C11 -= (kp + 4. / 3. * mu);
        C22 -= (kp + 4. / 3. * mu);
        C33 -= (kp + 4. / 3. * mu);
        C12 -= (kp - 2. / 3. * mu);
        C23 -= (kp - 2. / 3. * mu);
        C13 -= (kp - 2. / 3. * mu);
        C44 -= mu;
        C55 -= mu;
        C66 -= mu;
        attenuation = attBuilder->createAttenuation(QKp, QMu, kp, mu,
                                                    weightsCG4, elastic1D);
        C11 += (kp + 4. / 3. * mu);
        C22 += (kp + 4. / 3. * mu);
        C33 += (kp + 4. / 3. * mu);
        C12 += (kp - 2. / 3. * mu);
        C23 += (kp - 2. / 3. * mu);
        C13 += (kp - 2. / 3. * mu);
        C44 += mu;
        C55 += mu;
        C66 += mu;
    } else {
        // 1D or 3D
        elastic1D = op1D_3D::regularize1D<eigen::DMatXN>
        ({
            ref(C11), ref(C12), ref(C13), ref(C14), ref(C15), ref(C16),
            ref(C22), ref(C23), ref(C24), ref(C25), ref(C26),
            ref(C33), ref(C34), ref(C35), ref(C36),
            ref(C44), ref(C45), ref(C46),
            ref(C55), ref(C56),
            ref(C66)});
    }
    
    // elastic
    if (elastic1D) {
        return std::make_unique<Anisotropic>
        (attenuation,
         toPP(C11), toPP(C12), toPP(C13), toPP(C14), toPP(C15), toPP(C16),
         toPP(C22), toPP(C23), toPP(C24), toPP(C25), toPP(C26),
         toPP(C33), toPP(C34), toPP(C35), toPP(C36),
         toPP(C44), toPP(C45), toPP(C46),
         toPP(C55), toPP(C56),
         toPP(C66));
    } else {
        return std::make_unique<Anisotropic>
        (attenuation,
         C11, C12, C13, C14, C15, C16,
         C22, C23, C24, C25, C26,
         C33, C34, C35, C36,
         C44, C45, C46,
         C55, C56,
         C66);
    }
}

// transversely isotropic
std::unique_ptr<Elastic> Material::
createTISO(const std::unique_ptr<const AttBuilder> &attBuilder,
           const eigen::DRow4 &weightsCG4) const {
    using op1D_3D::toPP;
    using std::ref;
    
    // density and velocity
    const eigen::DMatXN &rho = getElemental("RHO");
    const eigen::DMatXN &vpv = getElemental("VPV");
    const eigen::DMatXN &vph = getElemental("VPH");
    const eigen::DMatXN &vsv = getElemental("VSV");
    const eigen::DMatXN &vsh = getElemental("VSH");
    
    // A, C, F, L, N
    eigen::DMatXN A, C, L, N;
    op1D_3D::times(rho, vph.array().square().matrix(), A);
    op1D_3D::times(rho, vpv.array().square().matrix(), C);
    op1D_3D::times(rho, vsv.array().square().matrix(), L);
    op1D_3D::times(rho, vsh.array().square().matrix(), N);
    // F = eta * (A - 2 L)
    eigen::DMatXN eta = getElemental("ETA");
    op1D_3D::regularize1D<eigen::DMatXN>({ref(A), ref(L), ref(eta)});
    eigen::DMatXN F = eta.cwiseProduct(A - 2. * L);
    
    // attenuation
    std::unique_ptr<Attenuation> attenuation = nullptr;
    bool elastic1D = false;
    if (attBuilder) {
        // 1D or 3D
        eigen::DMatXN QKp = getElemental("QKAPPA");
        eigen::DMatXN QMu = getElemental("QMU");
        elastic1D = op1D_3D::regularize1D<eigen::DMatXN>
        ({ref(A), ref(C), ref(F), ref(L), ref(N), ref(QKp), ref(QMu)});
        
        // build attenuation
        eigen::DMatXN kp = (4. * A + C + 4. * F - 4. * N) / 9.;
        eigen::DMatXN mu = (A + C - 2. * F + 6. * L + 5. * N) / 15.;
        A -= (kp + 4. / 3. * mu);
        C -= (kp + 4. / 3. * mu);
        F -= (kp - 2. / 3. * mu);
        L -= mu;
        N -= mu;
        attenuation = attBuilder->createAttenuation(QKp, QMu, kp, mu,
                                                    weightsCG4, elastic1D);
        A += (kp + 4. / 3. * mu);
        C += (kp + 4. / 3. * mu);
        F += (kp - 2. / 3. * mu);
        L += mu;
        N += mu;
    } else {
        // 1D or 3D
        elastic1D = op1D_3D::regularize1D<eigen::DMatXN>
        ({ref(A), ref(C), ref(F), ref(L), ref(N)});
    }
    
    // elastic
    if (elastic1D) {
        return std::make_unique<TransverselyIsotropic>
        (attenuation, toPP(A), toPP(C), toPP(F), toPP(L), toPP(N));
    } else {
        return std::make_unique<TransverselyIsotropic>
        (attenuation, A, C, F, L, N);
    }
}

// isotropic
std::unique_ptr<Elastic> Material::
createIsotropic(const std::unique_ptr<const AttBuilder> &attBuilder,
                const eigen::DRow4 &weightsCG4) const {
    using op1D_3D::toPP;
    using std::ref;
    
    // rho, vp, vs
    const eigen::DMatXN &rho = getElemental("RHO");
    const eigen::DMatXN &vp = getElemental("VP");
    const eigen::DMatXN &vs = getElemental("VS");
    
    // lambda, mu
    eigen::DMatXN lambda, mu;
    op1D_3D::times(rho, vp.array().square().matrix(), lambda);
    op1D_3D::times(rho, vs.array().square().matrix(), mu);
    op1D_3D::addTo(-2. * mu, lambda);
    
    // attenuation
    std::unique_ptr<Attenuation> attenuation = nullptr;
    bool elastic1D = false;
    if (attBuilder) {
        // 1D or 3D
        eigen::DMatXN QKp = getElemental("QKAPPA");
        eigen::DMatXN QMu = getElemental("QMU");
        elastic1D = op1D_3D::regularize1D<eigen::DMatXN>
        ({ref(lambda), ref(mu), ref(QKp), ref(QMu)});
        
        // build attenuation
        eigen::DMatXN kp = lambda + (2. / 3.) * mu;
        attenuation = attBuilder->createAttenuation(QKp, QMu, kp, mu,
                                                    weightsCG4, elastic1D);
        lambda = kp - (2. / 3.) * mu;
    } else {
        // 1D or 3D
        elastic1D = op1D_3D::regularize1D<eigen::DMatXN>
        ({ref(lambda), ref(mu)});
    }
    
    // elastic
    if (elastic1D) {
        return std::make_unique<Isotropic>
        (attenuation, toPP(lambda), toPP(mu));
    } else {
        return std::make_unique<Isotropic>
        (attenuation, lambda, mu);
    }
}
