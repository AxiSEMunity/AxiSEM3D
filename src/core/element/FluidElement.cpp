//
//  FluidElement.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/2/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  fluid spectral element

#include "FluidElement.hpp"
#include "FluidPoint.hpp"
#include "fft.hpp"
// output
#include "mapPPvsN.hpp"
// measure
#include "timer.hpp"

using spectral::nPEM;

// constructor
FluidElement::FluidElement(int quadTag, std::unique_ptr<const GradQuad> &grad,
                           std::unique_ptr<const PRT> &prt,
                           std::unique_ptr<const Acoustic> &acoustic,
                           const std::array<std::shared_ptr<FluidPoint>,
                           spectral::nPEM> &points):
Element(quadTag, grad, prt), mAcoustic(acoustic.release()),
mInFourier((mPRT ? mPRT->is1D() : true) && mAcoustic->is1D()),
mPoints(points) {
    // construct derived
    constructDerived();
}

// copy constructor
FluidElement::FluidElement(const FluidElement &other):
Element(other), mAcoustic(std::make_unique<Acoustic>(*(other.mAcoustic))),
mInFourier((mPRT ? mPRT->is1D() : true) && mAcoustic->is1D()),
mPoints(other.mPoints) {
    // construct derived
    constructDerived();
}

// construct derived
void FluidElement::constructDerived() {
    // point set
    pointSet(mInFourier);
    
    // check compatibility
    mAcoustic->checkCompatibility(mNr, mInFourier);
    
    // workspace
    if (sStrainSpherical_CD.rows() < mNr) {
        expandWorkspace(mNr);
    }
    
    // report request to FFT
    if (!mInFourier) {
        fft::gFFT_N3.addNR(mNr);
    }
}

// type info
std::string FluidElement::typeInfo() const {
    std::string info = "FluidElement";
    if (mPRT) {
        if (mPRT->is1D()) {
            info = info + "$PRT1D";
        } else {
            info = info + "$PRT3D";
        }
    }
    if (mAcoustic->is1D()) {
        info = info + "$Acoustic1D";
    } else {
        info = info + "$Acoustic3D";
    }
    return info;
}


/////////////////////////// point ///////////////////////////
// get point
Point &FluidElement::getPoint(int ipnt) const {
    return *(mPoints[ipnt]);
}


/////////////////////////// time loop ///////////////////////////
// collect displacement from points
void FluidElement::
collectDisplFromPoints(eigen::vec_ar1_CMatPP_RM &displElem) const {
    for (int ipol = 0; ipol < spectral::nPED; ipol++) {
        for (int jpol = 0; jpol < spectral::nPED; jpol++) {
            mPoints[ipol * spectral::nPED + jpol]->
            scatterDisplToElement(displElem, mNu_1, ipol, jpol);
        }
    }
}

// displacement to stiffness
void FluidElement::
displToStiff(const eigen::vec_ar1_CMatPP_RM &displElem,
             eigen::vec_ar1_CMatPP_RM &stiffElem) const {
    // displ => strain
    mGradQuad->computeGrad3(displElem,
                            sStrainSpherical_FR, mNu_1);
    
    // strain => stress
    if (mPRT) {
        mTransform->transformSPZ_RTZ3(sStrainSpherical_FR, mNu_1);
        if (mInFourier) {
            mPRT->sphericalToUndulated3_FR(sStrainSpherical_FR,
                                           sStrainUndulated_FR, mNu_1);
            mAcoustic->strainToStress_FR(sStrainUndulated_FR,
                                         sStressUndulated_FR, mNu_1);
            mPRT->undulatedToSpherical3_FR(sStressUndulated_FR,
                                           sStressSpherical_FR, mNu_1);
        } else {
            fft::gFFT_N3.computeC2R(sStrainSpherical_FR,
                                    sStrainSpherical_CD, mNr);
            mPRT->sphericalToUndulated3_CD(sStrainSpherical_CD,
                                           sStrainUndulated_CD, mNr);
            mAcoustic->strainToStress_CD(sStrainUndulated_CD,
                                         sStressUndulated_CD, mNr);
            mPRT->undulatedToSpherical3_CD(sStressUndulated_CD,
                                           sStressSpherical_CD, mNr);
            fft::gFFT_N3.computeR2C(sStressSpherical_CD,
                                    sStressSpherical_FR, mNr);
        }
        mTransform->transformRTZ_SPZ3(sStressSpherical_FR, mNu_1);
    } else {
        if (mInFourier) {
            mAcoustic->strainToStress_FR(sStrainSpherical_FR,
                                         sStressSpherical_FR, mNu_1);
        } else {
            fft::gFFT_N3.computeC2R(sStrainSpherical_FR,
                                    sStrainSpherical_CD, mNr);
            mAcoustic->strainToStress_CD(sStrainSpherical_CD,
                                         sStressSpherical_CD, mNr);
            fft::gFFT_N3.computeR2C(sStressSpherical_CD,
                                    sStressSpherical_FR, mNr);
        }
    }
    
    // stress => stiffness
    mGradQuad->computeQuad3(sStressSpherical_FR,
                            stiffElem, mNu_1);
}

// add stiffness to points
// allow a derived class to change stiffElem (no const)
void FluidElement::
addStiffToPoints(eigen::vec_ar1_CMatPP_RM &stiffElem) const {
    for (int ipol = 0; ipol < spectral::nPED; ipol++) {
        for (int jpol = 0; jpol < spectral::nPED; jpol++) {
            mPoints[ipol * spectral::nPED + jpol]->
            gatherStiffFromElement(stiffElem, ipol, jpol);
        }
    }
}

// compute stiffness term
void FluidElement::computeStiff() const {
    // collect displacement from points
    collectDisplFromPoints(sDisplSpherical_FR);
    
    // displacement to stiffness
    displToStiff(sDisplSpherical_FR, sStiffSpherical_FR);
    
    // add stiffness to points
    addStiffToPoints(sStiffSpherical_FR);
}

// measure cost
double FluidElement::measure(int count) const {
    // random displacement
    for (int ipnt = 0; ipnt < nPEM; ipnt++) {
        mPoints[ipnt]->randomDispl();
    }
    // measure
    SimpleTimer tm;
    tm.start();
    for (int irep = 0; irep < count; irep++) {
        computeStiff();
    }
    tm.pause();
    // reset to zero
    for (int ipnt = 0; ipnt < nPEM; ipnt++) {
        mPoints[ipnt]->resetToZero();
    }
    // return total, not divided by count
    return tm.elapsedTotal();
}


/////////////////////////// source ///////////////////////////
// prepare pressure source
void FluidElement::preparePressureSource() const {
    for (int ipnt = 0; ipnt < nPEM; ipnt++) {
        mPoints[ipnt]->preparePressureSource();
    }
}

// add pressure source
void FluidElement::addPressureSource(const eigen::CMatXN &pressure,
                                     int nu_1_pressure) const {
    for (int ipnt = 0; ipnt < nPEM; ipnt++) {
        mPoints[ipnt]->addPressureSource(pressure, nu_1_pressure, ipnt);
    }
}


/////////////////////////// wavefield output ///////////////////////////
// prepare wavefield output
void FluidElement::
prepareWavefieldOutput(const channel::fluid::ChannelOptions &chops,
                       bool enforceCoordTransform) {
    // pressure
    if (chops.mNeedBufferP) {
        for (int ipnt = 0; ipnt < nPEM; ipnt++) {
            mPoints[ipnt]->preparePressureOutput();
        }
    }
    
    // delta
    if (chops.mNeedBufferD) {
        for (int ipnt = 0; ipnt < nPEM; ipnt++) {
            mPoints[ipnt]->prepareDeltaOutput();
        }
    }
    
    // coord
    if (enforceCoordTransform) {
        bool needRTZ = (chops.mWCS == channel::WavefieldCS::RTZ);
        if (chops.mNeedBufferU && displInRTZ() != needRTZ) {
            createCoordTransform();
        }
    }
}

// chi field
void FluidElement::getChiField(eigen::CMatXN &chi) const {
    // collect displacement from points
    collectDisplFromPoints(sDisplSpherical_FR);
    
    // convert to flattened
    mapPPvsN::PP2N(sDisplSpherical_FR, chi, mNu_1);
}

// displ field
void FluidElement::getDisplField(eigen::CMatXN3 &displ, bool needRTZ) const {
    // collect displacement from points
    collectDisplFromPoints(sDisplSpherical_FR);
    
    ///////////////////// compute stress /////////////////////
    // displ => strain
    mGradQuad->computeGrad3(sDisplSpherical_FR,
                            sStrainSpherical_FR, mNu_1);
    
    // strain => stress
    if (mPRT) {
        mTransform->transformSPZ_RTZ3(sStrainSpherical_FR, mNu_1);
        if (mInFourier) {
            mPRT->sphericalToUndulated3_FR(sStrainSpherical_FR,
                                           sStrainUndulated_FR, mNu_1);
            mAcoustic->strainToStress_FR(sStrainUndulated_FR,
                                         sStressUndulated_FR, mNu_1);
        } else {
            fft::gFFT_N3.computeC2R(sStrainSpherical_FR,
                                    sStrainSpherical_CD, mNr);
            mPRT->sphericalToUndulated3_CD(sStrainSpherical_CD,
                                           sStrainUndulated_CD, mNr);
            mAcoustic->strainToStress_CD(sStrainUndulated_CD,
                                         sStressUndulated_CD, mNr);
            fft::gFFT_N3.computeR2C(sStressUndulated_CD,
                                    sStressUndulated_FR, mNr);
        }
        if (!needRTZ) {
            mTransform->transformRTZ_SPZ3(sStressUndulated_FR, mNu_1);
        }
        // convert to flattened
        mapPPvsN::PP2N(sStressUndulated_FR, displ, mNu_1);
    } else {
        if (mInFourier) {
            mAcoustic->strainToStress_FR(sStrainSpherical_FR,
                                         sStressSpherical_FR, mNu_1);
        } else {
            fft::gFFT_N3.computeC2R(sStrainSpherical_FR,
                                    sStrainSpherical_CD, mNr);
            mAcoustic->strainToStress_CD(sStrainSpherical_CD,
                                         sStressSpherical_CD, mNr);
            fft::gFFT_N3.computeR2C(sStressSpherical_CD,
                                    sStressSpherical_FR, mNr);
        }
        if (needRTZ) {
            mTransform->transformSPZ_RTZ3(sStressSpherical_FR, mNu_1);
        }
        // convert to flattened
        mapPPvsN::PP2N(sStressSpherical_FR, displ, mNu_1);
    }
}

// pressure field
void FluidElement::getPressureField(eigen::CMatXN &pressure) const {
    for (int ipnt = 0; ipnt < nPEM; ipnt++) {
        mPoints[ipnt]->scatterPressureToElement(pressure, mNu_1, ipnt);
    }
}

// delta field
void FluidElement::getDeltaField(eigen::CMatXN &delta) const {
    for (int ipnt = 0; ipnt < nPEM; ipnt++) {
        mPoints[ipnt]->scatterDeltaToElement(delta, mNu_1, ipnt);
    }
}
