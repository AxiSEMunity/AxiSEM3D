//
//  FieldArithmetic.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 8/14/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  arithmetic operations on elemental fields
//  * operations such as s_{ij} = C_{ijkl} e{kl}
//  * external calling will be agnostic to 1D or 3D

#ifndef FieldArithmetic_hpp
#define FieldArithmetic_hpp

#include "Property.hpp"

#define blockN(dim, nx) block(0, dim * N, nx, N)
#define flattenDiag(P) Eigen::Map<const RRowN>(P.data()).asDiagonal()

// operation case
// _1D_FR: 1D object used in Fourier space
// _1D_CD: 1D object used in Cardinal space
// _3D_CD: 3D object used in Cardinal space
enum class CaseFA {_1D_FR, _1D_CD, _3D_CD};

template <int PED, int N = PED * PED>
class FieldArithmetic {
    typedef Eigen::Matrix<numerical::Real, 1, N> RRowN;
    typedef Property<PED> PropertyN;
    
public:
    ///////////////// R = F * P /////////////////
    // CASE == _1D_FR
    template <CaseFA CASE, class TypeR, class TypeF> static
    inline typename std::enable_if<CASE == CaseFA::_1D_FR, void>::type
    FxP(int nx, TypeR &R, int dimR, const TypeF &F,
        int dimF, const PropertyN &P) {
        R[nx][dimR] = F[nx][dimF].cwiseProduct(P.get1D());
    }
    // CASE == _1D_CD
    template <CaseFA CASE, class TypeR, class TypeF> static
    inline typename std::enable_if<CASE == CaseFA::_1D_CD, void>::type
    FxP(int nx, TypeR &R, int dimR, const TypeF &F,
        int dimF, const PropertyN &P) {
        R.blockN(dimR, nx) = F.blockN(dimF, nx) * flattenDiag(P.get1D());
    }
    // CASE == _3D_CD
    template <CaseFA CASE, class TypeR, class TypeF> static
    inline typename std::enable_if<CASE == CaseFA::_3D_CD, void>::type
    FxP(int nx, TypeR &R, int dimR, const TypeF &F,
        int dimF, const PropertyN &P) {
        R.blockN(dimR, nx) = F.blockN(dimF, nx).cwiseProduct(P.get3D());
    }
    
    
    ///////////////// R = F1 * P1 + F2 * P2 /////////////////
    // CASE == _1D_FR
    template <CaseFA CASE, class TypeR, class TypeF> static
    inline typename std::enable_if<CASE == CaseFA::_1D_FR, void>::type
    F1xP1_F2xP2(int nx, TypeR &R, int dimR, const TypeF &F,
                int dimF1, const PropertyN &P1,
                int dimF2, const PropertyN &P2) {
        R[nx][dimR] = (F[nx][dimF1].cwiseProduct(P1.get1D()) +
                       F[nx][dimF2].cwiseProduct(P2.get1D()));
    }
    // CASE == _1D_CD
    template <CaseFA CASE, class TypeR, class TypeF> static
    inline typename std::enable_if<CASE == CaseFA::_1D_CD, void>::type
    F1xP1_F2xP2(int nx, TypeR &R, int dimR, const TypeF &F,
                int dimF1, const PropertyN &P1,
                int dimF2, const PropertyN &P2) {
        R.blockN(dimR, nx) = (F.blockN(dimF1, nx) * flattenDiag(P1.get1D()) +
                              F.blockN(dimF2, nx) * flattenDiag(P2.get1D()));
    }
    // CASE == _3D_CD
    template <CaseFA CASE, class TypeR, class TypeF> static
    inline typename std::enable_if<CASE == CaseFA::_3D_CD, void>::type
    F1xP1_F2xP2(int nx, TypeR &R, int dimR, const TypeF &F,
                int dimF1, const PropertyN &P1,
                int dimF2, const PropertyN &P2) {
        R.blockN(dimR, nx) = (F.blockN(dimF1, nx).cwiseProduct(P1.get3D()) +
                              F.blockN(dimF2, nx).cwiseProduct(P2.get3D()));
    }
    
    
    ///////////////// R = F1 * P1 + F2 * P2 + F3 * P3 /////////////////
    // CASE == _1D_FR
    template <CaseFA CASE, class TypeR, class TypeF> static
    inline  typename std::enable_if<CASE == CaseFA::_1D_FR, void>::type
    F1xP1_F2xP2_F3xP3(int nx, TypeR &R, int dimR, const TypeF &F,
                      int dimF1, const PropertyN &P1,
                      int dimF2, const PropertyN &P2,
                      int dimF3, const PropertyN &P3) {
        R[nx][dimR] = (F[nx][dimF1].cwiseProduct(P1.get1D()) +
                       F[nx][dimF2].cwiseProduct(P2.get1D()) +
                       F[nx][dimF3].cwiseProduct(P3.get1D()));
    }
    // CASE == _1D_CD
    template <CaseFA CASE, class TypeR, class TypeF> static
    inline typename std::enable_if<CASE == CaseFA::_1D_CD, void>::type
    F1xP1_F2xP2_F3xP3(int nx, TypeR &R, int dimR, const TypeF &F,
                      int dimF1, const PropertyN &P1,
                      int dimF2, const PropertyN &P2,
                      int dimF3, const PropertyN &P3) {
        R.blockN(dimR, nx) = (F.blockN(dimF1, nx) * flattenDiag(P1.get1D()) +
                              F.blockN(dimF2, nx) * flattenDiag(P2.get1D()) +
                              F.blockN(dimF3, nx) * flattenDiag(P3.get1D()));
    }
    // CASE == _3D_CD
    template <CaseFA CASE, class TypeR, class TypeF> static
    inline typename std::enable_if<CASE == CaseFA::_3D_CD, void>::type
    F1xP1_F2xP2_F3xP3(int nx, TypeR &R, int dimR, const TypeF &F,
                      int dimF1, const PropertyN &P1,
                      int dimF2, const PropertyN &P2,
                      int dimF3, const PropertyN &P3) {
        R.blockN(dimR, nx) = (F.blockN(dimF1, nx).cwiseProduct(P1.get3D()) +
                              F.blockN(dimF2, nx).cwiseProduct(P2.get3D()) +
                              F.blockN(dimF3, nx).cwiseProduct(P3.get3D()));
    }
    
    
    ///////////////// R = R1 + F * P /////////////////
    // CASE == _1D_FR
    template <CaseFA CASE, class TypeR, class TypeF> static
    inline typename std::enable_if<CASE == CaseFA::_1D_FR, void>::type
    R1_FxP(int nx, TypeR &R, int dimR, int dimR1,
           const TypeF &F, int dimF, const PropertyN &P) {
        R[nx][dimR] = R[nx][dimR1] + F[nx][dimF].cwiseProduct(P.get1D());
    }
    // CASE == _1D_CD
    template <CaseFA CASE, class TypeR, class TypeF> static
    inline typename std::enable_if<CASE == CaseFA::_1D_CD, void>::type
    R1_FxP(int nx, TypeR &R, int dimR, int dimR1,
           const TypeF &F, int dimF, const PropertyN &P) {
        R.blockN(dimR, nx) = (R.blockN(dimR1, nx) +
                              F.blockN(dimF, nx) * flattenDiag(P.get1D()));
    }
    // CASE == _3D_CD
    template <CaseFA CASE, class TypeR, class TypeF> static
    inline typename std::enable_if<CASE == CaseFA::_3D_CD, void>::type
    R1_FxP(int nx, TypeR &R, int dimR, int dimR1,
           const TypeF &F, int dimF, const PropertyN &P) {
        R.blockN(dimR, nx) = (R.blockN(dimR1, nx) +
                              F.blockN(dimF, nx).cwiseProduct(P.get3D()));
    }
    
    
    ///////////////// R = (F1 + F2 + F3) * P /////////////////
    // CASE == _1D_FR
    template <CaseFA CASE, class TypeR, class TypeF> static
    inline typename std::enable_if<CASE == CaseFA::_1D_FR, void>::type
    F123xP(int nx, TypeR &R, int dimR, const TypeF &F,
           int dimF1, int dimF2, int dimF3, const PropertyN &P) {
        R[nx][dimR] = (F[nx][dimF1] +
                       F[nx][dimF2] +
                       F[nx][dimF3]).cwiseProduct(P.get1D());
    }
    // CASE == _1D_CD
    template <CaseFA CASE, class TypeR, class TypeF> static
    inline typename std::enable_if<CASE == CaseFA::_1D_CD, void>::type
    F123xP(int nx, TypeR &R, int dimR, const TypeF &F,
           int dimF1, int dimF2, int dimF3, const PropertyN &P) {
        R.blockN(dimR, nx) = (F.blockN(dimF1, nx) +
                              F.blockN(dimF2, nx) +
                              F.blockN(dimF3, nx)) * flattenDiag(P.get1D());
    }
    // CASE == _3D_CD
    template <CaseFA CASE, class TypeR, class TypeF> static
    inline typename std::enable_if<CASE == CaseFA::_3D_CD, void>::type
    F123xP(int nx, TypeR &R, int dimR, const TypeF &F,
           int dimF1, int dimF2, int dimF3, const PropertyN &P) {
        R.blockN(dimR, nx) = (F.blockN(dimF1, nx) +
                              F.blockN(dimF2, nx) +
                              F.blockN(dimF3, nx)).cwiseProduct(P.get3D());
    }
    
    
    ///////////////// R = F /////////////////
    // CASE == _1D_FR
    template <CaseFA CASE, class TypeR, class TypeF> static
    inline typename std::enable_if<CASE == CaseFA::_1D_FR, void>::type
    F(int nx, TypeR &R, int dimR, const TypeF &F, int dimF) {
        R[nx][dimR] = F[nx][dimF];
    }
    // CASE != _1D_FR
    template <CaseFA CASE, class TypeR, class TypeF> static
    inline typename std::enable_if<CASE != CaseFA::_1D_FR, void>::type
    F(int nx, TypeR &R, int dimR, const TypeF &F, int dimF) {
        R.blockN(dimR, nx) = F.blockN(dimF, nx);
    }
    
    
    ///////////// R = (F1a + F1b) * P1 + F2 * P2 + F3 * P3 /////////////
    // CASE == _1D_FR
    template <CaseFA CASE, class TypeR, class TypeF> static
    inline typename std::enable_if<CASE == CaseFA::_1D_FR, void>::type
    F1abxP1_F2xP2_F3xP3(int nx, TypeR &R, int dimR, const TypeF &F,
                        int dimF1a, int dimF1b, const PropertyN &P1,
                        int dimF2, const PropertyN &P2,
                        int dimF3, const PropertyN &P3) {
        R[nx][dimR] = ((F[nx][dimF1a] + F[nx][dimF1b])
                       .cwiseProduct(P1.get1D()) +
                       F[nx][dimF2].cwiseProduct(P2.get1D()) +
                       F[nx][dimF3].cwiseProduct(P3.get1D()));
    }
    // CASE == _1D_CD
    template <CaseFA CASE, class TypeR, class TypeF> static
    inline typename std::enable_if<CASE == CaseFA::_1D_CD, void>::type
    F1abxP1_F2xP2_F3xP3(int nx, TypeR &R, int dimR, const TypeF &F,
                        int dimF1a, int dimF1b, const PropertyN &P1,
                        int dimF2, const PropertyN &P2,
                        int dimF3, const PropertyN &P3) {
        R.blockN(dimR, nx) = ((F.blockN(dimF1a, nx) + F.blockN(dimF1b, nx))
                              * flattenDiag(P1.get1D())
                              + F.blockN(dimF2, nx) * flattenDiag(P2.get1D())
                              + F.blockN(dimF3, nx) * flattenDiag(P3.get1D()));
    }
    // CASE == _3D_CD
    template <CaseFA CASE, class TypeR, class TypeF> static
    inline typename std::enable_if<CASE == CaseFA::_3D_CD, void>::type
    F1abxP1_F2xP2_F3xP3(int nx, TypeR &R, int dimR, const TypeF &F,
                        int dimF1a, int dimF1b, const PropertyN &P1,
                        int dimF2, const PropertyN &P2,
                        int dimF3, const PropertyN &P3) {
        R.blockN(dimR, nx) = ((F.blockN(dimF1a, nx) + F.blockN(dimF1b, nx))
                              .cwiseProduct(P1.get3D()) +
                              F.blockN(dimF2, nx).cwiseProduct(P2.get3D()) +
                              F.blockN(dimF3, nx).cwiseProduct(P3.get3D()));
    }
    
    
    ///////////////// R = F1 + F2 /////////////////
    // CASE == _1D_FR
    template <CaseFA CASE, class TypeR, class TypeF> static
    inline typename std::enable_if<CASE == CaseFA::_1D_FR, void>::type
    F1_F2(int nx, TypeR &R, int dimR, const TypeF &F,
          int dimF1, int dimF2) {
        R[nx][dimR] = F[nx][dimF1] + F[nx][dimF2];
    }
    // CASE != _1D_FR
    template <CaseFA CASE, class TypeR, class TypeF> static
    inline typename std::enable_if<CASE != CaseFA::_1D_FR, void>::type
    F1_F2(int nx, TypeR &R, int dimR, const TypeF &F,
          int dimF1, int dimF2) {
        R.blockN(dimR, nx) = F.blockN(dimF1, nx) + F.blockN(dimF2, nx);
    }
    
    
    ///////////////// R = R1 * P1 + F2 * P2 /////////////////
    // CASE == _1D_FR
    template <CaseFA CASE, class TypeR, class TypeF> static
    inline typename std::enable_if<CASE == CaseFA::_1D_FR, void>::type
    R1xP1_F2xP2(int nx, TypeR &R, int dimR, int dimR1, const PropertyN &P1,
                const TypeF &F, int dimF2, const PropertyN &P2) {
        R[nx][dimR] = (R[nx][dimR1].cwiseProduct(P1.get1D()) +
                       F[nx][dimF2].cwiseProduct(P2.get1D()));
    }
    // CASE == _1D_CD
    template <CaseFA CASE, class TypeR, class TypeF> static
    inline typename std::enable_if<CASE == CaseFA::_1D_CD, void>::type
    R1xP1_F2xP2(int nx, TypeR &R, int dimR, int dimR1, const PropertyN &P1,
                const TypeF &F, int dimF2, const PropertyN &P2) {
        R.blockN(dimR, nx) = (R.blockN(dimR1, nx) * flattenDiag(P1.get1D()) +
                              F.blockN(dimF2, nx) * flattenDiag(P2.get1D()));
    }
    // CASE == _3D_CD
    template <CaseFA CASE, class TypeR, class TypeF> static
    inline typename std::enable_if<CASE == CaseFA::_3D_CD, void>::type
    R1xP1_F2xP2(int nx, TypeR &R, int dimR, int dimR1, const PropertyN &P1,
                const TypeF &F, int dimF2, const PropertyN &P2) {
        R.blockN(dimR, nx) = (R.blockN(dimR1, nx).cwiseProduct(P1.get3D()) +
                              F.blockN(dimF2, nx).cwiseProduct(P2.get3D()));
    }
    
    
    ///////////////// R = F1 * P1 + F2 * P2 + ... + F6 * P6 /////////////////
    // CASE == _1D_FR
    template <CaseFA CASE, class TypeR, class TypeF> static
    inline typename std::enable_if<CASE == CaseFA::_1D_FR, void>::type
    F1xP1_to_F6xP6(int nx, TypeR &R, int dimR, const TypeF &F,
                   int dimF1, const PropertyN &P1,
                   int dimF2, const PropertyN &P2,
                   int dimF3, const PropertyN &P3,
                   int dimF4, const PropertyN &P4,
                   int dimF5, const PropertyN &P5,
                   int dimF6, const PropertyN &P6) {
        R[nx][dimR] = (F[nx][dimF1].cwiseProduct(P1.get1D()) +
                       F[nx][dimF2].cwiseProduct(P2.get1D()) +
                       F[nx][dimF3].cwiseProduct(P3.get1D()) +
                       F[nx][dimF4].cwiseProduct(P4.get1D()) +
                       F[nx][dimF5].cwiseProduct(P5.get1D()) +
                       F[nx][dimF6].cwiseProduct(P6.get1D()));
    }
    // CASE == _1D_CD
    template <CaseFA CASE, class TypeR, class TypeF> static
    inline typename std::enable_if<CASE == CaseFA::_1D_CD, void>::type
    F1xP1_to_F6xP6(int nx, TypeR &R, int dimR, const TypeF &F,
                   int dimF1, const PropertyN &P1,
                   int dimF2, const PropertyN &P2,
                   int dimF3, const PropertyN &P3,
                   int dimF4, const PropertyN &P4,
                   int dimF5, const PropertyN &P5,
                   int dimF6, const PropertyN &P6) {
        R.blockN(dimR, nx) = (F.blockN(dimF1, nx) * flattenDiag(P1.get1D()) +
                              F.blockN(dimF2, nx) * flattenDiag(P2.get1D()) +
                              F.blockN(dimF3, nx) * flattenDiag(P3.get1D()) +
                              F.blockN(dimF4, nx) * flattenDiag(P4.get1D()) +
                              F.blockN(dimF5, nx) * flattenDiag(P5.get1D()) +
                              F.blockN(dimF6, nx) * flattenDiag(P6.get1D()));
    }
    // CASE == _3D_CD
    template <CaseFA CASE, class TypeR, class TypeF> static
    inline typename std::enable_if<CASE == CaseFA::_3D_CD, void>::type
    F1xP1_to_F6xP6(int nx, TypeR &R, int dimR, const TypeF &F,
                   int dimF1, const PropertyN &P1,
                   int dimF2, const PropertyN &P2,
                   int dimF3, const PropertyN &P3,
                   int dimF4, const PropertyN &P4,
                   int dimF5, const PropertyN &P5,
                   int dimF6, const PropertyN &P6) {
        R.blockN(dimR, nx) = (F.blockN(dimF1, nx).cwiseProduct(P1.get3D()) +
                              F.blockN(dimF2, nx).cwiseProduct(P2.get3D()) +
                              F.blockN(dimF3, nx).cwiseProduct(P3.get3D()) +
                              F.blockN(dimF4, nx).cwiseProduct(P4.get3D()) +
                              F.blockN(dimF5, nx).cwiseProduct(P5.get3D()) +
                              F.blockN(dimF6, nx).cwiseProduct(P6.get3D()));
    }
};

// specialization for all GLL points
namespace faN {
    // property on all GLL points
    typedef Property<spectral::nPED> PropertyN;
    
    // operations on all GLL points
    typedef FieldArithmetic<spectral::nPED> FieldArithmeticN;
}

#endif /* FieldArithmetic_hpp */
