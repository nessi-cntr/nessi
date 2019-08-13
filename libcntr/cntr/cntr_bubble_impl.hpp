#ifndef CNTR_BUBBLE_IMPL_H
#define CNTR_BUBBLE_IMPL_H

#include "cntr_bubble_decl.hpp"
#include "cntr_herm_matrix_timestep_view_decl.hpp"

namespace cntr {

///////////////////////////////////////////////////////////////////////////////////

//  BUBBLE 1 :  C(t,t') = ii * A(t,t') * B(t',t)
//  BUBBLE 2 :  C(t,t') = ii * A(t,t') * B(t,t')

///////////////////////////////////////////////////////////////////////////////////

// -------------  Auxiliary routines: -------------

/// @private
// Matsubara:
// C_{c1,c2}(tau) = - A_{a1,a2}(tau) * B_{b2,b1}(-tau)
//                = - A_{a1,a2}(tau) * B_{b2,b1}(beta-tau) (no cc needed !!)
template <typename T>
void get_bubble_1_mat(std::complex<T> *cmat, int sc, int c1, int c2, std::complex<T> *amat,
                      int sa, int a1, int a2, std::complex<T> *bmat, int sb, int b1, int b2,
                      int sigb, int ntau) {
    int m;
    int a12 = a1 * sa + a2, sa2 = sa * sa;
    int b21 = b2 * sb + b1, sb2 = sb * sb;
    int c12 = c1 * sc + c2, sc2 = sc * sc;
    double sig = -1.0 * sigb;
    for (m = 0; m <= ntau; m++)
        cmat[c12 + m * sc2] = sig * amat[m * sa2 + a12] * bmat[(ntau - m) * sb2 + b21];
}
/// @private
template <typename T>
void
get_bubble_1_timestep(int tstp, std::complex<T> *cret, std::complex<T> *ctv,
                      std::complex<T> *cles, int sc, int c1, int c2, std::complex<T> *aret,
                      std::complex<T> *atv, std::complex<T> *ales, std::complex<T> *accret,
                      std::complex<T> *acctv, std::complex<T> *accles, int sa, int a1,
                      int a2, std::complex<T> *bret, std::complex<T> *btv,
                      std::complex<T> *bles, std::complex<T> *bccret, std::complex<T> *bcctv,
                      std::complex<T> *bccles, int sb, int b1, int b2, int sigb, int ntau) {
    int m;
    int a12 = a1 * sa + a2, a21 = a2 * sa + a1, sa2 = sa * sa;
    int b12 = b1 * sb + b2, b21 = b2 * sb + b1, sb2 = sb * sb;
    int c12 = c1 * sc + c2, sc2 = sc * sc;
    std::complex<T> msigb = -1.0 * sigb, ii = std::complex<T>(0, 1.0);
    std::complex<T> bgtr21_tt1, cgtr_tt1, cles_tt1, ales12_tt1, bgtr21_t1t, agtr12_tt1,
        bles21_t1t, bvt21, bles21_tt1;
    for (m = 0; m <= tstp; m++) {
        // Bles_{21}(tstp,m) = - Bccles_{12}(m,tstp)^*
        bles21_tt1 = -conj(bccles[m * sb2 + b12]);
        // Bgtr_{21}(tstp,m) = Bret_{21}(tstp,m) - bles21_tt1;
        bgtr21_tt1 = bret[m * sb2 + b21] + bles21_tt1;
        // bgtr21_t1t = - Bccgtr_{12}(t,t1)^*
        //            = - [ Bccret_{12}(t,t1) + Bccles_{12}(t,t1)  ]^*
        //            = - Bccret_{12}(t,t1)^* + Bles_{21}(t1,t)
        bles21_t1t = bles[m * sb2 + b21];
        bgtr21_t1t = -conj(bccret[m * sb2 + b12]) + bles21_t1t;
        // Ales_{12}(tstp,m) = -Accles_{21}(m,tstp)^*
        ales12_tt1 = -conj(accles[sa2 * m + a21]);
        // Agtr_{a1,a2}(tstp,m) = Aret_{a1,a2}(tstp,m) - Ales_{a1,a2}(tstp,m)
        agtr12_tt1 = aret[m * sa2 + a12] + ales12_tt1;
        // Cgtr_{12}(tstp,m) = ii * Agtr_{12}(tstp,m)*Bles_{21}(m,tstp)
        cgtr_tt1 = ii * agtr12_tt1 * bles21_t1t;
        // Cles_{12}(tstp,m) = ii * Ales_{12}(tstp,m)*Bgtr_{21}(m,tstp)
        cles_tt1 = ii * ales12_tt1 * bgtr21_t1t;
        // Cret_{12}(tstp,m) = Cgtr_{12}(tstp,m) - Cles_{12}(tstp,m)
        cret[m * sc2 + c12] = cgtr_tt1 - cles_tt1;
        // Cles_{12}(m,tstp) = ii * Ales_{12}(m,tstp)*Bgtr_{21}(tstp,m)
        cles[m * sc2 + c12] = ii * ales[m * sa2 + a12] * bgtr21_tt1;
    }
    for (m = 0; m <= ntau; m++) {
        bvt21 = msigb * conj(bcctv[(ntau - m) * sb2 + b12]);
        ctv[m * sc2 + c12] = ii * atv[m * sa2 + a12] * bvt21;
    }
}

//  BUBBLE 2 :
//  C(t,t') = ii * A(t,t') * B(t,t')

/// @private
// Matsubara:
// C_{c1,c2}(tau) = - A_{a1,a2}(tau) * B_{b1,b2}(tau)
template <typename T>
void get_bubble_2_mat(std::complex<T> *cmat, int sc, int c1, int c2, std::complex<T> *amat,
                      int sa, int a1, int a2, std::complex<T> *bmat, int sb, int b1, int b2,
                      int ntau) {
    int m;
    int a12 = a1 * sa + a2, sa2 = sa * sa;
    int b12 = b1 * sb + b2, sb2 = sb * sb;
    int c12 = c1 * sc + c2, sc2 = sc * sc;
    for (m = 0; m <= ntau; m++)
        cmat[c12 + m * sc2] = -amat[m * sa2 + a12] * bmat[m * sb2 + b12];
}
/// @private
template <typename T>
void
get_bubble_2_timestep(int tstp, std::complex<T> *cret, std::complex<T> *ctv,
                      std::complex<T> *cles, int sc, int c1, int c2, std::complex<T> *aret,
                      std::complex<T> *atv, std::complex<T> *ales, std::complex<T> *accret,
                      std::complex<T> *acctv, std::complex<T> *accles, int sa, int a1,
                      int a2, std::complex<T> *bret, std::complex<T> *btv,
                      std::complex<T> *bles, std::complex<T> *bccret, std::complex<T> *bcctv,
                      std::complex<T> *bccles, int sb, int b1, int b2, int ntau) {
    int m;
    int a12 = a1 * sa + a2, a21 = a2 * sa + a1, sa2 = sa * sa;
    int b12 = b1 * sb + b2, b21 = b2 * sb + b1, sb2 = sb * sb;
    int c12 = c1 * sc + c2, sc2 = sc * sc;
    std::complex<T> ii = std::complex<T>(0, 1.0);
    std::complex<T> bgtr12_tt1, agtr12_tt1, cgtr_tt1, cles_tt1;
    for (m = 0; m <= tstp; m++) {
        bgtr12_tt1 = bret[m * sb2 + b12] - conj(bccles[m * sb2 + b21]);
        agtr12_tt1 = aret[m * sa2 + a12] - conj(accles[m * sa2 + a21]);
        // Cgtr_{12}(tstp,m) = ii * Agtr_{12}(tstp,m)*Bles_{21}(m,tstp)
        cgtr_tt1 = ii * agtr12_tt1 * bgtr12_tt1;
        // Cles_{12}(tstp,m) = ii * Ales_{12}(tstp,m)*Bles_{12}(tstp,m)
        cles_tt1 = ii * conj(accles[m * sa2 + a21]) * conj(bccles[m * sb2 + b21]);
        // Cret_{12}(tstp,m) = Cgtr_{12}(tstp,m) - Cles_{12}(tstp,m)
        cret[m * sc2 + c12] = cgtr_tt1 - cles_tt1;
        // Cles_{12}(m,tstp) = ii * Ales_{12}(m,tstp)*Bgtr_{12}(m,tstp)
        cles[m * sc2 + c12] = ii * ales[m * sa2 + a12] * bles[m * sb2 + b12];
    }
    for (m = 0; m <= ntau; m++) {
        ctv[m * sc2 + c12] = ii * atv[m * sa2 + a12] * btv[m * sb2 + b12];
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
////// MAIN FUNCTIONS

/// @private
// A,B hermitian, orbital dimension > 1:
// C_{c1,c2}(t1,t2) = ii * A_{a1,a2}(t1,t2) * B_{b2,b1}(t2,t1)

template <typename T>
void Bubble1(int tstp, herm_matrix_timestep_view<T> &C, int c1, int c2,
             herm_matrix_timestep_view<T> &A, herm_matrix_timestep_view<T> &Acc, int a1,
             int a2, herm_matrix_timestep_view<T> &B, herm_matrix_timestep_view<T> &Bcc,
             int b1, int b2) {
    int ntau = C.ntau();
    assert(ntau == A.ntau());
    assert(ntau == B.ntau());
    assert(ntau == Acc.ntau());
    assert(ntau == Bcc.ntau());
    assert(tstp == C.tstp_);
    assert(tstp == B.tstp_);
    assert(tstp == A.tstp_);
    assert(tstp == Acc.tstp_);
    assert(tstp == Bcc.tstp_);
    assert(a1 <= A.size1());
    assert(a2 <= A.size1());
    assert(b1 <= B.size1());
    assert(b2 <= B.size1());
    assert(c1 <= C.size1());
    assert(c2 <= C.size1());
    assert(Acc.size1() ==  A.size1());
    assert(Bcc.size1() == B.size1());

    if (tstp == -1) {
        get_bubble_1_mat(C.matptr(0), C.size1(), c1, c2, A.matptr(0), A.size1(), a1, a2,
                         B.matptr(0), B.size1(), b1, b2, B.sig(), ntau);
    } else {
        get_bubble_1_timestep(tstp, C.retptr(0), C.tvptr(0), C.lesptr(0), C.size1(), c1, c2,
                              A.retptr(0), A.tvptr(0), A.lesptr(0), Acc.retptr(0),
                              Acc.tvptr(0), Acc.lesptr(0), A.size1(), a1, a2, B.retptr(0),
                              B.tvptr(0), B.lesptr(0), Bcc.retptr(0), Bcc.tvptr(0),
                              Bcc.lesptr(0), B.size1(), b1, b2, B.sig(), ntau);
    }
}
/// @private
template <typename T>
void Bubble1(int tstp, herm_matrix_timestep_view<T> &C, int c1, int c2,
             herm_matrix_timestep_view<T> &A, int a1, int a2,
             herm_matrix_timestep_view<T> &B, int b1, int b2) {
    return Bubble1(tstp, C, c1, c2, A, A, a1, a2, B, B, b1, b2);
}
/// @private
template <typename T>
void Bubble1(int tstp, herm_matrix_timestep_view<T> &C, herm_matrix_timestep_view<T> &A,
             herm_matrix_timestep_view<T> &Acc, herm_matrix_timestep_view<T> &B,
             herm_matrix_timestep_view<T> &Bcc) {
    return Bubble1(tstp, C, 0, 0, A, Acc, 0, 0, B, Bcc, 0, 0);
}
/// @private
template <typename T>
void Bubble1(int tstp, herm_matrix_timestep_view<T> &C, herm_matrix_timestep_view<T> &A,
             herm_matrix_timestep_view<T> &B) {
    return Bubble1(tstp, C, A, A, B, B);
}
/// @private
template <class GGC, class GGA, class GGB>
void Bubble1(int tstp, GGC &C, int c1, int c2, GGA &A, GGA &Acc, int a1, int a2, GGB &B,
             GGB &Bcc, int b1, int b2) {
    herm_matrix_timestep_view<typename GGC::scalar_type> ctmp(tstp, C);
    herm_matrix_timestep_view<typename GGA::scalar_type> atmp(tstp, A);
    herm_matrix_timestep_view<typename GGA::scalar_type> acctmp(tstp, Acc);
    herm_matrix_timestep_view<typename GGB::scalar_type> btmp(tstp, B);
    herm_matrix_timestep_view<typename GGB::scalar_type> bcctmp(tstp, Bcc);
    Bubble1(tstp, ctmp, c1, c2, atmp, acctmp, a1, a2, btmp, bcctmp, b1, b2);
}
/// @private
template <class GGC, class GGA, class GGB>
void Bubble1(int tstp, GGC &C, int c1, int c2, GGA &A, int a1, int a2, GGB &B, int b1,
             int b2) {
    herm_matrix_timestep_view<typename GGC::scalar_type> ctmp(tstp, C);
    herm_matrix_timestep_view<typename GGA::scalar_type> atmp(tstp, A);
    herm_matrix_timestep_view<typename GGB::scalar_type> btmp(tstp, B);
    Bubble1(tstp, ctmp, c1, c2, atmp, a1, a2, btmp, b1, b2);
}
/// @private
template <class GGC, class GGA, class GGB>
void Bubble1(int tstp, GGC &C, GGA &A, GGA &Acc, GGB &B, GGB &Bcc) {
    herm_matrix_timestep_view<typename GGC::scalar_type> ctmp(tstp, C);
    herm_matrix_timestep_view<typename GGA::scalar_type> atmp(tstp, A);
    herm_matrix_timestep_view<typename GGA::scalar_type> acctmp(tstp, Acc);
    herm_matrix_timestep_view<typename GGB::scalar_type> btmp(tstp, B);
    herm_matrix_timestep_view<typename GGB::scalar_type> bcctmp(tstp, Bcc);
    Bubble1(tstp, ctmp, atmp, acctmp, btmp, bcctmp);
}
/// @private
template <class GGC, class GGA, class GGB> void Bubble1(int tstp, GGC &C, GGA &A, GGB &B) {
    herm_matrix_timestep_view<typename GGC::scalar_type> ctmp(tstp, C);
    herm_matrix_timestep_view<typename GGA::scalar_type> atmp(tstp, A);
    herm_matrix_timestep_view<typename GGB::scalar_type> btmp(tstp, B);
    Bubble1(tstp, ctmp, atmp, btmp);
}

/// @private
// A,B hermitian, orbital dimension > 1:
// C_{c1,c2}(t1,t2) = ii * A_{a1,a2}(t1,t2) * B_{b1,b2}(t1,t2)
template <typename T>
void Bubble2(int tstp, herm_matrix_timestep_view<T> &C, int c1, int c2,
             herm_matrix_timestep_view<T> &A, herm_matrix_timestep_view<T> &Acc, int a1,
             int a2, herm_matrix_timestep_view<T> &B, herm_matrix_timestep_view<T> &Bcc,
             int b1, int b2) {
    int ntau = C.ntau();

    assert(ntau == A.ntau());
    assert(ntau == B.ntau());
    assert(ntau == Acc.ntau());
    assert(ntau == Bcc.ntau());
    assert(tstp == C.tstp_);
    assert(tstp == B.tstp_);
    assert(tstp == A.tstp_);
    assert(tstp == Acc.tstp_);
    assert(tstp == Bcc.tstp_);
    assert(a1 <= A.size1());
    assert(a2 <= A.size1());
    assert(b1 <= B.size1());
    assert(b2 <= B.size1());
    assert(c1 <= C.size1());
    assert(c2 <= C.size1());
    assert(Acc.size1() ==  A.size1());
    assert(Bcc.size1() == B.size1());

    if (tstp == -1) {
        get_bubble_2_mat(C.matptr(0), C.size1(), c1, c2, A.matptr(0), A.size1(), a1, a2,
                         B.matptr(0), B.size1(), b1, b2, ntau);
    } else {
        get_bubble_2_timestep(tstp, C.retptr(0), C.tvptr(0), C.lesptr(0), C.size1(), c1, c2,
                              A.retptr(0), A.tvptr(0), A.lesptr(0), Acc.retptr(0),
                              Acc.tvptr(0), Acc.lesptr(0), A.size1(), a1, a2, B.retptr(0),
                              B.tvptr(0), B.lesptr(0), Bcc.retptr(0), Bcc.tvptr(0),
                              Bcc.lesptr(0), B.size1(), b1, b2, ntau);
    }
}
/// @private
template <typename T>
void Bubble2(int tstp, herm_matrix_timestep_view<T> &C, int c1, int c2,
             herm_matrix_timestep_view<T> &A, int a1, int a2,
             herm_matrix_timestep_view<T> &B, int b1, int b2) {
    return Bubble2(tstp, C, c1, c2, A, A, a1, a2, B, B, b1, b2);
}
/// @private
template <typename T>
void Bubble2(int tstp, herm_matrix_timestep_view<T> &C, herm_matrix_timestep_view<T> &A,
             herm_matrix_timestep_view<T> &Acc, herm_matrix_timestep_view<T> &B,
             herm_matrix_timestep_view<T> &Bcc) {
    return Bubble2(tstp, C, 0, 0, A, Acc, 0, 0, B, Bcc, 0, 0);
}
/// @private
template <typename T>
void Bubble2(int tstp, herm_matrix_timestep_view<T> &C, herm_matrix_timestep_view<T> &A,
             herm_matrix_timestep_view<T> &B) {
    return Bubble2(tstp, C, A, A, B, B);
}
/// @private
template <class GGC, class GGA, class GGB>
void Bubble2(int tstp, GGC &C, int c1, int c2, GGA &A, GGA &Acc, int a1, int a2, GGB &B,
             GGB &Bcc, int b1, int b2) {
    herm_matrix_timestep_view<typename GGC::scalar_type> ctmp(tstp, C);
    herm_matrix_timestep_view<typename GGA::scalar_type> atmp(tstp, A);
    herm_matrix_timestep_view<typename GGA::scalar_type> acctmp(tstp, Acc);
    herm_matrix_timestep_view<typename GGB::scalar_type> btmp(tstp, B);
    herm_matrix_timestep_view<typename GGB::scalar_type> bcctmp(tstp, Bcc);
    Bubble2(tstp, ctmp, c1, c2, atmp, acctmp, a1, a2, btmp, bcctmp, b1, b2);
}
/// @private
template <class GGC, class GGA, class GGB>
void Bubble2(int tstp, GGC &C, int c1, int c2, GGA &A, int a1, int a2, GGB &B, int b1,
             int b2) {
    herm_matrix_timestep_view<typename GGC::scalar_type> ctmp(tstp, C);
    herm_matrix_timestep_view<typename GGA::scalar_type> atmp(tstp, A);
    herm_matrix_timestep_view<typename GGB::scalar_type> btmp(tstp, B);
    Bubble2(tstp, ctmp, c1, c2, atmp, a1, a2, btmp, b1, b2);
}
/// @private
template <class GGC, class GGA, class GGB>
void Bubble2(int tstp, GGC &C, GGA &A, GGA &Acc, GGB &B, GGB &Bcc) {
    herm_matrix_timestep_view<typename GGC::scalar_type> ctmp(tstp, C);
    herm_matrix_timestep_view<typename GGA::scalar_type> atmp(tstp, A);
    herm_matrix_timestep_view<typename GGA::scalar_type> acctmp(tstp, Acc);
    herm_matrix_timestep_view<typename GGB::scalar_type> btmp(tstp, B);
    herm_matrix_timestep_view<typename GGB::scalar_type> bcctmp(tstp, Bcc);
    Bubble2(tstp, ctmp, atmp, acctmp, btmp, bcctmp);
}
/// @private
template <class GGC, class GGA, class GGB> void Bubble2(int tstp, GGC &C, GGA &A, GGB &B) {
    herm_matrix_timestep_view<typename GGC::scalar_type> ctmp(tstp, C);
    herm_matrix_timestep_view<typename GGA::scalar_type> atmp(tstp, A);
    herm_matrix_timestep_view<typename GGB::scalar_type> btmp(tstp, B);
    Bubble2(tstp, ctmp, atmp, btmp);
}

} // namespace cntr

#endif  // CNTR_BUBBLE_IMPL_H
