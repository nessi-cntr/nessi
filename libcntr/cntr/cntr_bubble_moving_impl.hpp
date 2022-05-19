#ifndef CNTR_BUBBLE_MOVING_IMPL_H
#define CNTR_BUBBLE_MOVING_IMPL_H

#include "cntr_bubble_moving_decl.hpp"
#include "cntr_herm_matrix_timestep_moving_view_decl.hpp"

namespace cntr {
  // Bubbles for moving propagators
  // C(t-tstp,t-tstp-m)= ii * A(t-tstp,t-tstp-m) * B(t-tstp-m,t-tstp), l=0...tc
  // where aret[m] -> A(t-tstp,t-tstp-m) etc.
/// @private
  template<typename T>
  void get_bubble_1_timestep_moving(int tc,std::complex<T> *cret,std::complex<T> *cles,int sc,int c1,int c2,std::complex<T> *aret,std::complex<T> *ales,std::complex<T> *accret,std::complex<T> *accles,int sa,int a1,int a2,std::complex<T> *bret,std::complex<T> *bles,std::complex<T> *bccret,std::complex<T> *bccles,int sb,int b1,int b2,int sigb){
    int m;
    int a12=a1*sa+a2,a21=a2*sa+a1,sa2=sa*sa;
    int b12=b1*sb+b2,b21=b2*sb+b1,sb2=sb*sb;
    int c12=c1*sc+c2,sc2=sc*sc;
    std::complex<T> msigb=-1.0*sigb,ii=std::complex<T>(0,1.0);
    std::complex<T> bgtr21_tt1,cgtr_tt1,cles_tt1,ales12_tt1,bgtr21_t1t,agtr12_tt1,bles21_t1t;
    for(m=0;m<=tc;m++){
      //Ales{12}(t,t-m)
      ales12_tt1 = ales[m*sa2+a12];
      //Agtr{12}(t,t-m) = Ales{12}(t,t-m) + Aret{12}(t,t-m)
      agtr12_tt1 = ales12_tt1 + aret[m*sa2+a12];
      //Bles{21}(t-m,t)
      bles21_t1t = -conj(bccles[m*sb2+b12]);
      //Bgtr{21}(t-m,t) =
      bgtr21_t1t = bles21_t1t-conj(bccret[m*sb2+b12]);
      cgtr_tt1 = ii*agtr12_tt1*bles21_t1t;
      cles_tt1 = ii*ales12_tt1*bgtr21_t1t;
      cret[m*sc2+c12] = cgtr_tt1-cles_tt1;
      cles[m*sc2+c12] = ii*ales12_tt1*bgtr21_t1t;
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////// MAIN FUNCTIONS

  // A,B hermitian, orbital dimension > 1:
  // C_{c1,c2}(t1,t2) = ii * A_{a1,a2}(t1,t2) * B_{b2,b1}(t2,t1)
/// @private
  template <typename T>
  void Bubble1_moving(int tstp, herm_matrix_timestep_moving_view<T> &C, int c1, int c2,
		      herm_matrix_timestep_moving_view<T> &A, herm_matrix_timestep_moving_view<T> &Acc, int a1,
		      int a2, herm_matrix_timestep_moving_view<T> &B, herm_matrix_timestep_moving_view<T> &Bcc,
		      int b1, int b2) {
    int tc = C.tc();
    // assert(ntau == A.ntau());
    // assert(ntau == B.ntau());
    // assert(ntau == Acc.ntau());
    // assert(ntau == Bcc.ntau());
    // assert(tstp == C.tstp_);
    // assert(tstp == B.tstp_);
    // assert(tstp == A.tstp_);
    // assert(tstp == Acc.tstp_);
    // assert(tstp == Bcc.tstp_);
    assert(a1 <= A.size1());
    assert(a2 <= A.size1());
    assert(b1 <= B.size1());
    assert(b2 <= B.size1());
    assert(c1 <= C.size1());
    assert(c2 <= C.size1());
    assert(Acc.size1() ==  A.size1());
    assert(Bcc.size1() == B.size1());

    get_bubble_1_timestep_moving(tc,C.retptr(0),C.lesptr(0),C.size1(),c1,c2,
  				 A.retptr(0),A.lesptr(0),Acc.retptr(0),Acc.lesptr(0),A.size1(),a1,a2,
  				 B.retptr(0),B.lesptr(0),Bcc.retptr(0),Bcc.lesptr(0),B.size1(),b1,b2,
  				 B.sig());
  
  }
  /// @private
  template <typename T>
  void Bubble1_moving(int tstp, herm_matrix_timestep_moving_view<T> &C, int c1, int c2,
		      herm_matrix_timestep_moving_view<T> &A, int a1, int a2,
		      herm_matrix_timestep_moving_view<T> &B, int b1, int b2) {
    return Bubble1_moving(tstp, C, c1, c2, A, A, a1, a2, B, B, b1, b2);
  }
  /// @private
  template <typename T>
  void Bubble1_moving(int tstp, herm_matrix_timestep_moving_view<T> &C, herm_matrix_timestep_moving_view<T> &A,
		      herm_matrix_timestep_moving_view<T> &Acc, herm_matrix_timestep_moving_view<T> &B,
		      herm_matrix_timestep_moving_view<T> &Bcc) {
    return Bubble1_moving(tstp, C, 0, 0, A, Acc, 0, 0, B, Bcc, 0, 0);
  }
  /// @private
  template <typename T>
  void Bubble1_moving(int tstp, herm_matrix_timestep_moving_view<T> &C, herm_matrix_timestep_moving_view<T> &A,
		      herm_matrix_timestep_moving_view<T> &B) {
    return Bubble1_moving(tstp, C, A, A, B, B);
  }

  template <class GGC, class GGA, class GGB>
  void Bubble1_moving(int tstp, GGC &C, int c1, int c2, GGA &A, GGA &Acc, int a1, int a2, GGB &B,
		      GGB &Bcc, int b1, int b2) {
    herm_matrix_timestep_moving_view<typename GGC::scalar_type> ctmp(tstp, C);
    herm_matrix_timestep_moving_view<typename GGA::scalar_type> atmp(tstp, A);
    herm_matrix_timestep_moving_view<typename GGA::scalar_type> acctmp(tstp, Acc);
    herm_matrix_timestep_moving_view<typename GGB::scalar_type> btmp(tstp, B);
    herm_matrix_timestep_moving_view<typename GGB::scalar_type> bcctmp(tstp, Bcc);
    Bubble1_moving(tstp, ctmp, c1, c2, atmp, acctmp, a1, a2, btmp, bcctmp, b1, b2);
  }
  template <class GGC, class GGA, class GGB>
  void Bubble1_moving(int tstp, GGC &C, int c1, int c2, GGA &A, int a1, int a2, GGB &B, int b1,
		      int b2) {
    herm_matrix_timestep_moving_view<typename GGC::scalar_type> ctmp(tstp, C);
    herm_matrix_timestep_moving_view<typename GGA::scalar_type> atmp(tstp, A);
    herm_matrix_timestep_moving_view<typename GGB::scalar_type> btmp(tstp, B);
    Bubble1_moving(tstp, ctmp, c1, c2, atmp, a1, a2, btmp, b1, b2);
  }
  template <class GGC, class GGA, class GGB>
  void Bubble1_moving(int tstp, GGC &C, GGA &A, GGA &Acc, GGB &B, GGB &Bcc) {
    herm_matrix_timestep_moving_view<typename GGC::scalar_type> ctmp(tstp, C);
    herm_matrix_timestep_moving_view<typename GGA::scalar_type> atmp(tstp, A);
    herm_matrix_timestep_moving_view<typename GGA::scalar_type> acctmp(tstp, Acc);
    herm_matrix_timestep_moving_view<typename GGB::scalar_type> btmp(tstp, B);
    herm_matrix_timestep_moving_view<typename GGB::scalar_type> bcctmp(tstp, Bcc);
    Bubble1_moving(tstp, ctmp, atmp, acctmp, btmp, bcctmp);
  }
  template <class GGC, class GGA, class GGB> void Bubble1_moving(int tstp, GGC &C, GGA &A, GGB &B) {
    herm_matrix_timestep_moving_view<typename GGC::scalar_type> ctmp(tstp, C);
    herm_matrix_timestep_moving_view<typename GGA::scalar_type> atmp(tstp, A);
    herm_matrix_timestep_moving_view<typename GGB::scalar_type> btmp(tstp, B);
    Bubble1_moving(tstp, ctmp, atmp, btmp);
  }

/// @private
  template<typename T>
  void get_bubble_2_timestep_moving(int tc,std::complex<T> *cret,std::complex<T> *cles,int sc,int c1,int c2,std::complex<T> *aret,std::complex<T> *ales,std::complex<T> *accret,std::complex<T> *accles,int sa,int a1,int a2,std::complex<T> *bret,std::complex<T> *bles,std::complex<T> *bccret,std::complex<T> *bccles,int sb,int b1,int b2){
    int m;
    int a12=a1*sa+a2,a21=a2*sa+a1,sa2=sa*sa;
    int b12=b1*sb+b2,b21=b2*sb+b1,sb2=sb*sb;
    int c12=c1*sc+c2,sc2=sc*sc;
    std::complex<T> ii=std::complex<T>(0,1.0);
    std::complex<T> bgtr12_tt1,agtr12_tt1,cgtr_tt1,cles_tt1;
    for(m=0;m<=tc;m++){
      bgtr12_tt1 = bret[m*sb2+b12] + bles[m*sb2+b12];
      agtr12_tt1 = aret[m*sa2+a12] + ales[m*sa2+a21];
      // Cgtr_{12}(tstp,m) = ii * Agtr_{12}(tstp,m)*Bles_{21}(m,tstp)
      cgtr_tt1 = ii * agtr12_tt1 * bgtr12_tt1;
      // Cles_{12}(tstp,m) = ii * Ales_{12}(tstp,m)*Bles_{12}(tstp,m)
      cles_tt1 = ii * ales[m*sa2+a21]*bles[m*sb2+b21];
      // Cret_{12}(tstp,m) = Cgtr_{12}(tstp,m) - Cles_{12}(tstp,m)
      cret[m*sc2+c12] = cgtr_tt1-cles_tt1;
      // Cles_{12}(m,tstp) = ii * Ales_{12}(m,tstp)*Bgtr_{12}(m,tstp)
      cles[m*sc2+c12] = ii* ales[m*sa2+a12] * bles[m*sb2+b12]; 
    }
  }
  /// @private
  template <typename T>
  void Bubble2_moving(int tstp, herm_matrix_timestep_moving_view<T> &C, int c1, int c2,
		      herm_matrix_timestep_moving_view<T> &A, herm_matrix_timestep_moving_view<T> &Acc, int a1,
		      int a2, herm_matrix_timestep_moving_view<T> &B, herm_matrix_timestep_moving_view<T> &Bcc,
		      int b1, int b2) {
    int tc=C.tc();
    // assert(ntau == A.ntau());
    // assert(ntau == B.ntau());
    // assert(ntau == Acc.ntau());
    // assert(ntau == Bcc.ntau());
    // assert(tstp == C.tstp_);
    // assert(tstp == B.tstp_);
    // assert(tstp == A.tstp_);
    // assert(tstp == Acc.tstp_);
    // assert(tstp == Bcc.tstp_);
    assert(a1 <= A.size1());
    assert(a2 <= A.size1());
    assert(b1 <= B.size1());
    assert(b2 <= B.size1());
    assert(c1 <= C.size1());
    assert(c2 <= C.size1());
    assert(Acc.size1() ==  A.size1());
    assert(Bcc.size1() == B.size1());
    get_bubble_2_timestep_moving(tc,C.retptr(0),C.lesptr(0),C.size1(),c1,c2,
  				 A.retptr(0),A.lesptr(0),Acc.retptr(0),Acc.lesptr(0),A.size1(),a1,a2,
  				 B.retptr(0),B.lesptr(0),Bcc.retptr(0),Bcc.lesptr(0),B.size1(),b1,b2);

  }
  /// @private
  template <typename T>
  void Bubble2_moving(int tstp, herm_matrix_timestep_moving_view<T> &C, int c1, int c2,
		      herm_matrix_timestep_moving_view<T> &A, int a1, int a2,
		      herm_matrix_timestep_moving_view<T> &B, int b1, int b2) {
    return Bubble2_moving(tstp, C, c1, c2, A, A, a1, a2, B, B, b1, b2);
  }
  /// @private
  template <typename T>
  void Bubble2_moving(int tstp, herm_matrix_timestep_moving_view<T> &C, herm_matrix_timestep_moving_view<T> &A,
		      herm_matrix_timestep_moving_view<T> &Acc, herm_matrix_timestep_moving_view<T> &B,
		      herm_matrix_timestep_moving_view<T> &Bcc) {
    return Bubble2_moving(tstp, C, 0, 0, A, Acc, 0, 0, B, Bcc, 0, 0);
  }
  /// @private
  template <typename T>
  void Bubble2_moving(int tstp, herm_matrix_timestep_moving_view<T> &C, herm_matrix_timestep_moving_view<T> &A,
		      herm_matrix_timestep_moving_view<T> &B) {
    return Bubble2_moving(tstp, C, A, A, B, B);
  }

  template <class GGC, class GGA, class GGB>
  void Bubble2_moving(int tstp, GGC &C, int c1, int c2, GGA &A, GGA &Acc, int a1, int a2, GGB &B,
		      GGB &Bcc, int b1, int b2) {
    herm_matrix_timestep_moving_view<typename GGC::scalar_type> ctmp(tstp, C);
    herm_matrix_timestep_moving_view<typename GGA::scalar_type> atmp(tstp, A);
    herm_matrix_timestep_moving_view<typename GGA::scalar_type> acctmp(tstp, Acc);
    herm_matrix_timestep_moving_view<typename GGB::scalar_type> btmp(tstp, B);
    herm_matrix_timestep_moving_view<typename GGB::scalar_type> bcctmp(tstp, Bcc);
    Bubble2_moving(tstp, ctmp, c1, c2, atmp, acctmp, a1, a2, btmp, bcctmp, b1, b2);
  }
  template <class GGC, class GGA, class GGB>
  void Bubble2_moving(int tstp, GGC &C, int c1, int c2, GGA &A, int a1, int a2, GGB &B, int b1,
		      int b2) {
    herm_matrix_timestep_moving_view<typename GGC::scalar_type> ctmp(tstp, C);
    herm_matrix_timestep_moving_view<typename GGA::scalar_type> atmp(tstp, A);
    herm_matrix_timestep_moving_view<typename GGB::scalar_type> btmp(tstp, B);
    Bubble2_moving(tstp, ctmp, c1, c2, atmp, a1, a2, btmp, b1, b2);
  }
  template <class GGC, class GGA, class GGB>
  void Bubble2_moving(int tstp, GGC &C, GGA &A, GGA &Acc, GGB &B, GGB &Bcc) {
    herm_matrix_timestep_moving_view<typename GGC::scalar_type> ctmp(tstp, C);
    herm_matrix_timestep_moving_view<typename GGA::scalar_type> atmp(tstp, A);
    herm_matrix_timestep_moving_view<typename GGA::scalar_type> acctmp(tstp, Acc);
    herm_matrix_timestep_moving_view<typename GGB::scalar_type> btmp(tstp, B);
    herm_matrix_timestep_moving_view<typename GGB::scalar_type> bcctmp(tstp, Bcc);
    Bubble2_moving(tstp, ctmp, atmp, acctmp, btmp, bcctmp);
  }
  template <class GGC, class GGA, class GGB> void Bubble2_moving(int tstp, GGC &C, GGA &A, GGB &B) {
    herm_matrix_timestep_moving_view<typename GGC::scalar_type> ctmp(tstp, C);
    herm_matrix_timestep_moving_view<typename GGA::scalar_type> atmp(tstp, A);
    herm_matrix_timestep_moving_view<typename GGB::scalar_type> btmp(tstp, B);
    Bubble2_moving(tstp, ctmp, atmp, btmp);
  }  
}//namespace cntr


#endif
