#include "cntr.hpp"

namespace cntr{
template <typename T> class moving_function;
template <typename T> class moving_herm;
template <typename T> class moving_herm_pseudo;
template <typename T> class moving_herm_timestep;


#define MOVING_HERM_ASSERT_LEVEL 1

/*//////////////////////////////////////////////////////////
moving_herm stores the Greenfunctions in a range

Gret(t,t-s), Gles(t,t-s) for t0>=t>=t0-nt,  0<=s<=tc

Storage order such that
G.ret_[t1] + t2*element_size_  -> Gret(t0-t1,t0-t1-t2)
G.les_[t1] + t2*element_size_  -> Gles(t0-t1,t0-t1-t2)

t0 is not stored
///////////////////////////////////////////////////////////*/
template <typename T> class moving_herm {
public:
    typedef std::complex<T> cplx;
    /* construction, destruction */
    moving_herm();
    ~moving_herm();
    moving_herm(int tc,int nt,int size1,int sig);
    //moving_herm(int tc,int size1,int sig);
    moving_herm(const moving_herm &g);
    moving_herm & operator=(const moving_herm &g);
    void clear(void);
    void resize(int tc,int nt,int size1);
    /* access size etc ... */
    int element_size(void) const{ return element_size_;}
    int size1(void) const{ return size1_;}
    int size2(void) const{ return size2_;}
    int tc(void) const{ return tc_;}
    int nt(void) const{ return nt_;}
    int sig(void) const{ return sig_;}
    // raw pointer to elements ... to be used with care
    inline cplx * lesptr(int i,int j){return les_[i] + j*element_size_;}  // points to Gret(t0-i,t0-i-j)
    inline cplx * retptr(int i,int j){return ret_[i] + j*element_size_;}  // points to Gles(t0-i,t0-i-j)
    // reading basic and derived elements to any Matrix type:
    // the get_... address time-arguments "relative to t0"  (i,j) == (t0-i,t0-i-j)
    // and work only for 0 <= i,j <= tc
    template<class Matrix> void get_les(int i,int j,Matrix &M);
    template<class Matrix> void get_gtr(int i,int j,Matrix &M);
    template<class Matrix> void get_ret(int i,int j,Matrix &M);
    // these will adress only (0,0) element for dim>1:
    inline void get_les(int i,int j,cplx &x);
    inline void get_gtr(int i,int j,cplx &x);
    inline void get_ret(int i,int j,cplx &x);
    //  -sig*ii*Gles(t0-i,t0-i)
    cplx density_matrix(int i);
    template<class Matrix> void density_matrix(int tstp,Matrix &M);
    // writing basic elements (also relative to t0)
    template<class Matrix> void set_les(int i,int j,Matrix &M);
    template<class Matrix> void set_ret(int i,int j,Matrix &M);
    inline void set_les(int i,int j,cplx x);
    inline void set_ret(int i,int j,cplx x);
    //FUNCTIONS BELOW OPERATE BY DEFAULT ON TIMESTEP 0
    // raw pointer to elements ... to be used with care
    inline cplx *lesptr(int j){return les_[0] + j*element_size_;}  // points to Gret(t0-i,t0-i-j)
    inline cplx *retptr(int j){return ret_[0] + j*element_size_;}  // points to Gles(t0-i,t0-i-j)
    template<class Matrix> void get_les(int j,Matrix &M){get_les(0,j,M);}
    template<class Matrix> void get_gtr(int j,Matrix &M){get_gtr(0,j,M);}
    template<class Matrix> void get_ret(int j,Matrix &M){get_ret(0,j,M);}
    inline void get_les(int j,cplx &x){ get_les(0,j,x);}
    inline void get_gtr(int j,cplx &x){ get_gtr(0,j,x);}
    inline void get_ret(int j,cplx &x){ get_ret(0,j,x);}
    cplx density_matrix(void){density_matrix(0);}
    template<class Matrix> void density_matrix(Matrix &M){density_matrix(0,M);}
    template<class Matrix> void set_les(int j,Matrix &M){set_les(0,j,M);}
    template<class Matrix> void set_ret(int j,Matrix &M){set_ret(0,j,M);}
    inline void set_les(int j,cplx x){set_les(0,j,x);}
    inline void set_ret(int j,cplx x){set_ret(0,j,x);}
    // INPUT/OUTPUT
    void print_to_file(const char *file,int precision=16);
    void read_from_file(const char *file);
    // DATA exchange with HERM_MATRIX
    // read data to slice i (relative to t0)
    void set_timestep(int i,herm_matrix_timestep_view<T> &g,herm_matrix_timestep_view<T> &gcc);
    void set_timestep(int i,int tstp,herm_matrix<T> &g,herm_matrix<T> &gcc);
    void set_timestep(herm_matrix_timestep_view<T> &g,herm_matrix_timestep_view<T> &gcc){set_timestep(0,g,gcc);}
    void set_timestep(herm_matrix<T> &g,herm_matrix<T> &gcc){set_timestep(0,g,gcc);}
    // Full read
    void set_from_G_backward(int tstp,herm_matrix<T> &g,herm_matrix<T> &gcc);
    void set_from_G_backward(int tstp,herm_matrix<T> &g);
    ////////
    void clear_timestep(int i);
    ///
    void clear_timestep(void){clear_timestep(0);}
    void get_timestep(int i,moving_herm_timestep<T> &g);
    void get_timestep(moving_herm_timestep<T> &g){get_timestep(0,g,0);}
    ///
    void set_timestep(int i,moving_herm<T> &g,int j);
    void set_timestep(int i,moving_herm<T> &g){set_timestep(i,g,0);}
    void set_timestep(int i,moving_herm_timestep<T> &g){set_timestep(i,g,0);}
    void set_timestep(moving_herm<T> &g){set_timestep(0,g,0);}
    void set_timestep(moving_herm_timestep<T> &g){set_timestep(0,g,0);}
    ///
    void incr_timestep(int i,moving_herm<T> &g,int j,cplx alpha);
    void incr_timestep(moving_herm<T> &g,int i,cplx alpha){incr_timestep(0,g,i,alpha);}
    void incr_timestep(int i,moving_herm<T> &g,cplx alpha){incr_timestep(i,g,0,alpha);}
    void incr_timestep(int i,moving_herm_timestep<T> &g,cplx alpha){incr_timestep(i,g,0,alpha);}
    void incr_timestep(moving_herm<T> &g,cplx alpha){incr_timestep(0,g,0,alpha);}
    void incr_timestep(moving_herm_timestep<T> &g,cplx alpha){incr_timestep(0,g,0,alpha);}
    // only for timestep 0
    void left_multiply(moving_function<T> &g,T weight);
    void right_multiply(moving_function<T> &g,T weight);
    void forward(void);
 private:
	   cplx* data_;
	   cplx** les_;
	   cplx** ret_;
       int tc_;
       int nt_;
	   int size1_;
	   int size2_;
	   int element_size_;
	   int sig_; // Bose = +1, Fermi =-1
};

// timestep is moving_herm with nt=0

template <typename T> class moving_herm_timestep : public moving_herm<T> {
public:
    typedef std::complex<T> cplx;
    /* construction, destruction */
    moving_herm_timestep() : moving_herm<T>::moving_herm() {};
    ~moving_herm_timestep(){};
    moving_herm_timestep(int tc,int size1,int sig) : moving_herm<T>::moving_herm(tc,0,size1,sig) {};
    moving_herm_timestep(const moving_herm_timestep &g) : moving_herm<T>::moving_herm(g) {};
    moving_herm_timestep & operator=(const moving_herm_timestep &g) {
        moving_herm<T>::operator=(g);
        return *this;
    }
    void resize(int tc,int size1){
        moving_herm<T>::resize(tc,0,size1);
    }
};
///////////////////////////////////////////////////////////*/
template <typename T> class moving_function {
public:
    typedef std::complex<T> cplx;
    /* construction, destruction */
    moving_function();
    ~moving_function();
    moving_function(int tc,int size1=1);
    moving_function(const moving_function &g);
    moving_function & operator=(const moving_function &g);
    void clear(void);
    void resize(int tc,int size1);
    /* access size etc ... */
    int element_size(void) const{ return element_size_;}
    int size1(void) const{ return size1_;}
    int size2(void) const{ return size2_;}
    int tc(void) const{ return tc_;}
    void set_t0(int t0);
    // raw pointer to elements ... to be used with care
    inline cplx * ptr(int i){return value_[i];}  // points to Gret(t0-i,t0-i-j)
    template<class EigenMatrix> void set_value(int i,EigenMatrix &M);
	template<class EigenMatrix> void get_value(int i,EigenMatrix &M);
	cplx & operator[](int i){return *ptr(i);} // useful only for size=1
	cplx & operator[](int i) const{return *ptr(i);} // useful only for size=1
    // INPUT/OUTPUT
    void print_to_file(const char *file,int precision=16);
    void read_from_file(const char *file);
    // DATA exchange with HERM_MATRIX
    void forward(void);
 private:
	   cplx* data_;
	   cplx** value_;
	   int tc_;
       int t0_;
	   int size1_;
	   int size2_;
	   int element_size_;
};
//#########################################################################################
// works with pseudo and usual ret, because it does not address gtr
template < typename T>
void vie2_timestep(moving_herm<T> &G,moving_herm<T> &F,moving_herm<T> &Fcc,moving_herm<T> &Q, integration::Integrator<T> &I, T h);
template < typename T>
void dyson_timestep(moving_herm<T> &G,moving_herm<T> &Sigma,moving_function<T> &eps,T mu, integration::Integrator<T> &I, T h);

template < typename T>
void vie2_pseudo_timestep(moving_herm<T> &G,moving_herm<T> &F,moving_herm<T> &Fcc,moving_herm<T> &Q, integration::Integrator<T> &I, T h){
vie2_timestep(G,F,Fcc,Q,I,h);
}
template < typename T>
void dyson_pseudo_timestep(moving_herm<T> &G,moving_herm<T> &Sigma,moving_function<T> &eps,T lam0, integration::Integrator<T> &I, T h){
    dyson_timestep(G,Sigma,eps,lam0,I,h);
}

//#########################################################################################
//  k-th order polynomila extrapolate of realtime functions (G^ret, G^vt, G^les)
template <typename T>
void extrapolate_timestep(moving_herm<T> &G,integration::Integrator<T> &I);
template <typename T>
T distance_norm2(int j1,moving_herm<T> &g1,int j2,moving_herm<T> &g2);
template <typename T>
T distance_norm2(int j,moving_herm<T> &g1,moving_herm_timestep<T> &g2){distance_norm2(j,g1,0,g2);}
template <typename T>
T distance_norm2(moving_herm_timestep<T> &g1,moving_herm_timestep<T> &g2){distance_norm2(0,g1,0,g2);}
template <typename T>
T distance_norm2(moving_herm_timestep<T> &g1,int j,moving_herm<T> &g2){distance_norm2(0,g1,j,g2);}


/*#########################################################################################
#
#   ...   useful routines to compute diagrams
#
#   BUBBLE1:  C_{c1,c2}(t1,t2) = ii * A_{a1,a2}(t1,t2) * B_{b2,b1}(t2,t1) 
#
#   BUBBLE2:  C_{c1,c2}(t1,t2) = ii * A_{a1,a2}(t1,t2) * B_{b1,b2}(t1,t2) 
#     
#   template GGA,GGB,GGC can be herm_matrix,herm_matrix_timestep,herm_matrix_timestep_view
#
#   (implementation in cntr_diagram_utilities.hpp, easy to add further calls with timestep 
#    instead of green functions)
#########################################################################################*/
// In versions where arguments (orbital indices, timespteps) are missing, they default to 0
//
template<typename T>
void MovBubble1(int tstp_C, moving_herm<T>  &C,int c1,int c2, int tstp_A,moving_herm<T>  &A,moving_herm<T>  &Acc,int a1,int a2, int tstp_B,moving_herm<T>  &B,moving_herm<T>  &Bcc,int b1,int b2);
template<typename T>
void MovBubble1(moving_herm<T>  &C,int c1,int c2,moving_herm<T>  &A, moving_herm<T>  &Acc,int a1,int a2,moving_herm<T>  &B, moving_herm<T>  &Bcc,int b1,int b2);
template<typename T>
void MovBubble1(moving_herm<T>  &C,int c1,int c2,moving_herm<T>  &A,int a1,int a2,moving_herm<T>  &B,int b1,int b2);
template<typename T>
void MovBubble1(moving_herm<T>  &C,moving_herm<T>  &A,moving_herm<T>  &B);
template<typename T>
void MovBubble1(moving_herm<T> &C,moving_herm<T>  &A, moving_herm<T>  &Acc,moving_herm<T>  &B, moving_herm<T>  &Bcc);

template<typename T>
void MovBubble2(int tstp_C, moving_herm<T>  &C,int c1,int c2, int tstp_A,moving_herm<T>  &A,moving_herm<T>  &Acc,int a1,int a2, int tstp_B,moving_herm<T>  &B,moving_herm<T>  &Bcc,int b1,int b2);
template<typename T>
void MovBubble2(moving_herm<T>  &C,int c1,int c2,moving_herm<T>  &A, moving_herm<T> &Acc,int a1,int a2,moving_herm<T>  &B, moving_herm<T>  &Bcc,int b1,int b2);
template<typename T>
void MovBubble2(moving_herm<T>  &C,int c1,int c2,moving_herm<T>  &A,int a1,int a2,moving_herm<T>  &B,int b1,int b2);
template<typename T>
void MovBubble2(moving_herm<T>  &C,moving_herm<T>  &A,moving_herm<T>  &B);
template<typename T>
void MovBubble2(moving_herm<T> &C,moving_herm<T>  &A, moving_herm<T>  &Acc,moving_herm<T>  &B, moving_herm<T>  &Bcc);


}
#include "moving_herm_member.hpp"
#include "moving_function_member.hpp"
#include "moving_herm_vie2.hpp"
#include "moving_herm_dyson.hpp"
#include "moving_herm_diagram_utilities.hpp"
#include "moving_herm_utils.hpp"




