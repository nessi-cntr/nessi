#ifndef CNTR_FULL_MATRIX_DECL
#define CNTR_FULL_MATRIX_DECL

#include "cntr_global_settings.hpp"
#include "cntr_herm_pseudo_decl.hpp"


namespace cntr {

template <typename T> class full_matrix {
 public:
	  typedef std::complex<T> cplx;
	  /* construction, destruction */
		  full_matrix();
		  ~full_matrix();
		  full_matrix(int nt,int ntau,int size1=1,int sig=-1);
		  full_matrix(const full_matrix &g);
                  full_matrix(const herm_matrix<T> &g);
		  full_matrix & operator=(const full_matrix &g);
	  /* resize */
		  void resize_discard(int nt,int ntau,int size1);
		  void resize_nt(int nt);
		  void resize(int nt,int ntau,int size1);
		  void clear(void);
	  /* access size etc ... */
		  int element_size(void){ return element_size_;}
		  int size1(void){ return size1_;}
		  int size2(void){ return size2_;}
		  int ntau(void){ return ntau_;}
		  int nt(void){ return nt_;}
		  int sig(void){ return sig_;}
		  int element_size(void) const{ return element_size_;}
		  int size1(void) const{ return size1_;}
		  int size2(void) const{ return size2_;}
		  int ntau(void)const{ return ntau_;}
		  int nt(void) const{ return nt_;}
		  int sig(void) const{ return sig_;}
		  void set_sig(int sig){sig_=sig;}
	  /* conversion from other types: pseudoparticle GF */
	  // raw pointer to elements ... to be used with care
		  inline cplx * lesptr(int i,int j);
		  inline cplx * gtrptr(int i,int j);
		  inline cplx * tvptr(int i,int j);
		  inline cplx * vtptr(int i,int j);
		  inline cplx * matptr(int i);
	  // reading basic elements to any Matrix type
		  template<class Matrix> void get_les(int i,int j,Matrix &M);
		  template<class Matrix> void get_gtr(int i,int j,Matrix &M);
		  template<class Matrix> void get_ret(int i,int j,Matrix &M);
		  template<class Matrix> void get_vt(int i,int j,Matrix &M);
		  template<class Matrix> void get_tv(int i,int j,Matrix &M);
		  template<class Matrix> void get_mat(int i,Matrix &M);
		  template<class Matrix> void get_matminus(int i,Matrix &M);
      // reading complex numbers:
	  // these will adress only (0,0) element for dim>1:
		  inline void get_les(int i,int j,cplx &x);
		  inline void get_gtr(int i,int j,cplx &x);
		  inline void get_ret(int i,int j,cplx &x);
		  inline void get_vt(int i,int j,cplx &x);
		  inline void get_tv(int i,int j,cplx &x);
		  inline void get_mat(int i,cplx &x);
		  inline void get_matminus(int i,cplx &x);
	  // writing basic elements:
		  template<class Matrix> void set_les(int i,int j,Matrix &M);
		  template<class Matrix> void set_gtr(int i,int j,Matrix &M);
		  template<class Matrix> void set_vt(int i,int j,Matrix &M);
		  template<class Matrix> void set_tv(int i,int j,Matrix &M);
		  template<class Matrix> void set_mat(int i,Matrix &M);
		  inline void set_les(int i,int j,cplx x);
		  inline void set_gtr(int i,int j,cplx x);
		  inline void set_vt(int i,int j,cplx x);
		  inline void set_tv(int i,int j,cplx x);
		  inline void set_mat(int i,cplx x);
	  // INPUT/OUTPUT
	      void print_to_file(const char *file,int precision=16);
	      void read_from_file(const char *file);
	  // simple operations on the timesteps (n<0: Matsubara)
		  void set_timestep_zero(int tstp);
		  void set_timestep(int tstp,full_matrix<T> &g1);
		  void set_timestep(int tstp,herm_matrix<T> &g1,herm_matrix<T> &gcc);
		  void set_timestep(int tstp,herm_matrix<T> &g1){set_timestep(tstp,g1,g1);}
		  void set_timestep(int tstp,herm_pseudo<T> &g1,herm_pseudo<T> &gcc);
		  void set_timestep(int tstp,herm_pseudo<T> &g1){set_timestep(tstp,g1,g1);}
		  // assumes hermital symmetry:
		  void set_timestep(int tstp,herm_matrix_timestep<T> &g,herm_matrix_timestep<T> &gcc);
		  void set_timestep(int tstp,herm_matrix_timestep<T> &g1){set_timestep(tstp,g1,g1);}
		  void get_timestep(int tstp,herm_matrix_timestep<T> &g1);
	 private:
	   cplx* gtr_;
	   cplx* les_;
	   cplx* tv_;
	   cplx* vt_;
	   cplx* mat_;
	   int nt_;
	   int ntau_;
	   int size1_;
	   int size2_;
	   int element_size_;
	   int sig_; // Bose = +1, Fermi =-1
};

} //namespace cntr

#endif  // CNTR_FULL_MATRIX_DECL
