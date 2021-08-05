#ifndef CNTR_FULL_MATRIX_IMPL
#define CNTR_FULL_MATRIX_IMPL

#include "cntr_full_member_decl.hpp"
//#include "cntr_exception.hpp"
//#include "cntr_elements.hpp"
//#include "cntr_function_decl.hpp"

namespace cntr {
/*
THIS FILE IS PART OF cntr.hpp, NAMESPACE cntr
CONTENT: MEMBER FUNCTIONS OF cntr::full_matrix
*/
/* #######################################################################################
#
#   CONSTRUCTION/DESTRUCTION
#
########################################################################################*/
template <typename T> full_matrix<T>::full_matrix(){
   gtr_=0;
   les_=0;
   vt_=0;
   tv_=0;
   mat_=0;
   ntau_=0;
   nt_=0;
   size1_=0;
   size2_=0;
   element_size_=0;
   sig_=-1;
}
template <typename T> full_matrix<T>::~full_matrix(){
   delete [] gtr_;
   delete [] les_;
   delete [] vt_;
   delete [] tv_;
   delete [] mat_;
}
template <typename T> full_matrix<T>::full_matrix(int nt,int ntau,int size1,int sig){
   assert(size1>=0 && nt>=-1 && sig*sig==1 && ntau>=0);
   nt_=nt;
   ntau_=ntau;
   sig_=sig;
   size1_=size1;
   size2_=size1;
   element_size_=size1*size1;
   if(size1>0){
	   mat_ = new cplx [(ntau_+1)*element_size_];
	   memset(mat_, 0, sizeof(cplx)*(ntau_+1)*element_size_);
   }else{
	   mat_=0;
   }
   if(nt>=0 && size1>0){
	   gtr_ = new cplx [(nt_+1)*(nt_+1)*element_size_];
	   les_ = new cplx [(nt_+1)*(nt_+1)*element_size_];
	   tv_ = new cplx [(nt_+1)*(ntau_+1)*element_size_];
	   vt_ = new cplx [(nt_+1)*(ntau_+1)*element_size_];
	   memset(gtr_, 0, sizeof(cplx)*(nt_+1)*(nt_+1)*element_size_);
	   memset(les_, 0, sizeof(cplx)*(nt_+1)*(nt_+1)*element_size_);
	   memset(vt_, 0, sizeof(cplx)*(nt_+1)*(ntau_+1)*element_size_);
	   memset(tv_, 0, sizeof(cplx)*(nt_+1)*(ntau_+1)*element_size_);
   }else{
	   gtr_=0;
	   les_=0;
	   vt_=0;
	   tv_=0;
   }
}
template <typename T> unsigned long RAM_full_matrix(int nt,int ntau,int dim){
    size_t storage=0;
 	int element_size=dim*dim;
 	storage+=sizeof(typename full_matrix<T>::cplx)*(ntau+1)*element_size;
 	storage+=sizeof(typename full_matrix<T>::cplx)*(nt+1)*(nt+1)*element_size;
 	storage+=sizeof(typename full_matrix<T>::cplx)*(nt+1)*(nt+1)*element_size;
 	storage+=sizeof(typename full_matrix<T>::cplx)*(nt+1)*(ntau+1)*element_size;
	storage+=sizeof(typename full_matrix<T>::cplx)*(nt+1)*(ntau+1)*element_size;
    return storage;
}
template <typename T> full_matrix<T>::full_matrix(const full_matrix &g){
   nt_=g.nt_;
   ntau_=g.ntau_;
   sig_=g.sig_;
   size1_=g.size1_;
   size2_=g.size1_;
   element_size_=size1_*size1_;
   if(size1_>0){
	   mat_ = new cplx [(ntau_+1)*element_size_];
	   memcpy(mat_, g.mat_, sizeof(cplx)*(ntau_+1)*element_size_);
   }else{
	   mat_=0;
   }
   if(nt_>=0 && size1_>0){
	   gtr_ = new cplx [(nt_+1)*(nt_+1)*element_size_];
	   les_ = new cplx [(nt_+1)*(nt_+1)*element_size_];
	   tv_ = new cplx [(nt_+1)*(ntau_+1)*element_size_];
	   vt_ = new cplx [(nt_+1)*(ntau_+1)*element_size_];
	   memcpy(gtr_, g.gtr_, sizeof(cplx)*(nt_+1)*(nt_+1)*element_size_);
	   memcpy(les_, g.les_, sizeof(cplx)*(nt_+1)*(nt_+1)*element_size_);
	   memcpy(vt_, g.vt_, sizeof(cplx)*(nt_+1)*(ntau_+1)*element_size_);
	   memcpy(tv_, g.tv_, sizeof(cplx)*(nt_+1)*(ntau_+1)*element_size_);
   }else{
	   gtr_=0;
	   les_=0;
	   vt_=0;
	   tv_=0;
   }
}
template <typename T> full_matrix<T>::full_matrix(const herm_matrix<T> &g){

  typedef Eigen::Matrix<std::complex<T>,
			Eigen::Dynamic, Eigen::Dynamic> matrix_type;

  nt_=g.nt();
  ntau_=g.ntau();
  sig_=g.sig();
  size1_=g.size1();
  size2_=g.size1();
  element_size_=size1_*size1_;
   if(size1_>0){
	   mat_ = new cplx [(ntau_+1)*element_size_];
	   //memcpy(mat_, g.mat_, sizeof(cplx)*(ntau_+1)*element_size_);

	   for(int tau = 0; tau <= ntau_; tau++) {
	     matrix_type mat;
	     g.get_mat(tau, mat);
	     this->set_mat(tau, mat);
	   }

   }else{
	   mat_=0;
   }
   if(nt_>=0 && size1_>0){
	   gtr_ = new cplx [(nt_+1)*(nt_+1)*element_size_];
	   les_ = new cplx [(nt_+1)*(nt_+1)*element_size_];
	   tv_ = new cplx [(nt_+1)*(ntau_+1)*element_size_];
	   vt_ = new cplx [(nt_+1)*(ntau_+1)*element_size_];
	   //memcpy(gtr_, g.gtr_, sizeof(cplx)*(nt_+1)*(nt_+1)*element_size_);
	   //memcpy(les_, g.les_, sizeof(cplx)*(nt_+1)*(nt_+1)*element_size_);
	   //memcpy(vt_, g.vt_, sizeof(cplx)*(nt_+1)*(ntau_+1)*element_size_);
	   //memcpy(tv_, g.tv_, sizeof(cplx)*(nt_+1)*(ntau_+1)*element_size_);


	   for(int t1 = 0; t1 <= nt_; t1++) {
	   for(int t2 = 0; t2 <= nt_; t2++) {
	     matrix_type les;
	     g.get_les(t1, t2, les);
	     this->set_les(t1, t2, les);

	     matrix_type gtr;
	     g.get_gtr(t1, t2, gtr);
	     this->set_gtr(t1, t2, gtr);

	   }}

	   for(int t = 0; t <= nt_; t++) {
	   for(int tau = 0; tau <= ntau_; tau++) {
	     matrix_type vt;
	     g.get_vt(tau, t, vt);
	     this->set_vt(tau, t, vt);
	   }}

   }else{
	   gtr_=0;
	   les_=0;
	   vt_=0;
	   tv_=0;
   }
}
template <typename T> full_matrix<T> & full_matrix<T>::operator=(const full_matrix &g){
	if(this==&g) return *this;
	sig_=g.sig_;
	if( nt_!=g.nt_ || ntau_!=g.ntau_ || size1_!=g.size1_){
	   delete [] gtr_;
	   delete [] les_;
	   delete [] vt_;
	   delete [] tv_;
	   delete [] mat_;
	   nt_=g.nt_;
	   ntau_=g.ntau_;
	   size1_=g.size1_;
	   size2_=g.size1_;
	   element_size_=size1_*size1_;
	   if(size1_>0){
		   mat_ = new cplx [(ntau_+1)*element_size_];
	   }else{
		   mat_=0;
	   }
	   if(size1_>0 && nt_>=0){
		   gtr_ = new cplx [(nt_+1)*(nt_+1)*element_size_];
		   les_ = new cplx [(nt_+1)*(nt_+1)*element_size_];
		   tv_ = new cplx [(nt_+1)*(ntau_+1)*element_size_];
		   vt_ = new cplx [(nt_+1)*(ntau_+1)*element_size_];
	   }else{
		   gtr_=0;
		   les_=0;
		   vt_=0;
		   tv_=0;
	   }
	}
	if(size1_>0){
		   memcpy(mat_, g.mat_, sizeof(cplx)*(ntau_+1)*element_size_);
		   if(nt_>=0){
			   memcpy(gtr_, g.gtr_, sizeof(cplx)*(nt_+1)*(nt_+1)*element_size_);
			   memcpy(les_, g.les_, sizeof(cplx)*(nt_+1)*(nt_+1)*element_size_);
			   memcpy(vt_, g.vt_, sizeof(cplx)*(nt_+1)*(ntau_+1)*element_size_);
			   memcpy(tv_, g.tv_, sizeof(cplx)*(nt_+1)*(ntau_+1)*element_size_);
		   }
	}
	return *this;
}
/* #######################################################################################
#
#   RESIZE FUNCTIONS
#
########################################################################################*/
template <typename T> void full_matrix<T>::resize_discard(int nt,int ntau,int size1){
   assert(ntau>=0 && nt>=-1 && size1>=0);
   delete [] gtr_;
   delete [] les_;
   delete [] vt_;
   delete [] tv_;
   delete [] mat_;
   nt_=nt;
   ntau_=ntau;
   size1_=size1;
   size2_=size1;
   element_size_=size1*size1;
   if(size1>0){
	   mat_ = new cplx [(ntau_+1)*element_size_];
	   memset(mat_, 0, sizeof(cplx)*(ntau_+1)*element_size_);
   }else{
	   mat_=0;
   }
   if(nt_>=0 && size1_>0){
	   gtr_ = new cplx [(nt_+1)*(nt_+1)*element_size_];
	   les_ = new cplx [(nt_+1)*(nt_+1)*element_size_];
	   tv_ = new cplx [(nt_+1)*(ntau_+1)*element_size_];
	   vt_ = new cplx [(nt_+1)*(ntau_+1)*element_size_];
	   memset(gtr_, 0, sizeof(cplx)*(nt_+1)*(nt_+1)*element_size_);
	   memset(les_, 0, sizeof(cplx)*(nt_+1)*(nt_+1)*element_size_);
	   memset(vt_, 0, sizeof(cplx)*(nt_+1)*(ntau_+1)*element_size_);
	   memset(tv_, 0, sizeof(cplx)*(nt_+1)*(ntau_+1)*element_size_);
   }
   else{
	   gtr_=0;
	   les_=0;
	   vt_=0;
	   tv_=0;
   }
}
template <typename T> void full_matrix<T>::resize_nt(int nt){
   int ntold=nt_,i,j,l,nt1=(nt_>nt ? nt : nt_);
   cplx *gtr,*les,*vt,*tv;
   assert(nt>=-1);
   nt_=nt;
   if(size1_==0) return;
   if(nt_>=0){
	   gtr = new cplx [(nt_+1)*(nt_+1)*element_size_];
	   les = new cplx [(nt_+1)*(nt_+1)*element_size_];
	   tv = new cplx [(nt_+1)*(ntau_+1)*element_size_];
	   vt = new cplx [(nt_+1)*(ntau_+1)*element_size_];
	   memset(gtr, 0, sizeof(cplx)*(nt_+1)*(nt_+1)*element_size_);
	   memset(les, 0, sizeof(cplx)*(nt_+1)*(nt_+1)*element_size_);
	   memset(vt, 0, sizeof(cplx)*(nt_+1)*(ntau_+1)*element_size_);
	   memset(tv, 0, sizeof(cplx)*(nt_+1)*(ntau_+1)*element_size_);
	   for(i=0;i<=nt1;i++){
		  for(j=0;j<=nt1;j++){
			 for(l=0;l<element_size_;l++){
				gtr[element_size_*(i*(nt_+1)+j)+l]=gtr_[element_size_*(i*(ntold+1)+j)+l];
				les[element_size_*(i*(nt_+1)+j)+l]=les_[element_size_*(i*(ntold+1)+j)+l];
			 }
		  }
		  for(j=0;j<=ntau_;j++){
			 for(l=0;l<element_size_;l++){
				tv[element_size_*(i*(ntau_+1)+j)+l]=tv_[element_size_*(i*(ntau_+1)+j)+l];
				vt[element_size_*(j*(nt_+1)+i)+l]=vt_[element_size_*(j*(ntold+1)+i)+l];
			 }
		  }
	   }
   }else{
	   gtr=0;
	   les=0;
	   vt=0;
	   tv=0;
   }
   delete [] gtr_;
   delete [] les_;
   delete [] vt_;
   delete [] tv_;
   gtr_=gtr;
   les_=les;
   tv_=tv;
   vt_=vt;
}
template <typename T> void full_matrix<T>::resize(int nt,int ntau,int size1){
  //std::cout << "full_matrix<T>::" << __FUNCTION__ << " " << nt << " " << ntau << " " << size1 << std::endl;
  if(ntau==ntau_ && size1_==size1) resize_nt(nt);
  else resize_discard(nt,ntau,size1);
}
template <typename T> void full_matrix<T>::clear(void){
	if(size1_==0) return;
	memset(mat_, 0, sizeof(cplx)*(ntau_+1)*element_size_);
	if(nt_>=0){
		memset(gtr_, 0, sizeof(cplx)*(nt_+1)*(nt_+1)*element_size_);
		memset(les_, 0, sizeof(cplx)*(nt_+1)*(nt_+1)*element_size_);
		memset(vt_, 0, sizeof(cplx)*(nt_+1)*(ntau_+1)*element_size_);
		memset(tv_, 0, sizeof(cplx)*(nt_+1)*(ntau_+1)*element_size_);
	}
}
/* #######################################################################################
#
#   RAW POINTERS TO ELEMENTS
#
########################################################################################*/
template<typename T> inline std::complex<T> * full_matrix<T>::gtrptr(int t,int t1){
	assert( t>=0 && t1>=0 && t<=nt_ && t1<=nt_);
	return gtr_ + (t*(nt_+1)+t1)*element_size_;
}
template<typename T> inline std::complex<T> * full_matrix<T>::lesptr(int t,int t1){
	assert( t>=0 && t1>=0 && t<=nt_ && t1<=nt_);
	return les_ + (t*(nt_+1)+t1)*element_size_;
}
template<typename T> inline std::complex<T> * full_matrix<T>::tvptr(int t,int tau){
	assert( t>=0 && tau>=0 && t<=nt_ && tau<=ntau_);
	return tv_ + (t*(ntau_+1)+tau)*element_size_;
}
template<typename T> inline std::complex<T> * full_matrix<T>::vtptr(int tau,int t){
	assert( t>=0 && tau>=0 && tau<=ntau_ && t<=nt_);
	return vt_ + (tau*(nt_+1)+t)*element_size_;
}
template<typename T> inline std::complex<T> * full_matrix<T>::matptr(int tau){
	assert( tau>=0 && tau<=ntau_);
	return mat_ + tau*element_size_;
}
/* #######################################################################################
#
#   READING ELEMENTS TO ANY MATRIX TYPE
#   OR TO COMPLEX NUMBERS (then only the (0,0) element is addressed for dim>0)
#
########################################################################################*/
#define FULL_MATRIX_READ_ELEMENT_MATRIX {int r,s,dim=size1_;M.resize(dim,dim);for(r=0;r<dim;r++) for(s=0;s<dim;s++) M(r,s)=x[r*dim+s];}
template<typename T> template <class Matrix> void full_matrix<T>::get_gtr(int i,int j,Matrix &M){
  cplx * x=gtrptr(i,j);
  FULL_MATRIX_READ_ELEMENT_MATRIX
}
template<typename T> template <class Matrix> void full_matrix<T>::get_les(int i,int j,Matrix &M){
  cplx * x=lesptr(i,j);
  FULL_MATRIX_READ_ELEMENT_MATRIX
}
template<typename T> template <class Matrix> void full_matrix<T>::get_ret(int i,int j,Matrix &M){
  Matrix M1;
  get_gtr(i,j,M);
  get_les(i,j,M1);
  M -= M1;
}
template<typename T> template <class Matrix> void full_matrix<T>::get_vt(int i,int j,Matrix &M){
  cplx * x=vtptr(i,j);
  FULL_MATRIX_READ_ELEMENT_MATRIX
}
template<typename T> template <class Matrix> void full_matrix<T>::get_tv(int i,int j,Matrix &M){
  cplx * x=tvptr(i,j);
  FULL_MATRIX_READ_ELEMENT_MATRIX
}
template<typename T> template <class Matrix> void full_matrix<T>::get_mat(int i,Matrix &M){
  cplx * x=matptr(i);
  FULL_MATRIX_READ_ELEMENT_MATRIX
}
template<typename T> template <class Matrix> void full_matrix<T>::get_matminus(int i,Matrix &M){
  get_mat(ntau_-i,M);
  if(sig_==-1) M=-M;
}
template<typename T> inline void full_matrix<T>::get_gtr(int i,int j,cplx &x){ x=*gtrptr(i,j); }
template<typename T> inline void full_matrix<T>::get_les(int i,int j,cplx &x){ x=*lesptr(i,j); }
template<typename T> inline void full_matrix<T>::get_ret(int i,int j,cplx &x){ x= *gtrptr(i,j)-(*lesptr(i,j)); }
template<typename T> inline void full_matrix<T>::get_vt(int i,int j,cplx &x){ x=*vtptr(i,j); }
template<typename T> inline void full_matrix<T>::get_tv(int i,int j,cplx &x){ x=*tvptr(i,j); }
template<typename T> inline void full_matrix<T>::get_mat(int i,cplx &x){ x=*matptr(i); }
template<typename T> inline void full_matrix<T>::get_matminus(int i,cplx &x){ x=(1.0*sig_)*(*matptr(ntau_-i)); }
/* #######################################################################################
#
#   WRITING ELEMENTS FROM ANY MATRIX TYPE
#   OR FROM COMPLEX NUMBERS (then only the (0,0) element is addressed for dim>0)
#
########################################################################################*/
//#define FULL_MATRIX_SET_ELEMENT_MATRIX {int r,s,dim=size1_;assert(M.size1()==dim && M.size2()==dim);for(r=0;r<dim;r++) for(s=0;s<dim;s++) x[r*dim+s]=M(r,s);}
#define FULL_MATRIX_SET_ELEMENT_MATRIX {int r,s,dim=size1_;assert(M.rows()==dim && M.cols()==dim);for(r=0;r<dim;r++) for(s=0;s<dim;s++) x[r*dim+s]=M(r,s);}
template<typename T> template <class Matrix> void full_matrix<T>::set_les(int i,int j,Matrix &M){
  cplx * x=lesptr(i,j);
  FULL_MATRIX_SET_ELEMENT_MATRIX
}
template<typename T> template <class Matrix> void full_matrix<T>::set_gtr(int i,int j,Matrix &M){
  cplx * x=gtrptr(i,j);
  FULL_MATRIX_SET_ELEMENT_MATRIX
}
template<typename T> template <class Matrix> void full_matrix<T>::set_vt(int i,int j,Matrix &M){
  cplx * x=vtptr(i,j);
  FULL_MATRIX_SET_ELEMENT_MATRIX
}
template<typename T> template <class Matrix> void full_matrix<T>::set_tv(int i,int j,Matrix &M){
  cplx * x=tvptr(i,j);
  FULL_MATRIX_SET_ELEMENT_MATRIX
}
template<typename T> template <class Matrix> void full_matrix<T>::set_mat(int i,Matrix &M){
  cplx * x=matptr(i);
  FULL_MATRIX_SET_ELEMENT_MATRIX
}
template<typename T> inline void full_matrix<T>::set_les(int i,int j,cplx x){ *lesptr(i,j)=x; }
template<typename T> inline void full_matrix<T>::set_gtr(int i,int j,cplx x){ *gtrptr(i,j)=x; }
template<typename T> inline void full_matrix<T>::set_vt(int i,int j,cplx x){ *vtptr(i,j)=x; }
template<typename T> inline void full_matrix<T>::set_tv(int i,int j,cplx x){ *tvptr(i,j)=x; }
template<typename T> inline void full_matrix<T>::set_mat(int i,cplx x){ *matptr(i)=x; }
/* #######################################################################################
#
#   INPUT/OUTPUT FROM/TO FILES
#
########################################################################################*/
template<typename T> void full_matrix<T>::print_to_file(const char *file,int precision){
		int i,j,l,sg=element_size_;
		std::ofstream out;
		out.open(file,std::ios::out);
		out.precision(precision);
		out << "# " << nt_ << " " << ntau_ << " " << size1_ << " " << " " << sig_ << std::endl;
		for(j=0;j<=ntau_;j++){
		   out << "mat: " << j ;
		   for(l=0;l<sg;l++) out << " " << matptr(j)[l].real() << " " << matptr(j)[l].imag();
		   out << std::endl;
	    }
	    out << std::endl;
		if(nt_>=0){
			for(i=0;i<=nt_;i++){
			   for(j=0;j<=nt_;j++){
				  out << "gtr: " << i << " " << j;
				  for(l=0;l<sg;l++) out << " " << gtrptr(i,j)[l].real() << " " << gtrptr(i,j)[l].imag();
				  out << std::endl;
			   }
			   out << std::endl;
			}
			out << std::endl;
			for(i=0;i<=nt_;i++){
			   for(j=0;j<=nt_;j++){
				  out << "les: " << i << " " << j;
				  for(l=0;l<sg;l++) out << " " << lesptr(i,j)[l].real() << " " << lesptr(i,j)[l].imag();
				  out << std::endl;
			   }
			   out << std::endl;
			}
			out << std::endl;
			for(i=0;i<=nt_;i++){
			   for(j=0;j<=ntau_;j++){
				  out << "tv: " << i << " " << j;
				  for(l=0;l<sg;l++) out << " " << tvptr(i,j)[l].real() << " " << tvptr(i,j)[l].imag();
				  out << std::endl;
			   }
			   out << std::endl;
			}
			out << std::endl;
			for(i=0;i<=ntau_;i++){
			   for(j=0;j<=nt_;j++){
				  out << "vt: " << i << " " << j;
				  for(l=0;l<sg;l++) out << " " << vtptr(i,j)[l].real() << " " << vtptr(i,j)[l].imag();
				  out << std::endl;
			   }
			   out << std::endl;
			}
			out << std::endl;
        }
		out.close();
}
template< typename T> void full_matrix<T>::read_from_file(const char *file){
		int i,n,m,j,l,size1,sg,sig;
                double real, imag;
		std::string s;
		std::ifstream out;
		out.open(file,std::ios::in);
		if(!(out >> s >> n >> m >> size1 >> sig)){
		  std::cerr << "read G from file " << file << " error in file" << std::endl;
		  abort();
		}
		if(n>nt_ || m!=ntau_ || size1!=size1_) resize(n,m,size1);
		sig_=sig;
		sg=element_size_;
		for(j=0;j<=ntau_;j++){
		   out >> s >> s ;
		   for(l=0;l<sg;l++){
			 if(!( out >> real >> imag )){
			   std::cerr << "read G from file " << file << " error at mat (" << j << ")"<< std::endl;
			   abort();
			 }
			 matptr(j)[l] = std::complex<T>(real, imag);
		   }
		}
		if(n>=0){
			for(i=0;i<=n;i++){
			   for(j=0;j<=n;j++){
				   out >> s >> s >> s ;
				   for(l=0;l<sg;l++){
					 if(!( out >> real >> imag )){
					   std::cerr << "read G from file " << file << " error at gtr (" << i<< "," << j << ")"<< std::endl;
					   abort();
					 }
					 gtrptr(i,j)[l] = std::complex<T>(real, imag);
				   }
			   }
			}
			for(i=0;i<=n;i++){
			   for(j=0;j<=n;j++){
				   out >> s >> s >> s ;
				   for(l=0;l<sg;l++){
					 if(!( out >> real >> imag )){
					   std::cerr << "read G from file " << file << " error at les (" << i<< "," << j << ")"<< std::endl;
					   abort();
					 }
					 lesptr(i,j)[l] = std::complex<T>(real, imag);
				   }
			   }
			}
			for(i=0;i<=n;i++){
			   for(j=0;j<=ntau_;j++){
				   out >> s >> s >> s ;
				   for(l=0;l<sg;l++){
					 if(!( out >> real >> imag )){
					   std::cerr << "read G from file " << file << " error at tv (" << i<< "," << j << ")"<< std::endl;
					   abort();
					 }
					 tvptr(i,j)[l] = std::complex<T>(real, imag);
				   }
			   }
			}
			for(i=0;i<=ntau_;i++){
			   for(j=0;j<=n;j++){
				   out >> s >> s >> s ;
				   for(l=0;l<sg;l++){
					 if(!( out >> real >> imag )){
					   std::cerr << "read G from file " << file << " error at vt (" << i<< "," << j << ")"<< std::endl;
					   abort();
					 }
					 vtptr(i,j)[l] = std::complex<T>(real, imag);
				   }
			   }
			}
		}
		out.close();
}
/* #######################################################################################
#
#   SIMPLE OPERATIONS ON TIMESTEPS
#   NOTE: tstp IS A PHYSICAL TIME, tstp=-1 is the matsubara branch
#
########################################################################################*/
template<typename T> void full_matrix<T>::set_timestep_zero(int tstp){
   int i;
   assert(tstp>=-1 && tstp<=nt_);
   if(tstp==-1){
	  memset(matptr(0),0,sizeof(cplx)*(ntau_+1)*element_size_);
   }else{
      // the lines
	  memset(tvptr(tstp,0),0,sizeof(cplx)*(ntau_+1)*element_size_);
      memset(lesptr(tstp,0),0,sizeof(cplx)*(tstp+1)*element_size_);
	  memset(gtrptr(tstp,0),0,sizeof(cplx)*(tstp+1)*element_size_);
	  // and the columns
	  for(i=0;i<=ntau_;i++) element_set_zero<T,LARGESIZE>(size1_,vtptr(i,tstp));
	  for(i=0;i<tstp;i++){
		  element_set_zero<T,LARGESIZE>(size1_,lesptr(i,tstp));
		  element_set_zero<T,LARGESIZE>(size1_,gtrptr(i,tstp));
	  }
   }
}
template<typename T> void full_matrix<T>::set_timestep(int tstp,full_matrix<T> &g1){
   int i;
   assert(tstp>=-1 && tstp<=nt_ && tstp<=g1.nt());
   assert(g1.size1()==size1_);
   assert(g1.ntau()==ntau_);
   if(tstp==-1){
	  memcpy(matptr(0),g1.matptr(0),sizeof(cplx)*(ntau_+1)*element_size_);
   }else{
      // the lines
	  memcpy(tvptr(tstp,0),g1.tvptr(tstp,0),sizeof(cplx)*(ntau_+1)*element_size_);
      memcpy(lesptr(tstp,0),g1.lesptr(tstp,0),sizeof(cplx)*(tstp+1)*element_size_);
	  memcpy(gtrptr(tstp,0),g1.gtrptr(tstp,0),sizeof(cplx)*(tstp+1)*element_size_);
	  // and the columns
	  for(i=0;i<=ntau_;i++) element_set<T,LARGESIZE>(size1_,vtptr(i,tstp),g1.vtptr(i,tstp));
	  for(i=0;i<tstp;i++){
		  element_set<T,LARGESIZE>(size1_,lesptr(i,tstp),g1.lesptr(i,tstp));
		  element_set<T,LARGESIZE>(size1_,gtrptr(i,tstp),g1.gtrptr(i,tstp));
	  }
   }
}
template<typename T> void full_matrix<T>::set_timestep(int tstp,herm_matrix<T> &g,herm_matrix<T> &gcc){
   int i;
   cplx *ret_ti = new cplx [size1_*size1_];
   cplx *adv_it = new cplx [size1_*size1_];
   cplx *les_ti = new cplx [size1_*size1_];
   cplx *les_it = new cplx [size1_*size1_];
   assert(tstp>=-1 && tstp<=nt_ && tstp<=g.nt() && tstp<=gcc.nt());
   assert(g.size1()==size1_);
   assert(gcc.size1()==size1_);
   assert(g.ntau()==ntau_);
   assert(gcc.ntau()==ntau_);
   if(tstp==-1){
	  memcpy(matptr(0),g.matptr(0),sizeof(cplx)*(ntau_+1)*element_size_);
   }else{
	  memcpy(tvptr(tstp,0),g.tvptr(tstp,0),sizeof(cplx)*(ntau_+1)*element_size_);
      if(size1_==1){
			for(i=0;i<=tstp;i++){
			element_set<T,1>(size1_,ret_ti,g.retptr(tstp,i));
			element_set<T,1>(size1_,les_it,g.lesptr(i,tstp));
			element_minusconj<T,1>(size1_,les_ti,gcc.lesptr(i,tstp));
			element_minusconj<T,1>(size1_,adv_it,gcc.retptr(tstp,i)); // -adv
			element_set<T,1>(size1_,lesptr(tstp,i),les_ti);
			element_set<T,1>(size1_,lesptr(i,tstp),les_it);
			element_incr<T,1>(size1_,ret_ti,les_ti);
			element_incr<T,1>(size1_,adv_it,les_it);
			element_set<T,1>(size1_,gtrptr(tstp,i),ret_ti);
			element_set<T,1>(size1_,gtrptr(i,tstp),adv_it);
		  }
		  if(sig_==-1){
			  for(i=0;i<=ntau_;i++) element_conj<T,1>(size1_,vtptr(i,tstp),gcc.tvptr(tstp,ntau_-i));
		  }else{
			  for(i=0;i<=ntau_;i++) element_minusconj<T,1>(size1_,vtptr(i,tstp),gcc.tvptr(tstp,ntau_-i));
		  }
	  }else{
		  for(i=0;i<=tstp;i++){
			element_set<T,LARGESIZE>(size1_,ret_ti,g.retptr(tstp,i));
			element_set<T,LARGESIZE>(size1_,les_it,g.lesptr(i,tstp));
			element_minusconj<T,LARGESIZE>(size1_,les_ti,gcc.lesptr(i,tstp));
			element_minusconj<T,LARGESIZE>(size1_,adv_it,gcc.retptr(tstp,i)); // -adv
			element_set<T,LARGESIZE>(size1_,lesptr(tstp,i),les_ti);
			element_set<T,LARGESIZE>(size1_,lesptr(i,tstp),les_it);
			element_incr<T,LARGESIZE>(size1_,ret_ti,les_ti);
			element_incr<T,LARGESIZE>(size1_,adv_it,les_it);
			element_set<T,LARGESIZE>(size1_,gtrptr(tstp,i),ret_ti);
			element_set<T,LARGESIZE>(size1_,gtrptr(i,tstp),adv_it);
		  }
		  if(sig_==-1){
			  for(i=0;i<=ntau_;i++) element_conj<T,LARGESIZE>(size1_,vtptr(i,tstp),gcc.tvptr(tstp,ntau_-i));
		  }else{
			  for(i=0;i<=ntau_;i++) element_minusconj<T,LARGESIZE>(size1_,vtptr(i,tstp),gcc.tvptr(tstp,ntau_-i));
		  }
	  }
   }
   delete [] ret_ti;
   delete [] adv_it;
   delete [] les_ti;
   delete [] les_it;
}
template<typename T> void full_matrix<T>::set_timestep(int tstp,herm_pseudo<T> &g,herm_pseudo<T> &gcc){
   int i;
   cplx *ret_ti = new cplx [size1_*size1_];
   cplx *adv_it = new cplx [size1_*size1_];
   cplx *les_ti = new cplx [size1_*size1_];
   cplx *les_it = new cplx [size1_*size1_];
   assert(tstp>=-1 && tstp<=nt_ && tstp<=g.nt() && tstp<=gcc.nt());
   assert(g.size1()==size1_);
   assert(gcc.size1()==size1_);
   assert(g.ntau()==ntau_);
   assert(gcc.ntau()==ntau_);
   if(tstp==-1){
	  memcpy(matptr(0),g.matptr(0),sizeof(cplx)*(ntau_+1)*element_size_);
   }else{
      memcpy(tvptr(tstp,0),g.tvptr(tstp,0),sizeof(cplx)*(ntau_+1)*element_size_);
      if(size1_==1){
		  for(i=0;i<=tstp;i++){
			element_set<T,1>(size1_,ret_ti,g.retptr(tstp,i));
			element_set<T,1>(size1_,les_it,g.lesptr(i,tstp));
			element_minusconj<T,1>(size1_,les_ti,gcc.lesptr(i,tstp));
			element_minusconj<T,1>(size1_,adv_it,gcc.retptr(tstp,i)); // -adv
			element_set<T,1>(size1_,lesptr(tstp,i),les_ti);
			element_set<T,1>(size1_,lesptr(i,tstp),les_it);
			element_set<T,1>(size1_,gtrptr(tstp,i),ret_ti);
			element_set<T,1>(size1_,gtrptr(i,tstp),adv_it);
		  }
		  if(sig_==-1){
			  for(i=0;i<=ntau_;i++) element_conj<T,1>(size1_,vtptr(i,tstp),gcc.tvptr(tstp,ntau_-i));
		  }else{
			  for(i=0;i<=ntau_;i++) element_minusconj<T,1>(size1_,vtptr(i,tstp),gcc.tvptr(tstp,ntau_-i));
		  }
	  }else{
		  for(i=0;i<=tstp;i++){
			element_set<T,LARGESIZE>(size1_,ret_ti,g.retptr(tstp,i));
			element_set<T,LARGESIZE>(size1_,les_it,g.lesptr(i,tstp));
			element_minusconj<T,LARGESIZE>(size1_,les_ti,gcc.lesptr(i,tstp));
			element_minusconj<T,LARGESIZE>(size1_,adv_it,gcc.retptr(tstp,i)); // -adv
			element_set<T,LARGESIZE>(size1_,lesptr(tstp,i),les_ti);
			element_set<T,LARGESIZE>(size1_,lesptr(i,tstp),les_it);
			element_set<T,LARGESIZE>(size1_,gtrptr(tstp,i),ret_ti);
			element_set<T,LARGESIZE>(size1_,gtrptr(i,tstp),adv_it);
		  }
		  if(sig_==-1){
			  for(i=0;i<=ntau_;i++) element_conj<T,LARGESIZE>(size1_,vtptr(i,tstp),gcc.tvptr(tstp,ntau_-i));
		  }else{
			  for(i=0;i<=ntau_;i++) element_minusconj<T,LARGESIZE>(size1_,vtptr(i,tstp),gcc.tvptr(tstp,ntau_-i));
		  }
	  }
   }
   delete [] ret_ti;
   delete [] adv_it;
   delete [] les_ti;
   delete [] les_it;
}
template<typename T> void full_matrix<T>::set_timestep(int tstp,herm_matrix_timestep<T> &g,herm_matrix_timestep<T> &gcc){
   int i;
   cplx *ret_ti = new cplx [size1_*size1_];
   cplx *adv_it = new cplx [size1_*size1_];
   cplx *les_ti = new cplx [size1_*size1_];
   cplx *les_it = new cplx [size1_*size1_];
   assert(tstp>=-1 && tstp<=nt_ && tstp==g.tstp_ && tstp==gcc.tstp_);
   assert(g.size1_==size1_);
   assert(gcc.size1_==size1_);
   assert(g.ntau_==ntau_);
   assert(gcc.ntau_==ntau_);
   if(tstp==-1){
	  memcpy(matptr(0),g.matptr(0),sizeof(cplx)*(ntau_+1)*element_size_);
   }else{
	  memcpy(tvptr(tstp,0),g.tvptr(0),sizeof(cplx)*(ntau_+1)*element_size_);
      if(size1_==1){
			for(i=0;i<=tstp;i++){
			element_set<T,1>(size1_,ret_ti,g.retptr(i));
			element_set<T,1>(size1_,les_it,g.lesptr(i));
			element_minusconj<T,1>(size1_,les_ti,gcc.lesptr(i));
			element_minusconj<T,1>(size1_,adv_it,gcc.retptr(i)); // -adv
			element_set<T,1>(size1_,lesptr(i,tstp),les_it);
			element_set<T,1>(size1_,lesptr(tstp,i),les_ti);
			element_incr<T,1>(size1_,ret_ti,les_ti);
			element_incr<T,1>(size1_,adv_it,les_it);
			element_set<T,1>(size1_,gtrptr(i,tstp),adv_it);
			element_set<T,1>(size1_,gtrptr(tstp,i),ret_ti);
		  }
		  if(sig_==-1){
			  for(i=0;i<=ntau_;i++) element_conj<T,1>(size1_,vtptr(i,tstp),gcc.tvptr(ntau_-i));
		  }else{
			  for(i=0;i<=ntau_;i++) element_minusconj<T,1>(size1_,vtptr(i,tstp),gcc.tvptr(ntau_-i));
		  }
	  }else{
		  for(i=0;i<=tstp;i++){
			element_set<T,LARGESIZE>(size1_,ret_ti,g.retptr(i));
			element_set<T,LARGESIZE>(size1_,les_it,g.lesptr(i));
			element_minusconj<T,LARGESIZE>(size1_,les_ti,gcc.lesptr(i));
			element_minusconj<T,LARGESIZE>(size1_,adv_it,gcc.retptr(i)); // -adv
			element_set<T,LARGESIZE>(size1_,lesptr(i,tstp),les_it);
			element_set<T,LARGESIZE>(size1_,lesptr(tstp,i),les_ti);
			element_incr<T,LARGESIZE>(size1_,ret_ti,les_ti);
			element_incr<T,LARGESIZE>(size1_,adv_it,les_it);
			element_set<T,LARGESIZE>(size1_,gtrptr(i,tstp),adv_it);
			element_set<T,LARGESIZE>(size1_,gtrptr(tstp,i),ret_ti);
		  }
		  if(sig_==-1){
			  for(i=0;i<=ntau_;i++) element_conj<T,LARGESIZE>(size1_,vtptr(i,tstp),gcc.tvptr(ntau_-i));
		  }else{
			  for(i=0;i<=ntau_;i++) element_minusconj<T,LARGESIZE>(size1_,vtptr(i,tstp),gcc.tvptr(ntau_-i));
		  }
	  }
   }
   delete [] ret_ti;
   delete [] adv_it;
   delete [] les_ti;
   delete [] les_it;
}
template<typename T> void full_matrix<T>::get_timestep(int tstp,herm_matrix_timestep<T> &timestep){
   int i;
   cplx *x = new cplx [size1_*size1_];
   assert(tstp>=-1 && tstp<=nt_ && tstp==timestep.tstp_);
   assert(timestep.size1_==size1_);
   assert(timestep.ntau_==ntau_);
   if(tstp==-1){
	  memcpy(timestep.matptr(0),matptr(0),sizeof(cplx)*(ntau_+1)*element_size_);
   }else{
      memcpy(timestep.retptr(0),gtrptr(tstp,0),sizeof(cplx)*(tstp+1)*element_size_);
	  memcpy(timestep.tvptr(0),tvptr(tstp,0),sizeof(cplx)*(ntau_+1)*element_size_);
      if(size1_==1){
		  for(i=0;i<=tstp;i++){
			element_set<T,1>(size1_,x,lesptr(tstp,i));
			element_smul<T,1>(size1_,x,-1.0);
			element_incr<T,1>(size1_,timestep.retptr(i),x);
			element_conj<T,1>(size1_,timestep.lesptr(i),x);
		  }
	  }else{
		  for(i=0;i<=tstp;i++){
			element_set<T,LARGESIZE>(size1_,x,lesptr(tstp,i));
			element_smul<T,LARGESIZE>(size1_,x,-1.0);
			element_incr<T,LARGESIZE>(size1_,timestep.retptr(i),x);
			element_conj<T,LARGESIZE>(size1_,timestep.lesptr(i),x);
		  }
	  }
   }
   delete [] x;
}

} // namespace cntr

#endif  // CNTR_FULL_MATRIX_IMPL
