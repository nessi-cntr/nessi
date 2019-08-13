#ifndef CNTR_ELEMENTS_H
#define CNTR_ELEMENTS_H

#include "eigen_map.hpp"
#include "linalg.hpp"

namespace cntr {

// a matrix-matrix multiplication, like gemm in blas

#define LARGESIZE (-1) // Fall back to dynamic size
#define BLAS_SIZE 4    // use blas for larger matrices
// #define USE_BLAS
/* #######################################################################################
#
#   general matrix operations ...
#
########################################################################################*/
#define CPLX std::complex<T>
#define VALUE_TYPE T
/// @private
template<typename T,int DIM1,int DIM2>
inline void element_set_zero(int size1,int size2,CPLX * z){
  memset(z,0,sizeof(CPLX)*size1*size2);
}
/// @private
template<typename T,int DIM1,int DIM2>
inline void element_set(int size1,int size2,CPLX *z,CPLX *z1) {
  memcpy(z,z1,sizeof(CPLX)*size1*size2);
}
/// @private
template<typename T,int DIM1,int DIM2>
inline void element_swap(int size1,int size2,CPLX *z,CPLX *z1) {
  int nm=size1*size2,i;
  CPLX x;
  for(i=0;i<nm;i++){x=z[i];z[i]=z1[i];z1[i]=x;}
}
// hermitian conjugate
/// @private
template<typename T,int DIM1,int DIM2>
inline void element_conj(int size1,int size2,CPLX * z,CPLX *z1){
     // z1 is of size [ size2 x size1 ]
	 int l,m;
     CPLX temp;
     for(l=0;l<size1;l++){
       for(m=0;m<size2;m++){
         temp=z1[m*size1+l];
		 z[l*size2 + m ] = CPLX(temp.real(),-temp.imag());
       }
     }
}
/// @private
template<typename T,int DIM1,int DIM2>
inline void element_conj(int size1,int size2,CPLX * z){
   CPLX *z1 = new CPLX [size1*size2];
   element_set<T,DIM2,DIM1>(size2,size1,z1,z);
   element_conj<T,DIM1,DIM2>(size1,size2,z,z1);
   delete [] z1;
}
/// @private
template<typename T,int DIM1,int DIM2>
inline void element_minusconj(int size1,int size2,CPLX * z,CPLX *z1){
     // z1 is of size [ size2 x size1 ]
	 int l,m;
     CPLX temp; 
     for(l=0;l<size1;l++){
       for(m=0;m<size2;m++){
         temp=z1[m*size1+l];
		 z[l*size2 + m ] = CPLX(-temp.real(),temp.imag());
       }
     }
}
/// @private
template<typename T,int DIM1,int DIM2>
inline void element_minusconj(int size1,int size2,CPLX * z){
   CPLX *z1 = new CPLX [size1*size2];
   element_set<T,DIM2,DIM1>(size2,size1,z1,z);
   element_minusconj<T,DIM1,DIM2>(size1,size2,z,z1);
   delete [] z1;
}
/// @private
template<typename T,int DIM1,int DIM2>
inline T element_norm2(int size1,int size2,CPLX * z){
  int i,nm=size1*size2;
  T r=0.0;
  for(i=0;i<nm;i++) r += z[i].real()*z[i].real()+z[i].imag()*z[i].imag();
  return sqrt(r);
}
// scalar multiplication
/// @private
template<typename T,int DIM1,int DIM2>
inline void element_smul(int size1,int size2,CPLX *z,T z0) { 
  int nm=size1*size2,i;
  for(i=0;i<nm;i++) z[i] *= z0;
}
/// @private
template<typename T,int DIM1,int DIM2>
inline void element_smul(int size1,int size2,CPLX *z,CPLX z0) { 
  int nm=size1*size2,i;
  for(i=0;i<nm;i++) z[i] *= z0;
}
// multiplication
/// @private
template<typename T,int DIM1,int DIM2,int DIM3>
inline void element_mult(int size1,int size2,int size3,CPLX *z,CPLX *z1,CPLX *z2){
    element_map<DIM1, DIM2>(size1, size2, z) =
        element_map<DIM1, DIM3>(size1, size3, z1) *
        element_map<DIM3, DIM2>(size3, size2, z2);
}
// increment with product and other
/// @private
template<typename T,int DIM1,int DIM2>
inline void element_incr(int size1,int size2,CPLX *z,CPLX *z1) { 
    if (z != z1) {
        element_map<DIM1, DIM2>(size1, size2, z).noalias() +=
            element_map<DIM1, DIM2>(size1, size2, z1);
    } else {
        element_map<DIM1, DIM2>(size1, size2, z) +=
            element_map<DIM1, DIM2>(size1, size2, z1);
    }
}
/// @private
template<typename T,int DIM1,int DIM2>
inline void element_incr(int size1,int size2,CPLX *z,CPLX alpha,CPLX *z1) { 
    if (z != z1) {
        element_map<DIM1, DIM2>(size1, size2, z).noalias() +=
            alpha * element_map<DIM1, DIM2>(size1, size2, z1);
    } else {
        element_map<DIM1, DIM2>(size1, size2, z) +=
            alpha * element_map<DIM1, DIM2>(size1, size2, z1);
    }
}
/// @private
template<typename T,int DIM1,int DIM2,int DIM3>
inline void element_incr(int size1,int size2,int size3,CPLX *z,CPLX *z1,CPLX *z2){
    if (z != z1 && z != z2) {
        element_map<DIM1, DIM2>(size1, size2, z).noalias() +=
            element_map<DIM1, DIM3>(size1, size3, z1) *
            element_map<DIM3, DIM2>(size3, size2, z2);
    } else {
        element_map<DIM1, DIM2>(size1, size2, z) +=
            element_map<DIM1, DIM3>(size1, size3, z1) *
            element_map<DIM3, DIM2>(size3, size2, z2);
    }
}
/// @private
template<typename T,int DIM1,int DIM2,int DIM3>
inline void element_incr(int size1,int size2,int size3,CPLX *z,CPLX alpha,CPLX *z1,CPLX *z2){
    if (z != z1 && z != z2) {
        element_map<DIM1, DIM2>(size1, size2, z).noalias() +=
            alpha * (element_map<DIM1, DIM3>(size1, size3, z1) *
                     element_map<DIM3, DIM2>(size3, size2, z2));
    } else {
        element_map<DIM1, DIM2>(size1, size2, z) +=
            alpha * (element_map<DIM1, DIM3>(size1, size3, z1) *
                     element_map<DIM3, DIM2>(size3, size2, z2));
    }
}
// specialization for DIM=1, T=double
#undef CPLX
#define CPLX std::complex<double>
#undef VALUE_TYPE
#define VALUE_TYPE double
/// @private
template<> inline void element_set_zero<VALUE_TYPE,1,1>(int size1,int size2,CPLX * z){
  *z=0.0;
}
/// @private
template<> inline void element_set<VALUE_TYPE,1,1>(int size1,int size2,CPLX *z,CPLX *z1) {
  *z=*z1;
}
/// @private
template<> inline void element_swap<VALUE_TYPE,1,1>(int size1,int size2,CPLX *z,CPLX *z1) {
    CPLX x;
    x=*z;*z=*z1;*z1=x;
}
/// @private
template<> inline void element_conj<VALUE_TYPE,1,1>(int size1,int size2,CPLX * z,CPLX *z1){
     CPLX temp=*z1;
     //temp.imag()=-temp.imag();
     temp = std::conj( temp );
	 *z=temp;
}
/// @private
template<> inline void element_conj<VALUE_TYPE,1,1>(int size1,int size2,CPLX * z){
   CPLX temp=*z;
   //temp.imag()=-temp.imag();
   temp = std::conj( temp );
   *z=temp;
}
/// @private
template<> inline void element_minusconj<VALUE_TYPE,1,1>(int size1,int size2,CPLX * z,CPLX *z1){
     CPLX temp=*z1;
     //temp.real()=-temp.real();
     temp = -std::conj( temp );
	 *z=temp;
}
/// @private
template<> inline void element_minusconj<VALUE_TYPE,1,1>(int size1,int size2,CPLX * z){
   CPLX temp=*z;
   //temp.real()=-temp.real();
   temp = -std::conj( temp );
   *z=temp;
}
/// @private
template<> inline VALUE_TYPE element_norm2<VALUE_TYPE,1,1>(int size1,int size2,CPLX * z){
  CPLX x=*z;
  return sqrt(x.real()*x.real()+x.imag()*x.imag());
}
/// @private
template<> inline void element_smul<VALUE_TYPE,1,1>(int size1,int size2,CPLX *z,VALUE_TYPE alpha) {
  *z=*z*alpha;
}
/// @private
template<> inline void element_smul<VALUE_TYPE,1,1>(int size1,int size2,CPLX *z,CPLX alpha){
  *z=*z*alpha;
}
/// @private
template<> inline void element_mult<VALUE_TYPE,1,1,1>(int size1,int size2,int size3,CPLX *z,CPLX *z1,CPLX *z2){
  *z=(*z1)*(*z2);
}
/// @private
template<> inline void element_incr<VALUE_TYPE,1,1>(int size1,int size2,CPLX *z,CPLX *z1) {
  *z=*z+*z1;
}
/// @private
template<> inline void element_incr<VALUE_TYPE,1,1>(int size1,int size2,CPLX *z,CPLX alpha,CPLX *z1) {
  *z=*z+alpha*(*z1);
}
/// @private
template<> inline void element_incr<VALUE_TYPE,1,1,1>(int size1,int size2,int size3,CPLX *z,CPLX *z1,CPLX *z2){
  *z=*z+(*z1)*(*z2);
}
/// @private
template<> inline void element_incr<VALUE_TYPE,1,1,1>(int size1,int size2,int size3,CPLX *z,CPLX alpha,CPLX *z1,CPLX *z2){
   *z=*z+alpha*(*z1)*(*z2);
}
// specialization for DIM=1, T=float
#undef CPLX
#define CPLX std::complex<float>
#undef VALUE_TYPE
#define VALUE_TYPE float
/// @private
template<> inline void element_set_zero<VALUE_TYPE,1,1>(int size1,int size2,CPLX * z){
  *z=0.0;
}
/// @private
template<> inline void element_set<VALUE_TYPE,1,1>(int size1,int size2,CPLX *z,CPLX *z1) {
  *z=*z1;
}
/// @private
template<> inline void element_swap<VALUE_TYPE,1,1>(int size1,int size2,CPLX *z,CPLX *z1) {
  CPLX x;
  x=*z;*z=*z1;*z1=x;
}
/// @private
template<> inline void element_conj<VALUE_TYPE,1,1>(int size1,int size2,CPLX * z,CPLX *z1){
     CPLX temp=*z1;
     //temp.imag()=-temp.imag();
     temp = std::conj( temp );
	 *z=temp;
}
/// @private
template<> inline void element_conj<VALUE_TYPE,1,1>(int size1,int size2,CPLX * z){
   CPLX temp=*z;
   //temp.imag()=-temp.imag();
   temp = std::conj( temp );
   *z=temp;
}
/// @private
template<> inline void element_minusconj<VALUE_TYPE,1,1>(int size1,int size2,CPLX * z,CPLX *z1){
     CPLX temp=*z1;
     //temp.real()=-temp.real();
   temp = -std::conj( temp );
	 *z=temp;
}
/// @private
template<> inline void element_minusconj<VALUE_TYPE,1,1>(int size1,int size2,CPLX * z){
   CPLX temp=*z;
   //temp.real()=-temp.real();
   temp = -std::conj( temp );
   *z=temp;
}
/// @private
template<> inline VALUE_TYPE element_norm2<VALUE_TYPE,1,1>(int size1,int size2,CPLX * z){
  CPLX x=*z;
  return sqrt(x.real()*x.real()+x.imag()*x.imag());
}
/// @private
template<> inline void element_smul<VALUE_TYPE,1,1>(int size1,int size2,CPLX *z,VALUE_TYPE alpha) {
  *z=*z*alpha;
}
/// @private
template<> inline void element_smul<VALUE_TYPE,1,1>(int size1,int size2,CPLX *z,CPLX alpha) {
  *z=*z*alpha;
}
/// @private
template<> inline void element_mult<VALUE_TYPE,1,1,1>(int size1,int size2,int size3,CPLX *z,CPLX *z1,CPLX *z2){
  *z=(*z1)*(*z2);
}
/// @private
template<> inline void element_incr<VALUE_TYPE,1,1>(int size1,int size2,CPLX *z,CPLX *z1) {
  *z=*z+*z1;
}
/// @private
template<> inline void element_incr<VALUE_TYPE,1,1>(int size1,int size2,CPLX *z,CPLX alpha,CPLX *z1) {
  *z=*z+alpha*(*z1);
}
/// @private
template<> inline void element_incr<VALUE_TYPE,1,1,1>(int size1,int size2,int size3,CPLX *z,CPLX *z1,CPLX *z2){
  *z=*z+(*z1)*(*z2);
}
/// @private
template<> inline void element_incr<VALUE_TYPE,1,1,1>(int size1,int size2,int size3,CPLX *z,CPLX alpha,CPLX *z1,CPLX *z2){
   *z=*z+alpha*(*z1)*(*z2);
}

/* #######################################################################################
#
#   square matrix operations ... partly derived from the general ones,
#                                from which they inherit the specialization
#
########################################################################################*/
#undef CPLX
#define CPLX std::complex<T>
#undef VALUE_TYPE
#define VALUE_TYPE T
/// @private
template<typename T,int DIM>
inline void element_set_zero(int size1,CPLX * z){ 
   element_set_zero<T,DIM,DIM>(size1,size1,z);
}
/// @private
template<typename T,int DIM>
inline void element_set(int size1,CPLX *z,CPLX z0) { 
  element_set_zero<T,DIM>(size1,z);
  for(int i=0;i<size1;i++) z[i*size1+i]=z0;
}
/// @private
template<typename T,int DIM>
inline void element_set(int size1,CPLX *z,CPLX *z1) { 
  element_set<T,DIM,DIM>(size1,size1,z,z1);
}
/// @private
template<typename T,int DIM>
inline void element_swap(int size1,CPLX *z,CPLX *z1) { 
  element_swap<T,DIM,DIM>(size1,size1,z,z1);
}
/// @private
template<typename T,int DIM>
inline CPLX element_trace(int size1,CPLX * z){ 
  CPLX result=0.0;
  for(int i=0;i<size1;i++) result += z[i*size1+i];
  return result;
}
// hermitian conjugate
/// @private
template<typename T,int DIM>
inline void element_conj(int size1,CPLX * z,CPLX *z1){
     element_conj<T,DIM,DIM>(size1,size1,z,z1);
}
/// @private
template<typename T,int DIM>
inline void element_conj(int size1,CPLX * z){
   element_conj<T,DIM,DIM>(size1,size1,z);
}
/// @private
template<typename T,int DIM>
inline void element_minusconj(int size1,CPLX * z,CPLX *z1){
     element_minusconj<T,DIM,DIM>(size1,size1,z,z1);
}
/// @private
template<typename T,int DIM>
inline void element_minusconj(int size1,CPLX * z){
   element_minusconj<T,DIM,DIM>(size1,size1,z);
}
/// @private
template<typename T,int DIM>
inline T element_norm2(int size1,CPLX * z){
  return element_norm2<T,DIM,DIM>(size1,size1,z);
}
// scalar multiplication
/// @private
template<typename T,int DIM>
inline void element_smul(int size1,CPLX *z,T z0) { 
  element_smul<T,DIM,DIM>(size1,size1,z,z0);
}
/// @private
template<typename T,int DIM>
inline void element_smul(int size1,CPLX *z,CPLX z0) { 
  element_smul<T,DIM,DIM>(size1,size1,z,z0);
}
// multiplication
/// @private
template<typename T,int DIM>
inline void element_mult(int size1,CPLX *z,CPLX *z1,CPLX *z2){ 
   element_mult<T,DIM,DIM,DIM>(size1,size1,size1,z,z1,z2);
}
// increment with product etc.
/// @private
template<typename T,int DIM>
inline void element_incr(int size1,CPLX *z,CPLX z0) {
  for(int i=0;i<size1;i++) z[i*size1+i]+=z0;
}
/// @private
template<typename T,int DIM>
inline void element_incr(int size1,CPLX *z,CPLX *z1) {
    element_incr<T,DIM,DIM>(size1,size1,z,z1);
}
/// @private
template<typename T,int DIM>
inline void element_incr(int size1,CPLX *z,CPLX alpha,CPLX *z1) {
    element_incr<T,DIM,DIM>(size1,size1,z,alpha,z1);
}
/// @private
template<typename T,int DIM>
inline void element_incr(int size1,CPLX *z,CPLX *z1,CPLX *z2){
    element_incr<T,DIM,DIM,DIM>(size1,size1,size1,z,z1,z2);
}
/// @private
template<typename T,int DIM>
inline void element_incr(int size1,CPLX *z,CPLX alpha,CPLX *z1,CPLX *z2){
    element_incr<T,DIM,DIM,DIM>(size1,size1,size1,z,alpha,z1,z2);
}
// matrix inverse , use double precision because i am lazy to look
// for the float inverse routine
/// @private
template<typename T,int DIM>
inline void element_inverse(int size1,CPLX *z,CPLX *z1) {
   int i,nm=size1*size1;
   std::complex<double> *zd = new std::complex<double> [size1*size1];
   std::complex<double> *z1d = new std::complex<double> [size1*size1];
   for(i=0;i<nm;i++) z1d[i]= (std::complex<double>) z1[i];
   linalg::cplx_matrix_inverse(z1d,zd,size1);
   for(i=0;i<nm;i++) z[i]= (CPLX)zd[i];
   delete [] zd;
   delete [] z1d;
}
// solution of linear equation  M_ij X_j =  Q_i (right) and X_j M_ij =Q_i (left)
// make this more efficient at a later time, allowing for a reuse of the matrix M (?)
/// @private
template<typename T,int DIM>
inline void element_linsolve_right(int size1,int n,CPLX *X,CPLX *M,CPLX *Q){
  // this uses currently always double precision
  int size=size1,llen=size*n,l,s2=size*size;
  int p,q,m;
  std::complex<double> *mtemp = new std::complex<double> [llen*llen];
  std::complex<double> *qtemp = new std::complex<double> [n*s2];
  std::complex<double> *xtemp = new std::complex<double> [n*s2];
  // set up the matrix - reshuffle elements line by line
  for(l=0;l<n;l++){
	for(m=0;m<n;m++){
	  for(p=0;p<size;p++){
		for(q=0;q<size;q++){
		  mtemp[ (l*size + p)*llen + m*size+q] = (std::complex<double>) M[ (l*n+m)*s2+p*size+q];
		}
	  }
	}
  }
  for(l=0;l<n;l++){
	for(p=0;p<size;p++){
	   for(q=0;q<size;q++){
		 qtemp[ q*n*size + l*size+p] = (std::complex<double>) Q[ l*s2+p*size+q];
	  }
	}
  }
  linalg::cplx_sq_solve_many(mtemp,qtemp,xtemp,llen,size);
  for(l=0;l<n;l++){
	for(p=0;p<size;p++){
	   for(q=0;q<size;q++){
		 X[ l*s2+p*size+q]= (CPLX) xtemp[ q*n*size + l*size+p];
	  }
	}
  }
  delete [] xtemp;
  delete [] qtemp;
  delete [] mtemp;
}
/// @private
template<typename T,int DIM>
inline void element_linsolve_right(int size1,CPLX *X,CPLX *M,CPLX *Q){
   element_linsolve_right<T,DIM>(size1,1,X,M,Q);
}
/// @private
template<typename T,int DIM>
inline void element_linsolve_left(int size1,int n,CPLX *X,CPLX *M,CPLX *Q){
 int l,m,element_size=size1*size1;
 // transpose all elements
 for(l=0;l<n;l++) for(m=0;m<n;m++) element_conj<T,DIM>(size1, M + (l*n+m)*element_size);
 for(l=0;l<n;l++)  element_conj<T,DIM>(size1, Q + l*element_size);
 element_linsolve_right<T,DIM>(size1,n,X,M,Q);
 for(l=0;l<n;l++)  element_conj<T,DIM>(size1, X + l*element_size);
}
/// @private
template<typename T,int DIM>
inline void element_linsolve_left(int size1,CPLX *X,CPLX *M,CPLX *Q){
   element_linsolve_left<T,DIM>(size1,1,X,M,Q);
}

// specialization of square matrix funtions which are not derived from general
// for DIM=1,T=double
#undef CPLX
#define CPLX std::complex<double>
#undef VALUE_TYPE
#define VALUE_TYPE double
/// @private
template<> inline void element_set<VALUE_TYPE,1>(int size1,CPLX *z,CPLX z0) {
  *z=z0;
}
/// @private
template<> inline CPLX element_trace<VALUE_TYPE,1>(int size1,CPLX * z){
  return *z;
}
/// @private
template<> inline void element_incr<VALUE_TYPE,1>(int size1,CPLX *z,CPLX z0) {
  *z=*z+z0;
}
/// @private
template<> inline void element_linsolve_right<VALUE_TYPE,1>(int size1,CPLX *X,CPLX *M,CPLX *Q) {
   *X = *Q / (*M);
}
/// @private
template<> inline void element_linsolve_left<VALUE_TYPE,1>(int size1,CPLX *X,CPLX *M,CPLX *Q) {
   *X = *Q / (*M);
}
/// @private
template<> inline void element_inverse<VALUE_TYPE,1>(int size1,CPLX *z,CPLX *z1) {
  *z= 1.0/(*z1);
}
// specialization of square matrix funtions which are not derived from general
// for DIM=1,T=float
#undef CPLX
#define CPLX std::complex<float>
#undef VALUE_TYPE
#define VALUE_TYPE float
/// @private
template<> inline void element_set<VALUE_TYPE,1>(int size1,CPLX *z,CPLX z0) {
  *z=z0;
}
/// @private
template<> inline CPLX element_trace<VALUE_TYPE,1>(int size1,CPLX * z){
  return *z;
}
/// @private
template<> inline void element_incr<VALUE_TYPE,1>(int size1,CPLX *z,CPLX z0) {
  *z=*z+z0;
}
/// @private
template<> inline void element_linsolve_right<VALUE_TYPE,1>(int size1,CPLX *X,CPLX *M,CPLX *Q) {
   *X = *Q / (*M);
}
/// @private
template<> inline void element_linsolve_left<VALUE_TYPE,1>(int size1,CPLX *X,CPLX *M,CPLX *Q) {
   *X = *Q / (*M);
}
/// @private
template<> inline void element_inverse<VALUE_TYPE,1>(int size1,CPLX *z,CPLX *z1) {
  *z= CPLX(1.0,0.0)/(*z1);
}


#undef CPLX
#undef VALUE_TYPE

}  // namespace cntr

#endif  // CNTR_ELEMENTS_H
