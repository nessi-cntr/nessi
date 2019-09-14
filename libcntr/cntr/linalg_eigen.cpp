///////////////////////////////////////////////////////////////////////////////////////////////////////////
// interface to simple linear algebra routines:
// implementation using EIGEN library, ans a lot of copying back and forth
///////////////////////////////////////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include "eigen_typedef.h"
#include "cntr_global_settings.hpp"

namespace linalg {

// Consider using Eigen:Map instead of doing copying back and forth

void set_cdmatrix(int n,void *a,cdmatrix &A){
	cdouble *aa=(cdouble*)a;
    A=Eigen::Map<cdmatrix>(aa,n,n).transpose();
}
void set_dmatrix(int n,void *a,dmatrix &A){
	double *aa=(double*)a;
	A=Eigen::Map<dmatrix>(aa,n,n).transpose();
}
void set_cdvector(int n,void *a,cdvector &A){
	cdouble *aa=(cdouble*)a;
    A=Eigen::Map<cdvector>(aa,n);
}

void set_dvector(int n,void *a,dvector &A){
    double *aa=(double*)a;
    A=Eigen::Map<dvector>(aa,n);
}
    
void get_cdvector(int n,void *a,cdvector &A){
	int i;
	cdouble *aa=(cdouble*)a;
	for(i=0;i<n;i++) aa[i]=A(i);
}
void get_dvector(int n,void *a,dvector &A){
	int i;
	double *aa=(double*)a;
	for(i=0;i<n;i++) aa[i]=A(i);
}
void get_cdmatrix(int n,void *a,cdmatrix &A){
	int i,j;
	cdouble *aa=(cdouble*)a;
	for(i=0;i<n;i++) for(j=0;j<n;j++) aa[i*n+j]=A(i,j);
}
void get_dmatrix(int n,void *a, const dmatrix &A){
	int i,j;
	double *aa=(double*)a;
	for(i=0;i<n;i++) for(j=0;j<n;j++) aa[i*n+j]=A(i,j);
}

//old loop style

//void set_cdmatrix(int n,void *a,cdmatrix &A){
//    int i,j;
//    cdouble *aa=(cdouble*)a;
//    A.resize(n,n);
//    for(i=0;i<n;i++) for(j=0;j<n;j++) A(i,j)=aa[i*n+j];
//}
//void set_dmatrix(int n,void *a,dmatrix &A){
//    int i,j;
//    double *aa=(double*)a;
//    A.resize(n,n);
//    for(i=0;i<n;i++) for(j=0;j<n;j++) A(i,j)=aa[i*n+j];
//}
//void set_cdvector(int n,void *a,cdvector &A){
//    int i;
//    cdouble *aa=(cdouble*)a;
//    A.resize(n);	
//    for(i=0;i<n;i++) A(i)=aa[i];
//}
//void set_dvector(int n,void *a,dvector &A){
//	int i;
//	double *aa=(double*)a;
//	A.resize(n);	
//	for(i=0;i<n;i++) A(i)=aa[i];
//}
    
void cplx_sq_solve(void *a,void  *b,void *x,int n)
{
   cdmatrix A_eigen;
   cdvector b_eigen,X_eigen;   
   set_cdmatrix(n,a,A_eigen);
   set_cdvector(n,b,b_eigen);
   Eigen::FullPivLU<cdmatrix> lu(A_eigen);
   X_eigen=lu.solve(b_eigen);
   get_cdvector(n,x,X_eigen);
}
// solve d mXm problems, n=m*d 
void cplx_sq_solve_many(void *a,void  *b,void *x,int n,int d)
{
   int l;
   cdouble *b1=(cdouble*)b;
   cdouble *x1=(cdouble*)x;
   cdmatrix A_eigen;
   cdvector b_eigen,X_eigen;
   set_cdmatrix(n,a,A_eigen);
   Eigen::FullPivLU<cdmatrix> lu(A_eigen);
   for(l=0;l<d;l++){
	   set_cdvector(n,b1+n*l,b_eigen);
	   X_eigen=lu.solve(b_eigen);
	   get_cdvector(n,x1+n*l,X_eigen);
	}
}
void real_sq_solve(double *a,double  *b,double *x,int n)
{
   dmatrix A_eigen;
   dvector b_eigen,X_eigen;   
   set_dmatrix(n,a,A_eigen);
   set_dvector(n,b,b_eigen);
   Eigen::FullPivLU<dmatrix> lu(A_eigen);
   X_eigen=lu.solve(b_eigen);
   get_dvector(n,x,X_eigen);
}
void cplx_matrix_inverse(void *a,void *x,int n)
{
   cdmatrix A_eigen;
   cdmatrix X_eigen;   
   set_cdmatrix(n,a,A_eigen);
   Eigen::FullPivLU<cdmatrix> lu(A_eigen);
   X_eigen=lu.inverse();
   get_cdmatrix(n,x,X_eigen);
}
void linalg_matrix_inverse(double *a,double *x,int n)
{
   dmatrix A_eigen;
   dmatrix X_eigen;   
   set_dmatrix(n,a,A_eigen);
   Eigen::FullPivLU<dmatrix> lu(A_eigen);
   X_eigen=lu.inverse();
   get_dmatrix(n,x,X_eigen);
}
void eigen_hermv(int n,std::complex<double> *A,double *eval,std::complex<double> *evec){
	cdmatrix A_eigen;
	cdmatrix evec_eigen;   
	dvector eval_eigen;
	set_cdmatrix(n,A,A_eigen);
	Eigen::SelfAdjointEigenSolver<cdmatrix> eigensolver(A_eigen);
	eval_eigen=eigensolver.eigenvalues();
	evec_eigen=eigensolver.eigenvectors();
	get_cdmatrix(n,evec,evec_eigen);
	get_dvector(n,eval,eval_eigen);
}

//void QR_decomposition(double *aa,double *qq,double *rr,int n,int m) {
//
//  assert( n == m );
//
//  dmatrix mat_aa;
//  set_dmatrix(n, aa, mat_aa);
//  Eigen::ColPivHouseholderQR<dmatrix> qr(mat_aa);
//
//  dmatrix mat_r = qr.matrixQR().triangularView<Eigen::Upper>();
//  dmatrix mat_q = qr.householderQ();
//
//  get_dmatrix(n, qq, mat_q);
//  get_dmatrix(n, rr, mat_r);
//
//}

} // end namespace linalg
