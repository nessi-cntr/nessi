#ifndef NESS2_MATRIX_IMPL_HPP
#define NESS2_MATRIX_IMPL_HPP

#include <cstring>
#include <stdexcept>
#include <complex>
#include <cstdlib>
#include "ness2_matrix_decl.hpp"


namespace ness2 {

  namespace matrix{
    
    //////////////////////////////////////////////////////////////////////////////
    // some matrix helper functions, may be replaced by proper linarg calls
    void setZero(cplx *A,int a1,int a2){
      int dim=a1*a2;
      for (int j=0;j<dim;j++) A[j]=0.0;
    }
    bool isZero(cplx *M,int n,int m){
      int dim=n*m;
      bool ok=true;
      for (int j=0;j<dim;j++) if(M[j]!=0.0) ok=false;
      return ok;      
    }
    
    void setId(cplx *A,int a1){
      int dim=a1*a1;
      for (int j=0;j<dim;j++) A[j]=0.0;
      for (int j=0;j<a1;j++) A[j*a1+j]=1.0;
    }
    void setDiag(cplx *A,int a1,cplx x){
      int dim=a1*a1;
      for (int j=0;j<dim;j++) A[j]=0.0;
      for (int j=0;j<a1;j++) A[j*a1+j]=x;
    }
    cplx trace(cplx *A,int a1){
      cplx res=0.0;
      for (int j=0;j<a1;j++) res+=A[j*a1+j];
      return res;
    }
    void setRandom(cplx *A,int a1,int a2){
      int dim=a1*a2;
      for (int j=0;j<dim;j++){
	double x=((double) rand())/RAND_MAX-0.5;
	double y=((double) rand())/RAND_MAX-0.5;
	A[j]=cplx(x,y);
      }
    }
    void set_minus_hc(cplx *A,int a1,cplx *B,int b1){
      for(int i=0;i<a1;i++){
	for(int j=0;j<a1;j++){
	  A[i*a1+j]=-std::conj(B[j*a1+i]);
	}
      }
    }
    void set_hc(cplx *A,int a1,cplx *B,int b1){
      for(int i=0;i<a1;i++){
	for(int j=0;j<a1;j++){
	  A[i*a1+j]=std::conj(B[j*a1+i]);
	}
      }
    }
    
    void print_matrix(cplx *A,int a1,int a2,int prec){
      std::cout.precision(prec);
      for(int j1=0;j1<a1;j1++){
	for(int j2=0;j2<a2;j2++) std::cout << A[j1*a2+j2] << "\t";
	std::cout << std::endl;
      }
    }
    void incr(cplx *A,int a1,int a2,cplx *B,int b1,int b2,cplx x){
      int dim=a1*a2;
      for (int j=0;j<dim;j++) A[j]+=x*B[j];
    }
    void set(cplx *A,int a1,int a2,cplx *B,int b1,int b2){
      int dim=a1*a2;
      for (int j=0;j<dim;j++) A[j]=B[j];
    }
    void from_eigen(cplx *A,int a1,int a2,cdmatrix &M){
      // sizes must match
      for(int j1=0;j1<a1;j1++){
	for(int j2=0;j2<a2;j2++){
	  A[j1*a2+j2]=M(j1,j2);
	}
      }
    }
    cdmatrix to_eigen(cplx *A,int a1,int a2){
      cdmatrix M(a1,a2);
      for(int j1=0;j1<a1;j1++){
	for(int j2=0;j2<a2;j2++){
	  M(j1,j2)=A[j1*a2+j2];
	}
      }
      return M;
    }


    
    double diff_norm2(cplx *A,int a1,int a2,cplx *B,int b1,int b2){
      int dim=a1*a2;
      double s1=0;
      for (int j=0;j<dim;j++){
	cplx tmp=A[j]-B[j];
	s1+=tmp.real()*tmp.real()+tmp.imag()*tmp.imag();
      }
      return sqrt(s1);
    }
    double norm2(cplx *A,int a1,int a2){
      int dim=a1*a2;
      double s1=0;
      for (int j=0;j<dim;j++){
	cplx tmp=A[j];
	s1+=tmp.real()*tmp.real()+tmp.imag()*tmp.imag();
      }
      return sqrt(s1);
    }
    void smul(cplx *A,int a1,int a2,cplx x){
      int dim=a1*a2;
      for (int j=0;j<dim;j++) A[j] *= x;
    }
    void matrix_vector(cplx *C,cplx *A,int a1,int a2,cplx *B){
      for(int i=0;i<a1;i++){
	cplx tmp=0.0;
	for(int j=0;j<a2;j++){
	  tmp += A[i*a2+j]*B[j];
	}
	C[i] = tmp;
      }
    }
    
    void incr_matrix_mult(cplx *C,int c1,int c2,cplx *A,int a1,int a2,cplx *B,int b1,int b2){
      for(int i=0;i<c1;i++){
	for(int j=0;j<c2;j++){
	  cplx tmp=0.0;
	  for(int l=0;l<a2;l++) tmp += A[i*a2+l]*B[l*b2+j];
	  C[i*c2+j] += tmp;
	}
      }
    }
    void incr_matrix_mult(cplx *C,int c1,int c2,cplx *A,int a1,int a2,cplx *B,int b1,int b2,cplx x1,cplx x2){
      // C = x1*C + x2*(A*B)
      for(int i=0;i<c1;i++){
	for(int j=0;j<c2;j++){
	  cplx tmp=0.0;
	  for(int l=0;l<a2;l++) tmp += A[i*a2+l]*B[l*b2+j];
	  C[i*c2+j] = x1*C[i*c2+j]+x2*tmp;
	}
      }
    }
    void matrix_mult(cplx *C,int c1,int c2,cplx *A,int a1,int a2,cplx *B,int b1,int b2){
      for(int i=0;i<c1;i++){
	for(int j=0;j<c2;j++){
	  cplx tmp=0.0;
	  for(int l=0;l<a2;l++) tmp += A[i*a2+l]*B[l*b2+j];
	  C[i*c2+j] = tmp;
	}
      }
    }
    
    void tensor_prod(cplx *C,int c1,int c2,cplx *A,int a1,int a2,cplx *B,int b1,int b2){
      for(int ia=0;ia<a1;ia++){
	for(int ib=0;ib<b1;ib++){
	  int ic=ia*b1+ib;
	  for(int ja=0;ja<a2;ja++){
	    for(int jb=0;jb<b2;jb++){
	      int jc=ja*b2+jb;
	      C[ic*c2+jc]=A[ia*a2+ja]*B[ib*b2+jb];
	    }
	  }
	}
      }
    }

    ///////////////////////////////////////////////////////////////
    // interface to eigen: copied from nessi, to make a standalone code:
    //////////////////////////////////////////////////////////////
    // Consider using Eigen:Map instead of doing copying back and forth

    void set_cdmatrix(int n,void *a,cdmatrix &A){
      cplx *aa=(cplx*)a;
      A=Eigen::Map<cdmatrix>(aa,n,n).transpose();
    }
    void set_dmatrix(int n,void *a,dmatrix &A){
      double *aa=(double*)a;
      A=Eigen::Map<dmatrix>(aa,n,n).transpose();
    }
    void set_cdvector(int n,void *a,cdvector &A){
      cplx *aa=(cplx*)a;
      A=Eigen::Map<cdvector>(aa,n);
    }
    
    void set_dvector(int n,void *a,dvector &A){
      double *aa=(double*)a;
      A=Eigen::Map<dvector>(aa,n);
    }
    
    void get_cdvector(int n,void *a,cdvector &A){
      int i;
      cplx *aa=(cplx*)a;
      for(i=0;i<n;i++) aa[i]=A(i);
    }
    void get_dvector(int n,void *a,dvector &A){
      int i;
      double *aa=(double*)a;
      for(i=0;i<n;i++) aa[i]=A(i);
    }
    void get_cdmatrix(int n,void *a,cdmatrix &A){
      int i,j;
      cplx *aa=(cplx*)a;
      for(i=0;i<n;i++) for(j=0;j<n;j++) aa[i*n+j]=A(i,j);
    }
    void get_dmatrix(int n,void *a, const dmatrix &A){
      int i,j;
      double *aa=(double*)a;
      for(i=0;i<n;i++) for(j=0;j<n;j++) aa[i*n+j]=A(i,j);
    }
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
      cplx *b1=(cplx*)b;
      cplx *x1=(cplx*)x;
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
    
    // with entries of size*size
    void linsolve_right(int size1,int n,cplx *X,cplx *M,cplx *Q){
      int size=size1,llen=size*n,l,s2=size*size;
      int p,q,m;
      cplx *mtemp = new cplx [llen*llen];
      cplx *qtemp = new cplx [n*s2];
      cplx *xtemp = new cplx [n*s2];
      // set up the matrix - reshuffle elements line by line
      for(l=0;l<n;l++){
	for(m=0;m<n;m++){
	  for(p=0;p<size;p++){
	    for(q=0;q<size;q++){
	      mtemp[ (l*size + p)*llen + m*size+q] = M[ (l*n+m)*s2+p*size+q];
	    }
	  }
	}
      }
      for(l=0;l<n;l++){
	for(p=0;p<size;p++){
	  for(q=0;q<size;q++){
	    qtemp[ q*n*size + l*size+p] = Q[ l*s2+p*size+q];
	  }
	}
      }
      cplx_sq_solve_many(mtemp,qtemp,xtemp,llen,size);
      for(l=0;l<n;l++){
	for(p=0;p<size;p++){
	  for(q=0;q<size;q++){
	    X[ l*s2+p*size+q]= xtemp[ q*n*size + l*size+p];
	  }
	}
      }
      delete [] xtemp;
      delete [] qtemp;
      delete [] mtemp;
    }
       
  } // namespace
  
} // namespace ness2



#endif



