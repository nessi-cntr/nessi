#include "catch.hpp"
#include "cntr.hpp"
#include "vector"
#include <iostream>

#define CPLX_VEC vector<std::complex<double> >
#define DOUBLE_VEC vector<double>
#define CPLX std::complex<double>

using namespace linalg;
using namespace std;

//TO DO:
//

void set_cdvector(int n,void *a,cdvector &A){
	int i;
	cdouble *aa=(cdouble*)a;
	A.resize(n);
	for(i=0;i<n;i++) A(i)=aa[i];
}
void set_dvector(int n,void *a,dvector &A){
	int i;
	double *aa=(double*)a;
	A.resize(n);
	for(i=0;i<n;i++) A(i)=aa[i];
}

void set_cdmatrix(int n,void *a,cdmatrix &A){
	int i,j;
	cdouble *aa=(cdouble*)a;
	A.resize(n,n);
	for(i=0;i<n;i++) for(j=0;j<n;j++) A(i,j)=aa[i*n+j];
}
void set_dmatrix(int n,void *a,dmatrix &A){
	int i,j;
	double *aa=(double*)a;
	A.resize(n,n);
	for(i=0;i<n;i++) for(j=0;j<n;j++) A(i,j)=aa[i*n+j];
}

TEST_CASE("Linear Algebra","[Linear Algebra]"){
    
    double eps=1e-7;
    
    //CPLX_VEC cA,cb,cx;
    // DOUBLE_VEC dA,db,dx;
    
    
    int dim=3;
    int Msize=dim*dim;
    int d=2;
    complex<double> I(0.0,1.0);
    
    SECTION ("real_sq_solve"){
        
        dvector dx_linalg,dx_exact,dx_diff;
        
        double *dA;
        dA = new double[Msize];
        dA[0]=1.0; dA[1]=2.2; dA[2]=3.4;
        dA[3]=4.6; dA[4]=5.7; dA[5]=6.8;
        dA[6]=7.5; dA[7]=8.5; dA[8]=9.6;
        
        double *db;
        db = new double[dim];
        double *dx;
        dx = new double[dim];
        
        db[0]=1.1;db[1]=2.4;db[2]=3.5;
        dx_exact.resize(dim);
        dx_exact<<0.4208144796380182,
                 -0.687782805429876,
                 0.6447963800905043;
        
        //DO we need to specify the size of x before hand?
        real_sq_solve(dA,db,dx,dim);
        
        set_dvector(dim,dx,dx_linalg);
        
        dx_diff=dx_linalg-dx_exact;
        
        REQUIRE(dx_diff.norm()<eps);
        
        delete dA;
        delete db;
        delete dx;
        
    }
    
    
    SECTION ("complex_sq_solve"){
      cdvector cx_linalg,cx_exact,cx_diff;
        CPLX *cA;
        cA = new CPLX[Msize];
        CPLX I(0.0,1.0);
        cA[0]=1.0+0.2*I; cA[1]=2.2; cA[2]=3.4;
        cA[3]=4.6; cA[4]=5.7+1.0*I; cA[5]=6.8;
        cA[6]=7.5; cA[7]=8.5; cA[8]=9.6;
     
        CPLX *cb;
        cb = new CPLX[dim];
        CPLX *cx;
        cx = new CPLX[dim];
     
        cb[0]=1.1+0.1*I;cb[1]=2.4;cb[2]=3.5+0.3*I;
        cx_exact.resize(dim);
        cx_exact<<0.17752793262472366+0.01775980726928994*I,
        -0.1949605345360082+0.015227852835187461*I,
        0.3985109425906921+0.0038921558730451444*I;
     
        cplx_sq_solve(cA,cb,cx,dim);
     
        set_cdvector(dim,cx,cx_linalg);
     
        cx_diff=cx_linalg-cx_exact;
     
        REQUIRE(cx_diff.norm()<eps);
     
        delete cA;
        delete cb;
        delete cx;
     
    }
    
    SECTION ("complex_sq_solve_many"){
        
        cdvector cx_linalg,cx_exact,cx_diff;
        CPLX *cA;
        cA = new CPLX[Msize];
        CPLX I(0.0,1.0);
        cA[0]=1.0+0.2*I; cA[1]=2.2; cA[2]=3.4;
        cA[3]=4.6; cA[4]=5.7+1.0*I; cA[5]=6.8;
        cA[6]=7.5; cA[7]=8.5; cA[8]=9.6;
        
        CPLX *cb;
        cb = new CPLX[dim*d];
        CPLX *cx;
        cx = new CPLX[dim*d];
        
        cb[0]=1.1+0.1*I;cb[1]=2.4;cb[2]=3.5+0.3*I;
        cb[3]=0.3+1.1*I;cb[4]=2.4+1.2*I;cb[5]=1.2+0.7*I;
        
        cx_exact.resize(dim*d);
        cx_exact<<0.17752793262472366+0.01775980726928994*I,
        -0.1949605345360082+0.015227852835187461*I,
        0.3985109425906921+0.0038921558730451444*I,
        -0.13432028547481906+0.23896436702037188*I,
        0.369779042978215-1.5729342395547015*I,
        -0.0974708046097587+1.2789279462043928*I;
        
        cplx_sq_solve_many(cA,cb,cx,dim,d);
        
        set_cdvector(dim*d,cx,cx_linalg);
        
        cx_diff=cx_linalg-cx_exact;
        
        REQUIRE(cx_diff.norm()<eps);
        
        delete cA;
        delete cb;
        delete cx;
     
    }
    
    SECTION ("real_matrix_inverse"){
        
        dmatrix dx_linalg,dx_exact,dx_diff;
        
        double *dA;
        dA = new double[Msize];
        dA[0]=1.0; dA[1]=2.2; dA[2]=3.4;
        dA[3]=4.6; dA[4]=5.7; dA[5]=6.8;
        dA[6]=7.5; dA[7]=8.5; dA[8]=9.6;
        
        double *dx;
        dx = new double[Msize];
        
        dx_exact.resize(dim,dim);
        
        dx_exact<<6.968325791855188,-17.60180995475108,9.99999999999997,
                  -15.47511312217192,35.97285067873294,-19.999999999999947,
                  8.257918552036188, -18.099547511312174, 9.999999999999975;
        
        linalg_matrix_inverse(dA,dx,dim);
        
        set_dmatrix(dim,dx,dx_linalg);
        
        dx_diff=dx_linalg-dx_exact;
        
        REQUIRE(dx_diff.norm()<eps);
        
        delete dA;
        delete dx;
        
    }

    SECTION ("complex_matrix_inverse"){
        
        CPLX I(0.0,1.0);
        cdmatrix cx_linalg,cx_exact,cx_diff;
        CPLX *cA;
        cA = new CPLX[Msize];
        cA[0]=1.0+0.2*I; cA[1]=2.2; cA[2]=3.4;
        cA[3]=4.6; cA[4]=5.7+1.0*I; cA[5]=6.8;
        cA[6]=7.5; cA[7]=8.5; cA[8]=9.6;
     
        CPLX *cx;
        cx = new CPLX[Msize];
     
        cx_exact.resize(dim,dim);
        cx_exact<<-0.5434692749211176-0.2642089142264278*I,-0.06601716570752769+0.46161706554848714*I,0.2392408605773946-0.2334047643083186*I,
        -0.058040798642607985+0.40584328127913327*I,0.020998479292621497-0.9596997815397691*I,0.005682193353650102+0.5360511831376434*I,
        0.4759756614969322-0.15292719105983593*I,0.03298350716866412+0.4890958491119149*I,-0.08777136435796724-0.2922811796205812*I;
     
        cplx_matrix_inverse(cA,cx,dim);
     
        set_cdmatrix(dim,cx,cx_linalg);
     
        cx_diff=cx_linalg-cx_exact;
     
        REQUIRE(cx_diff.norm()<eps);
     
        delete cA;
        delete cx;
     
    }
    SECTION ("eigen_hermv"){
        
        CPLX I(0.0,1.0);
        dvector dval_linalg,dval_exact,dval_diff;
        cdmatrix cvec_linalg,cvec_exact,cvec_diff;
        CPLX *cA;
        cA = new CPLX[Msize];
        cA[0]=1.0; cA[1]=2.2+0.3*I; cA[2]=3.4;
        cA[3]=2.2-0.3*I; cA[4]=5.7; cA[5]=2.7+1.2*I;
        cA[6]=3.4; cA[7]=2.7-1.2*I; cA[8]=9.6;
        
        double *edval;
        CPLX *ecvec;
        ecvec = new CPLX[Msize];
        edval= new double[dim];
        
        dval_exact.resize(dim);
        
        dval_exact<< -0.5428967932641324,4.301152480220869,12.54174431304328;
        cvec_exact.resize(dim,dim);
        cvec_exact<<0.9280589674229315,0.1816887925154496,0.3251088058791392,
        -0.2232869827028977+0.12109020676634462*I,0.3415956672114408-0.7817384646267002*I,0.4464947724487101+0.0912133565995379*I,
        -0.2659825825942915-0.05865069413973232*I,-0.11360282329778541+0.47568997706332894*I,0.822764182049572-0.09841700478046957*I;
        
        eigen_hermv(dim,cA,edval,ecvec);
        
        set_cdmatrix(dim,ecvec,cvec_linalg);
        set_dvector(dim,edval,dval_linalg);
        
        cvec_diff=cvec_linalg-cvec_exact;
        dval_diff=dval_linalg-dval_exact;
        
        REQUIRE(dval_diff.norm()<eps);
        REQUIRE(cvec_diff.norm()<eps);
        
        delete edval;
        delete ecvec;
        delete cA;
        
    }
    
}
