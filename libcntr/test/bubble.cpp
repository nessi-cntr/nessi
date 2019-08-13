#include "catch.hpp"
#include "cntr.hpp"
#include "vector"
#include <iostream>

#define GREEN cntr::herm_matrix<double> 
#define CPLX std::complex<double>

using namespace std;

//TO DO:
// 3. make test for multi orbital case
///////////////Basic Functions////////////////
const complex<double> xj(0.0,1.0);

double f_fermi(double w0,double beta){
	
	return 1.0/(exp(beta*w0)+1.0);
}

double f_bose(double w0,double beta){
	
	return 1.0/(exp(beta*w0)-1.0);
}

///////////////Free Equilibrium Green's functions of 1x1  ////////////////
double G0_matsu(double w0, double beta, double tau,int sign){
	
    assert (abs(sign)==1);
    
    double dvalue;
    
    if(sign==-1){
        if(tau>=0){
            dvalue=-f_fermi(-w0,beta)*exp(-tau*w0);
        }else{
            dvalue=f_fermi(-w0,beta)*exp(-(beta+tau)*w0);
        }
    }else if (sign==1){
        if(tau>=0){
            dvalue=f_bose(-w0,beta)*exp(-tau*w0);
        }else{
            
            dvalue=f_bose(-w0,beta)*exp(-(beta+tau)*w0);
        }
    }
    
    return dvalue;
	
}

CPLX G0_les(double w0, double beta, double t1, double t2, int sign){
    
    assert (abs(sign)==1);
    
    CPLX cvalue;
	
    if(sign==-1){
        cvalue=xj*f_fermi(w0,beta)*exp(-xj*(t1-t2)*w0);
    }else if(sign==1){
        cvalue=-xj*f_bose(w0,beta)*exp(-xj*(t1-t2)*w0);
    }
    
    return cvalue;
	
}

CPLX G0_gtr(double w0, double beta, double t1, double t2, int sign){
    
	assert (abs(sign)==1);
    
    CPLX cvalue;
    if(sign==-1){
        cvalue=-xj*f_fermi(-w0,beta)*exp(-xj*(t1-t2)*w0);
    }else if(sign==1){
        cvalue=xj*f_bose(-w0,beta)*exp(-xj*(t1-t2)*w0);
    }
    
    return cvalue;
}

CPLX G0_tv(double w0, double beta, double t1,double tau2, int sign){
	assert (abs(sign)==1);
    
    CPLX cvalue;
    if(sign==-1){
        cvalue=xj*f_fermi(w0,beta)*exp(-xj*t1*w0)*exp(tau2*w0);
    }else if(sign==1){
        cvalue=-xj*f_bose(w0,beta)*exp(-xj*t1*w0)*exp(tau2*w0);
    }
    
    
    return cvalue;
}

CPLX G0_vt(double w0, double beta, double tau1,double t2, int sign){
	
    assert (abs(sign)==1);
    CPLX cvalue;
    
    if(sign==-1){
        cvalue=-xj*f_fermi(-w0,beta)*exp(-tau1*w0)*exp(xj*t2*w0);
    }else if(sign==1){
        cvalue=xj*f_bose(-w0,beta)*exp(-tau1*w0)*exp(xj*t2*w0);
    }
    
    return cvalue;
}


///////////////Free Equilibrium Green's functions of NxN  ////////////////
double G0_matsu(cdmatrix U, dvector D, double beta, double tau,int sign, int l1,int l2){
    assert(D.size()==U.rows());
    
    int size=D.size();
    cdmatrix tmp(size,size);
    tmp.setZero(size,size);
    
    for(int l=0;l<size;l++) {tmp(l,l)=G0_matsu(D(l),beta,tau,sign);}
    tmp=U.adjoint()*tmp*U;
    return real(tmp(l1,l2));
    
}

CPLX G0_les(cdmatrix U, dvector D, double beta, double t1, double t2, int sign, int l1,int l2){
    assert(D.size()==U.rows());
    
    int size=D.size();
    cdmatrix tmp(size,size);
    tmp.setZero(size,size);
    
    for(int l=0;l<size;l++) {tmp(l,l)=G0_les(D(l),beta,t1,t2,sign);}
    tmp=U.adjoint()*tmp*U;
    
    return tmp(l1,l2);
    
}

CPLX G0_gtr(cdmatrix U, dvector D, double beta, double t1, double t2, int sign, int l1,int l2){
    assert(D.size()==U.rows());
    
    int size=D.size();
    cdmatrix tmp(size,size);
    tmp.setZero(size,size);
    
    for(int l=0;l<size;l++) {tmp(l,l)=G0_gtr(D(l),beta,t1,t2,sign);}
    tmp=U.adjoint()*tmp*U;
    
    return tmp(l1,l2);
    
}

CPLX G0_tv(cdmatrix U, dvector D, double beta, double t1,double tau2, int sign, int l1,int l2){
    
    assert(D.size()==U.rows());
    
    int size=D.size();
    cdmatrix tmp(size,size);
    tmp.setZero(size,size);
    
    for(int l=0;l<size;l++) {tmp(l,l)=G0_tv(D(l),beta,t1,tau2,sign);}
    tmp=U.adjoint()*tmp*U;
    
    return tmp(l1,l2);
    
}

CPLX G0_vt(cdmatrix U, dvector D,double beta, double tau1,double t2, int sign, int l1,int l2){
    assert(D.size()==U.rows());
    
    int size=D.size();
    cdmatrix tmp(size,size);
    tmp.setZero(size,size);
    
    for(int l=0;l<size;l++) {tmp(l,l)=G0_vt(D(l),beta,tau1,t2,sign);}
    tmp=U.adjoint()*tmp*U;
    
    return tmp(l1,l2);
    
}

///////////////Set free equilibrium Green's functions of 1x1 ////////////////
void set_G0_ref_size1(GREEN & A_ref, double beta, double h, double wa, int sign_a){
    
    int nt_=A_ref.nt();
    int ntau_=A_ref.ntau();
    int tstp=-1;
    double dtau=beta/(double)ntau_;
    
    cdmatrix tmp(1,1);
    //Matsubara 
    for(int i=0;i<=ntau_;i++){
         double tau=dtau*(double)i;
        tmp(0,0)=G0_matsu(wa,beta,tau,sign_a);
        
        A_ref.set_mat(i,tmp);
    }
    
    //left-mixng
    for(int i2=0;i2<=ntau_;i2++){
        double tau2=dtau*(double)i2;
        for(int i1=0;i1<=nt_;i1++){
            double t1=h*(double)i1;
            tmp(0,0)=G0_tv(wa,beta,t1,tau2,sign_a);
            A_ref.set_tv(i1,i2,tmp);
        }
    }
    
    //lesser
    //[note:cntr only store i1<=i2]
    for(int i2=0;i2<=nt_;i2++){
        double t2=h*(double)i2;
        for(int i1=0;i1<=nt_;i1++){
            double t1=h*(double)i1;
            if(i1<=i2){
                tmp(0,0)=G0_les(wa,beta,t1,t2,sign_a);
                A_ref.set_les(i1,i2,tmp);
            }
        }
    }
    
    //reteraded
    for(int i2=0;i2<=nt_;i2++){
        double t2=h*(double)i2;
        for(int i1=0;i1<=nt_;i1++){
            double t1=h*(double)i1;
            if(i1>=i2){
                tmp(0,0)=G0_gtr(wa,beta,t1,t2,sign_a)-G0_les(wa,beta,t1,t2,sign_a);
                A_ref.set_ret(i1,i2,tmp);
            }
        }
    }
    
    
    
}

///////////////Set free equilibrium Green's functions of NxN ////////////////
// eps=U^\dagger Da U
void set_G0_ref_sizeN(GREEN & A_ref, double beta, double h, cdmatrix Ua, dvector Da, int sign_a){
    
    int size=A_ref.size1();
    assert(size==Ua.rows());
    assert(size==Da.size());
    
    int nt_=A_ref.nt();
    int ntau_=A_ref.ntau();
    double dtau=beta/(double)ntau_;
    
    cdmatrix tmp(size,size);
    
    //Matsubara 
    for(int i=0;i<=ntau_;i++){
        double tau=dtau*(double)i;
        tmp.setZero(size,size);
        
        for(int l=0;l<size;l++) {tmp(l,l)=G0_matsu(Da(l),beta,tau,sign_a);}
        tmp=Ua.adjoint()*tmp*Ua;
        A_ref.set_mat(i,tmp);
    }
    
    //left-mixng
    for(int i2=0;i2<=ntau_;i2++){
        double tau2=dtau*(double)i2;
        for(int i1=0;i1<=nt_;i1++){
            double t1=h*(double)i1;
            tmp.setZero(size,size);
            
            for(int l=0;l<size;l++) {tmp(l,l)=G0_tv(Da(l),beta,t1,tau2,sign_a);}
            tmp=Ua.adjoint()*tmp*Ua;
            A_ref.set_tv(i1,i2,tmp);
        }
    }
    //lesser
    //[note:cntr only store i1<=i2]
    for(int i2=0;i2<=nt_;i2++){
        double t2=h*(double)i2;
        for(int i1=0;i1<=nt_;i1++){
            double t1=h*(double)i1;
            tmp.setZero(size,size);
            
            if(i1<=i2){
                for(int l=0;l<size;l++) {tmp(l,l)=G0_les(Da(l),beta,t1,t2,sign_a);}
                tmp=Ua.adjoint()*tmp*Ua;
                A_ref.set_les(i1,i2,tmp);
            }
        }
    }
    
    //rteraded
    for(int i2=0;i2<=nt_;i2++){
        double t2=h*(double)i2;
        for(int i1=0;i1<=nt_;i1++){
            double t1=h*(double)i1;
            tmp.setZero(size,size);
            
            if(i1>=i2){
                for(int l=0;l<size;l++) {tmp(l,l)=G0_gtr(Da(l),beta,t1,t2,sign_a)-G0_les(Da(l),beta,t1,t2,sign_a);}
                tmp=Ua.adjoint()*tmp*Ua;
                A_ref.set_ret(i1,i2,tmp);
            }
        }
    }
    
    
}

///////////////Set reference bubble of 1x1 ////////////////

void get_bubble1_ref_size1(GREEN & C_ref, double beta, double h, double wa, double wb, int sign_a,int sign_b){
    
    assert (C_ref.sig()==sign_a*sign_b);
    
    int nt_=C_ref.nt();
    int ntau_=C_ref.ntau();
    double dtau=beta/(double)ntau_;
    
    cdmatrix tmp(1,1);

    //Matsubara 
    for(int i=0;i<=ntau_;i++){
        double tau=dtau*(double)i;
        //be care ful for i==0 
        if(i==0){
            if(sign_b==-1){
                tmp(0,0)=G0_matsu(wa,beta,tau,sign_a)*G0_matsu(wb,beta,beta,sign_b);
            }
            else if(sign_b==1){
                tmp(0,0)=-G0_matsu(wa,beta,tau,sign_a)*G0_matsu(wb,beta,beta,sign_b);
            }
        }else{
            tmp(0,0)=-G0_matsu(wa,beta,tau,sign_a)*G0_matsu(wb,beta,-tau,sign_b);
        }
        C_ref.set_mat(i,tmp);
    }

    //left-mixng
    for(int i2=0;i2<=ntau_;i2++){
        double tau2=dtau*(double)i2;
        for(int i1=0;i1<=nt_;i1++){
            double t1=h*(double)i1;
            tmp(0,0)=xj*G0_tv(wa,beta,t1,tau2,sign_a)*G0_vt(wb,beta,tau2,t1,sign_b);
            C_ref.set_tv(i1,i2,tmp);
        }
    }

    //lesser
    //[note:cntr only store i1<=i2]
    for(int i2=0;i2<=nt_;i2++){
        double t2=h*(double)i2;
        for(int i1=0;i1<=nt_;i1++){
            double t1=h*(double)i1;
            if(i1<=i2){
                tmp(0,0)=xj*G0_les(wa,beta,t1,t2,sign_a)*G0_gtr(wb,beta,t2,t1,sign_b);
                C_ref.set_les(i1,i2,tmp);
            }
        }
    }

    //rteraded
    for(int i2=0;i2<=nt_;i2++){
        double t2=h*(double)i2;
        for(int i1=0;i1<=nt_;i1++){
            double t1=h*(double)i1;
            if(i1>=i2){
                tmp(0,0)=xj*G0_gtr(wa,beta,t1,t2,sign_a)*G0_les(wb,beta,t2,t1,sign_b);
                tmp(0,0)-=xj*G0_les(wa,beta,t1,t2,sign_a)*G0_gtr(wb,beta,t2,t1,sign_b);
                C_ref.set_ret(i1,i2,tmp);
           }
        }
    }

}

void get_bubble2_ref_size1(GREEN & C_ref, double beta, double h, double wa, double wb, int sign_a,int sign_b){
    
    assert (C_ref.sig()==sign_a*sign_b);
    
    int nt_=C_ref.nt();
    int ntau_=C_ref.ntau();
    int tstp=-1;
    double dtau=beta/(double)ntau_;
    
    cdmatrix tmp(1,1);
    
    //Matsubara 
    for(int i=0;i<=ntau_;i++){
        double tau=dtau*(double)i;
        tmp(0,0)=-G0_matsu(wa,beta,tau,sign_a)*G0_matsu(wb,beta,tau,sign_b);
        C_ref.set_mat(i,tmp);
    }
    
    //left-mixng
    for(int i2=0;i2<=ntau_;i2++){
        double tau2=dtau*(double)i2;
        for(int i1=0;i1<=nt_;i1++){
            double t1=h*(double)i1;
            tmp(0,0)=xj*G0_tv(wa,beta,t1,tau2,sign_a)*G0_tv(wb,beta,t1,tau2,sign_b);
            C_ref.set_tv(i1,i2,tmp);
        }
    }
    
    //lesser
    //[note:cntr only store i1<=i2]
    for(int i2=0;i2<=nt_;i2++){
        double t2=h*(double)i2;
        for(int i1=0;i1<=nt_;i1++){
            double t1=h*(double)i1;
            if(i1<=i2){
                tmp(0,0)=xj*G0_les(wa,beta,t1,t2,sign_a)*G0_les(wb,beta,t1,t2,sign_b);
                C_ref.set_les(i1,i2,tmp);
            }
        }
    }
    
    //rteraded
    for(int i2=0;i2<=nt_;i2++){
        double t2=h*(double)i2;
        for(int i1=0;i1<=nt_;i1++){
            double t1=h*(double)i1;
            if(i1>=i2){
                tmp(0,0)=xj*G0_gtr(wa,beta,t1,t2,sign_a)*G0_gtr(wb,beta,t1,t2,sign_b);
                tmp(0,0)-=xj*G0_les(wa,beta,t1,t2,sign_a)*G0_les(wb,beta,t1,t2,sign_b);
                C_ref.set_ret(i1,i2,tmp);
            }
        }
    }
    
}

///////////////Set reference bubble of NxN ////////////////
//C_{ic1,ic2}(t,t')=iA_{ia1,ia2}(t,t')*B_{ib2,ib1}(t',t)
void get_bubble1_ref_sizeN(GREEN & C_ref, int ic1, int ic2, cdmatrix Ua, dvector Da, int sign_a, int ia1, int ia2, cdmatrix Ub, dvector Db, int sign_b, int ib1, int ib2,
                           double beta, double h){
    
    assert (C_ref.sig()==sign_a*sign_b);
    
    int nt_=C_ref.nt();
    int ntau_=C_ref.ntau();
    double dtau=beta/(double)ntau_;
    
    cdmatrix tmp;
    
    //Matsubara
    for(int i=0;i<=ntau_;i++){
        double tau=dtau*(double)i;
        
        C_ref.get_mat(i,tmp);
        
        if(i==0){
            if(sign_b==-1){
                tmp(ic1,ic2)=G0_matsu(Ua,Da,beta,tau,sign_a,ia1,ia2)*G0_matsu(Ub,Db,beta,beta,sign_b,ib2,ib1);
            }
            else if(sign_b==1){
                tmp(ic1,ic2)=-G0_matsu(Ua,Da,beta,tau,sign_a,ia1,ia2)*G0_matsu(Ub,Db,beta,beta,sign_b,ib2,ib1);
            }
        }else{
            tmp(ic1,ic2)=-G0_matsu(Ua,Da,beta,tau,sign_a,ia1,ia2)*G0_matsu(Ub,Db,beta,-tau,sign_b,ib2,ib1);
        }
        
        C_ref.set_mat(i,tmp);
    }
    
    //left-mixng
    for(int i2=0;i2<=ntau_;i2++){
        double tau2=dtau*(double)i2;
        for(int i1=0;i1<=nt_;i1++){
            double t1=h*(double)i1;
            
            C_ref.get_tv(i1,i2,tmp);
            tmp(ic1,ic2)=xj*G0_tv(Ua,Da,beta,t1,tau2,sign_a,ia1,ia2)*G0_vt(Ub,Db,beta,tau2,t1,sign_b,ib2,ib1);
            C_ref.set_tv(i1,i2,tmp);
        }
    }
    
    //lesser
    //[note:cntr only store i1<=i2]
    for(int i2=0;i2<=nt_;i2++){
        double t2=h*(double)i2;
        for(int i1=0;i1<=nt_;i1++){
            double t1=h*(double)i1;
            if(i1<=i2){
                C_ref.get_les(i1,i2,tmp);
                tmp(ic1,ic2)=xj*G0_les(Ua,Da,beta,t1,t2,sign_a,ia1,ia2)*G0_gtr(Ub,Db,beta,t2,t1,sign_b,ib2,ib1);
                C_ref.set_les(i1,i2,tmp);
            }
        }
    }
    
    //rteraded
    for(int i2=0;i2<=nt_;i2++){
        double t2=h*(double)i2;
        for(int i1=0;i1<=nt_;i1++){
            double t1=h*(double)i1;
            if(i1>=i2){
                C_ref.get_ret(i1,i2,tmp);
                tmp(ic1,ic2)=xj*G0_gtr(Ua,Da,beta,t1,t2,sign_a,ia1,ia2)*G0_les(Ub,Db,beta,t2,t1,sign_b,ib2,ib1);
                tmp(ic1,ic2)-=xj*G0_les(Ua,Da,beta,t1,t2,sign_a,ia1,ia2)*G0_gtr(Ub,Db,beta,t2,t1,sign_b,ib2,ib1);
                C_ref.set_ret(i1,i2,tmp);
            }
        }
    }
    
}

void get_bubble2_ref_sizeN(GREEN & C_ref, int ic1, int ic2, cdmatrix Ua, dvector Da, int sign_a, int ia1, int ia2, cdmatrix Ub, dvector Db, int sign_b, int ib1, int ib2,
                           double beta, double h){
    
    assert (C_ref.sig()==sign_a*sign_b);
    
    int nt_=C_ref.nt();
    int ntau_=C_ref.ntau();
    double dtau=beta/(double)ntau_;
    
    cdmatrix tmp;
    
    //Matsubara 
    for(int i=0;i<=ntau_;i++){
        double tau=dtau*(double)i;
        
        C_ref.get_mat(i,tmp);
        tmp(ic1,ic2)=-G0_matsu(Ua,Da,beta,tau,sign_a,ia1,ia2)*G0_matsu(Ub,Db,beta,tau,sign_b,ib1,ib2);
        C_ref.set_mat(i,tmp);
    }
    
    //left-mixng
    for(int i2=0;i2<=ntau_;i2++){
        double tau2=dtau*(double)i2;
        for(int i1=0;i1<=nt_;i1++){
            double t1=h*(double)i1;
            
            C_ref.get_tv(i1,i2,tmp);
            tmp(ic1,ic2)=xj*G0_tv(Ua,Da,beta,t1,tau2,sign_a,ia1,ia2)*G0_tv(Ub,Db,beta,t1,tau2,sign_b,ib1,ib2);
            C_ref.set_tv(i1,i2,tmp);
        }
    }
    
    //lesser
    //[note:cntr only store i1<=i2]
    for(int i2=0;i2<=nt_;i2++){
        double t2=h*(double)i2;
        for(int i1=0;i1<=nt_;i1++){
            double t1=h*(double)i1;
            if(i1<=i2){
                C_ref.get_les(i1,i2,tmp);
                tmp(ic1,ic2)=xj*G0_les(Ua,Da,beta,t1,t2,sign_a,ia1,ia2)*G0_les(Ub,Db,beta,t1,t2,sign_b,ib1,ib2);
                C_ref.set_les(i1,i2,tmp);
            }
        }
    }
    
    //rteraded
    for(int i2=0;i2<=nt_;i2++){
        double t2=h*(double)i2;
        for(int i1=0;i1<=nt_;i1++){
            double t1=h*(double)i1;
            if(i1>=i2){
                C_ref.get_ret(i1,i2,tmp);
                tmp(ic1,ic2)=xj*G0_gtr(Ua,Da,beta,t1,t2,sign_a,ia1,ia2)*G0_gtr(Ub,Db,beta,t1,t2,sign_b,ib1,ib2);
                tmp(ic1,ic2)-=xj*G0_les(Ua,Da,beta,t1,t2,sign_a,ia1,ia2)*G0_les(Ub,Db,beta,t1,t2,sign_b,ib1,ib2);
                C_ref.set_ret(i1,i2,tmp);
            }
        }
    }
    
}

/////////////////////////////////////////////
TEST_CASE("Bubble:size1","[Bubble:size1]"){
    
    
    int nt,ntau,kt,size,tstp;
	double wa,wb,h,beta,dtau;
    double eps,err;
    cdmatrix eps_a(1,1),eps_b(1,1);
    
    GREEN A_fermi,B_fermi,C_fermi;
    GREEN A_bose,B_bose,C_bose;
    
    GREEN C_fermi_ref;
    GREEN C_bose_ref;
    
    
    eps=1e-8;
    
    nt=200;
	ntau=500;
    beta=0.5;
    dtau=beta/(double)ntau;
    h=0.02;
    wa=1.123;
	wb=0.345;
    eps_a(0,0)=1.123;
    eps_b(0,0)=0.345;
	kt=5;
    
    size=1;
    

    SECTION ("FF bubbe1"){
        
        A_fermi=GREEN(nt,ntau,size,-1);
        B_fermi=GREEN(nt,ntau,size,-1);
        //cntr::green_from_H(A_fermi,0.0,eps_a,beta,h);
        set_G0_ref_size1(A_fermi, beta, h, wa, -1);
        //cntr::green_from_H(B_fermi,0.0,eps_b,beta,h);
        set_G0_ref_size1(B_fermi, beta, h, wb, -1);
        
        C_bose=GREEN(nt,ntau,size,1);
        C_bose_ref=GREEN(nt,ntau,size,1);
        
        err=0.0;
        get_bubble1_ref_size1(C_bose_ref, beta, h, wa, wb, -1, -1);

        for(tstp=-1;tstp<=nt;tstp++){
            cntr::Bubble1(tstp,C_bose,A_fermi,B_fermi);
            
            err+=cntr::distance_norm2(tstp,C_bose,C_bose_ref);
            //cout<<tstp<<" "<<cntr::distance_norm2(tstp,C_bose,C_bose_ref)<<endl;
            
        }
        REQUIRE(err<eps);
    }
    
    
    SECTION ("FF bubbe2"){
     
        A_fermi=GREEN(nt,ntau,size,-1);
        B_fermi=GREEN(nt,ntau,size,-1);
        //cntr::green_from_H(A_fermi,0.0,eps_a,beta,h);
        set_G0_ref_size1(A_fermi, beta, h, wa, -1);
        //cntr::green_from_H(B_fermi,0.0,eps_b,beta,h);
        set_G0_ref_size1(B_fermi, beta, h, wb, -1);
        
        C_bose=GREEN(nt,ntau,size,1);
        C_bose_ref=GREEN(nt,ntau,size,1);
        
        err=0.0;
        get_bubble2_ref_size1(C_bose_ref, beta, h, wa, wb, -1, -1);
        
        for(tstp=-1;tstp<=nt;tstp++){
            cntr::Bubble2(tstp,C_bose,A_fermi,B_fermi);
            
            err+=cntr::distance_norm2(tstp,C_bose,C_bose_ref);
            //cout<<tstp<<" "<<cntr::distance_norm2(tstp,C_bose,C_bose_ref)<<endl;
            
        }
        REQUIRE(err<eps);
     
     /*
     if(tstp==-1){
     for(int m=0;m<=ntau;m++){
     double tau=dtau*(double)m;
     CPLX ctmp,ctmp2;
     A_fermi.get_mat(m,ctmp);
     B_fermi.get_mat(ntau-m,ctmp2);
     cout<<"A,B martin:"<<tau<<" "<<ctmp<<" "<<ctmp2<<endl;
     
     cout<<"A,B heree:"<<tau<<" "<<G0_matsu(wa,beta,tau,-1)<<" "<<G0_matsu(wb,beta,beta-tau,-1)<<endl;
     
     C_bose.get_mat(m,ctmp);
     C_bose_ref.get_mat(m,ctmp2);
     cout<<"C:"<<tau<<" "<<ctmp<<" "<<ctmp2<<endl;
     }
     }*/
     
    }
    
    SECTION ("FB bubbe1"){
        
        A_fermi=GREEN(nt,ntau,size,-1);
        set_G0_ref_size1(A_fermi, beta, h, wa, -1);
        B_bose=GREEN(nt,ntau,size,1);
        set_G0_ref_size1(B_bose, beta, h, wb, 1);
        
        C_fermi=GREEN(nt,ntau,size,-1);
        C_fermi_ref=GREEN(nt,ntau,size,-1);
        
        err=0.0;
        
        get_bubble1_ref_size1(C_fermi_ref, beta, h, wa, wb, -1, 1);
        
        for(tstp=-1;tstp<=nt;tstp++){
            cntr::Bubble1(tstp,C_fermi,A_fermi,B_bose);
            
            err+=cntr::distance_norm2(tstp,C_fermi,C_fermi_ref);
            
        }
        REQUIRE(err<eps);
        
    }
    
    SECTION ("FB bubbe2"){
        
        A_fermi=GREEN(nt,ntau,size,-1);
        set_G0_ref_size1(A_fermi, beta, h, wa, -1);
        B_bose=GREEN(nt,ntau,size,1);
        set_G0_ref_size1(B_bose, beta, h, wb, 1);
        
        C_fermi=GREEN(nt,ntau,size,-1);
        C_fermi_ref=GREEN(nt,ntau,size,-1);
        
        err=0.0;
        
        get_bubble2_ref_size1(C_fermi_ref, beta, h, wa, wb, -1, 1);
        
        for(tstp=-1;tstp<=nt;tstp++){
            cntr::Bubble2(tstp,C_fermi,A_fermi,B_bose);
            
            err+=cntr::distance_norm2(tstp,C_fermi,C_fermi_ref);
            
        }
        REQUIRE(err<eps);
        
    }
    
    SECTION ("BF bubbe1"){
        
        A_bose=GREEN(nt,ntau,size,1);
        set_G0_ref_size1(A_bose, beta, h, wa, 1);
        B_fermi=GREEN(nt,ntau,size,-1);
        set_G0_ref_size1(B_fermi, beta, h, wb, -1);
        
        C_fermi=GREEN(nt,ntau,size,-1);
        C_fermi_ref=GREEN(nt,ntau,size,-1);
        
        err=0.0;
        
        get_bubble1_ref_size1(C_fermi_ref, beta, h, wa, wb, 1, -1);
        
        for(tstp=-1;tstp<=nt;tstp++){
            cntr::Bubble1(tstp,C_fermi,A_bose,B_fermi);
            
            err+=cntr::distance_norm2(tstp,C_fermi,C_fermi_ref);
            
        }
        REQUIRE(err<eps);
    }
    
    
    SECTION ("BF bubbe2"){
        
        A_bose=GREEN(nt,ntau,size,1);
        set_G0_ref_size1(A_bose, beta, h, wa, 1);
        B_fermi=GREEN(nt,ntau,size,-1);
        set_G0_ref_size1(B_fermi, beta, h, wb, -1);
        
        C_fermi=GREEN(nt,ntau,size,-1);
        C_fermi_ref=GREEN(nt,ntau,size,-1);
        
        err=0.0;
        
        get_bubble2_ref_size1(C_fermi_ref, beta, h, wa, wb, 1, -1);
        
        for(tstp=-1;tstp<=nt;tstp++){
            cntr::Bubble2(tstp,C_fermi,A_bose,B_fermi);
            
            err+=cntr::distance_norm2(tstp,C_fermi,C_fermi_ref);
            
        }
        REQUIRE(err<eps);
        
    }
    
    SECTION ("BB bubbe1"){
        
        A_bose=GREEN(nt,ntau,size,1);
        set_G0_ref_size1(A_bose, beta, h, wa, 1);
        B_bose=GREEN(nt,ntau,size,1);
        set_G0_ref_size1(B_bose, beta, h, wb, 1);
        
        C_bose=GREEN(nt,ntau,size,1);
        C_bose_ref=GREEN(nt,ntau,size,1);
        
        err=0.0;
        
        get_bubble1_ref_size1(C_bose_ref, beta, h, wa, wb, 1, 1);
        
        for(tstp=-1;tstp<=nt;tstp++){
            cntr::Bubble1(tstp,C_bose,A_bose,B_bose);
            
            err+=cntr::distance_norm2(tstp,C_bose,C_bose_ref);
            
        }
        REQUIRE(err<eps);
        
    }
    
    
    SECTION ("BB bubbe2"){
        
        A_bose=GREEN(nt,ntau,size,1);
        set_G0_ref_size1(A_bose, beta, h, wa, 1);
        B_bose=GREEN(nt,ntau,size,1);
        set_G0_ref_size1(B_bose, beta, h, wb, 1);
        
        C_bose=GREEN(nt,ntau,size,1);
        C_bose_ref=GREEN(nt,ntau,size,1);
        
        err=0.0;
        
        get_bubble2_ref_size1(C_bose_ref, beta, h, wa, wb, 1, 1);
        
        for(tstp=-1;tstp<=nt;tstp++){
            cntr::Bubble2(tstp,C_bose,A_bose,B_bose);
            
            err+=cntr::distance_norm2(tstp,C_bose,C_bose_ref);
            
        }
        REQUIRE(err<eps);
        
        
    }
    
}


TEST_CASE("Bubble:sizeN","[Bubble:sizeN]"){
    
    
    int nt_,ntau_,kt_,tstp_;
	double h_,beta_,dtau_;
    double eps_,err_;
    
    int size_=2;
    dvector Ea_(size_),Eb_(size_);
    Ea_<<1.123, 0.567;
    Eb_<<0.345, 0.876;
    
    cdmatrix Ua_(size_,size_),Ub_(size_,size_);
    
    Ua_<< cos(1.0), -sin(1.0),
         sin(1.0), cos(1.0);
    
    Ub_<< cos(2.0), -sin(2.0),
        sin(2.0), cos(2.0);
    
    GREEN A_,B_,C_;
    GREEN C_ref_;
    
    eps_=1e-8;
    
    nt_=50;
	ntau_=50;
    beta_=0.5;
    dtau_=beta_/(double)ntau_;
    h_=0.02;
	kt_=5;
    

    
    
    SECTION ("FF bubbe1"){
        
        int sign_a_=-1;
        int sign_b_=-1;
        int sign_c_=sign_a_*sign_b_;
        
        A_=GREEN(nt_,ntau_,size_,sign_a_);
        B_=GREEN(nt_,ntau_,size_,sign_b_);
        
        set_G0_ref_sizeN(A_, beta_, h_, Ua_,Ea_, sign_a_);
        set_G0_ref_sizeN(B_, beta_, h_, Ub_,Eb_, sign_b_);
        
        C_=GREEN(nt_,ntau_,size_,sign_c_);
        C_ref_=GREEN(nt_,ntau_,size_,sign_c_);
        
        err_=0.0;
        for(int ic1=0;ic1<size_;ic1++){
            for(int ic2=0;ic2<size_;ic2++){
            get_bubble1_ref_sizeN(C_ref_, ic1, ic2,  Ua_, Ea_, sign_a_, ic1, ic2, Ub_, Eb_, sign_b_, ic1, ic2,beta_, h_);
            }
        }
        
        for(tstp_=-1;tstp_<=nt_;tstp_++){
            
            for(int ic1=0;ic1<size_;ic1++){
                for(int ic2=0;ic2<size_;ic2++){
                    cntr::Bubble1(tstp_,C_,ic1,ic2,A_,A_,ic1,ic2,B_,B_,ic1,ic2);
                }
            }
            
            err_+=cntr::distance_norm2(tstp_,C_,C_ref_);
            //cout<<tstp<<" "<<cntr::distance_norm2(tstp,C_bose,C_bose_ref)<<endl;
            
        }
        REQUIRE(err_<eps_);
    }
    
    
    SECTION ("FF bubbe2"){
        
        int sign_a_=-1;
        int sign_b_=-1;
        int sign_c_=sign_a_*sign_b_;
        
        A_=GREEN(nt_,ntau_,size_,sign_a_);
        B_=GREEN(nt_,ntau_,size_,sign_b_);

        set_G0_ref_sizeN(A_, beta_, h_, Ua_,Ea_, sign_a_);
        set_G0_ref_sizeN(B_, beta_, h_, Ub_,Eb_, sign_b_);
                
        C_=GREEN(nt_,ntau_,size_,sign_c_);
        C_ref_=GREEN(nt_,ntau_,size_,sign_c_);
        
        err_=0.0;
        
        for(int ic1=0;ic1<size_;ic1++){
            for(int ic2=0;ic2<size_;ic2++){
                get_bubble2_ref_sizeN(C_ref_, ic1, ic2,  Ua_, Ea_, sign_a_, ic1, ic2, Ub_, Eb_, sign_b_, ic1, ic2,beta_, h_);
            }
        }
        
        for(tstp_=-1;tstp_<=nt_;tstp_++){
            
            for(int ic1=0;ic1<size_;ic1++){
                for(int ic2=0;ic2<size_;ic2++){
                    cntr::Bubble2(tstp_,C_,ic1,ic2,A_,A_,ic1,ic2,B_,B_,ic1,ic2);
                }
            }
            
            err_+=cntr::distance_norm2(tstp_,C_,C_ref_);
            
        }
        REQUIRE(err_<eps_);
    
    }
    
    SECTION ("FB bubbe1"){
        
        int sign_a_=-1;
        int sign_b_=+1;
        int sign_c_=sign_a_*sign_b_;
        
        A_=GREEN(nt_,ntau_,size_,sign_a_);
        B_=GREEN(nt_,ntau_,size_,sign_b_);
        
        set_G0_ref_sizeN(A_, beta_, h_, Ua_,Ea_, sign_a_);
        set_G0_ref_sizeN(B_, beta_, h_, Ub_,Eb_, sign_b_);
        
        C_=GREEN(nt_,ntau_,size_,sign_c_);
        C_ref_=GREEN(nt_,ntau_,size_,sign_c_);
        
        err_=0.0;
        
        for(int ic1=0;ic1<size_;ic1++){
            for(int ic2=0;ic2<size_;ic2++){
                get_bubble1_ref_sizeN(C_ref_, ic1, ic2,  Ua_, Ea_, sign_a_, ic1, ic2, Ub_, Eb_, sign_b_, ic1, ic2,beta_, h_);
            }
        }
        
        for(tstp_=-1;tstp_<=nt_;tstp_++){
            
            for(int ic1=0;ic1<size_;ic1++){
                for(int ic2=0;ic2<size_;ic2++){
                    cntr::Bubble1(tstp_,C_,ic1,ic2,A_,A_,ic1,ic2,B_,B_,ic1,ic2);
                }
            }
            
            err_+=cntr::distance_norm2(tstp_,C_,C_ref_);
            //cout<<tstp<<" "<<cntr::distance_norm2(tstp,C_bose,C_bose_ref)<<endl;
            
        }
        REQUIRE(err_<eps_);
        
    }
    
    SECTION ("FB bubbe2"){
        
        int sign_a_=-1;
        int sign_b_=+1;
        int sign_c_=sign_a_*sign_b_;
        
        A_=GREEN(nt_,ntau_,size_,sign_a_);
        B_=GREEN(nt_,ntau_,size_,sign_b_);
        
        set_G0_ref_sizeN(A_, beta_, h_, Ua_,Ea_, sign_a_);
        set_G0_ref_sizeN(B_, beta_, h_, Ub_,Eb_, sign_b_);
        
        C_=GREEN(nt_,ntau_,size_,sign_c_);
        C_ref_=GREEN(nt_,ntau_,size_,sign_c_);
        
        err_=0.0;
        
        for(int ic1=0;ic1<size_;ic1++){
            for(int ic2=0;ic2<size_;ic2++){
                get_bubble2_ref_sizeN(C_ref_, ic1, ic2,  Ua_, Ea_, sign_a_, ic1, ic2, Ub_, Eb_, sign_b_, ic1, ic2,beta_, h_);
            }
        }

        for(tstp_=-1;tstp_<=nt_;tstp_++){
            
            for(int ic1=0;ic1<size_;ic1++){
                for(int ic2=0;ic2<size_;ic2++){
                    cntr::Bubble2(tstp_,C_,ic1,ic2,A_,A_,ic1,ic2,B_,B_,ic1,ic2);
                }
            }
            
            err_+=cntr::distance_norm2(tstp_,C_,C_ref_);
            
        }
        REQUIRE(err_<eps_);
        
    }
    
    SECTION ("BF bubbe1"){
        
        int sign_a_=+1;
        int sign_b_=-1;
        int sign_c_=sign_a_*sign_b_;
        
        A_=GREEN(nt_,ntau_,size_,sign_a_);
        B_=GREEN(nt_,ntau_,size_,sign_b_);
        
        set_G0_ref_sizeN(A_, beta_, h_, Ua_,Ea_, sign_a_);
        set_G0_ref_sizeN(B_, beta_, h_, Ub_,Eb_, sign_b_);
        
        C_=GREEN(nt_,ntau_,size_,sign_c_);
        C_ref_=GREEN(nt_,ntau_,size_,sign_c_);
        
        err_=0.0;
        
        for(int ic1=0;ic1<size_;ic1++){
            for(int ic2=0;ic2<size_;ic2++){
                get_bubble1_ref_sizeN(C_ref_, ic1, ic2,  Ua_, Ea_, sign_a_, ic1, ic2, Ub_, Eb_, sign_b_, ic1, ic2,beta_, h_);
            }
        }
        
        for(tstp_=-1;tstp_<=nt_;tstp_++){
            
            for(int ic1=0;ic1<size_;ic1++){
                for(int ic2=0;ic2<size_;ic2++){
                    cntr::Bubble1(tstp_,C_,ic1,ic2,A_,A_,ic1,ic2,B_,B_,ic1,ic2);
                }
            }
            
            err_+=cntr::distance_norm2(tstp_,C_,C_ref_);
            //cout<<tstp<<" "<<cntr::distance_norm2(tstp,C_bose,C_bose_ref)<<endl;
            
        }
        REQUIRE(err_<eps_);
    }
    
    
    SECTION ("BF bubbe2"){
        
        int sign_a_=+1;
        int sign_b_=-1;
        int sign_c_=sign_a_*sign_b_;
        
        A_=GREEN(nt_,ntau_,size_,sign_a_);
        B_=GREEN(nt_,ntau_,size_,sign_b_);
        
        set_G0_ref_sizeN(A_, beta_, h_, Ua_,Ea_, sign_a_);
        set_G0_ref_sizeN(B_, beta_, h_, Ub_,Eb_, sign_b_);
        
        C_=GREEN(nt_,ntau_,size_,sign_c_);
        C_ref_=GREEN(nt_,ntau_,size_,sign_c_);
        
        err_=0.0;
        
        for(int ic1=0;ic1<size_;ic1++){
            for(int ic2=0;ic2<size_;ic2++){
                get_bubble2_ref_sizeN(C_ref_, ic1, ic2,  Ua_, Ea_, sign_a_, ic1, ic2, Ub_, Eb_, sign_b_, ic1, ic2,beta_, h_);
            }
        }
        
        for(tstp_=-1;tstp_<=nt_;tstp_++){
            
            for(int ic1=0;ic1<size_;ic1++){
                for(int ic2=0;ic2<size_;ic2++){
                    cntr::Bubble2(tstp_,C_,ic1,ic2,A_,A_,ic1,ic2,B_,B_,ic1,ic2);
                }
            }
            
            err_+=cntr::distance_norm2(tstp_,C_,C_ref_);
            
        }
        REQUIRE(err_<eps_);
        
    }
    
    SECTION ("BB bubbe1"){
        
        int sign_a_=+1;
        int sign_b_=+1;
        int sign_c_=sign_a_*sign_b_;
        
        A_=GREEN(nt_,ntau_,size_,sign_a_);
        B_=GREEN(nt_,ntau_,size_,sign_b_);
        
        set_G0_ref_sizeN(A_, beta_, h_, Ua_,Ea_, sign_a_);
        set_G0_ref_sizeN(B_, beta_, h_, Ub_,Eb_, sign_b_);
        
        C_=GREEN(nt_,ntau_,size_,sign_c_);
        C_ref_=GREEN(nt_,ntau_,size_,sign_c_);
        
        err_=0.0;
        
        for(int ic1=0;ic1<size_;ic1++){
            for(int ic2=0;ic2<size_;ic2++){
                get_bubble1_ref_sizeN(C_ref_, ic1, ic2,  Ua_, Ea_, sign_a_, ic1, ic2, Ub_, Eb_, sign_b_, ic1, ic2,beta_, h_);
            }
        }
        
        for(tstp_=-1;tstp_<=nt_;tstp_++){
            
            for(int ic1=0;ic1<size_;ic1++){
                for(int ic2=0;ic2<size_;ic2++){
                    cntr::Bubble1(tstp_,C_,ic1,ic2,A_,A_,ic1,ic2,B_,B_,ic1,ic2);
                }
            }
            
            err_+=cntr::distance_norm2(tstp_,C_,C_ref_);
            //cout<<tstp<<" "<<cntr::distance_norm2(tstp,C_bose,C_bose_ref)<<endl;
            
        }
        REQUIRE(err_<eps_);
        
    }
    
    
    SECTION ("BB bubbe2"){
        
        int sign_a_=+1;
        int sign_b_=+1;
        int sign_c_=sign_a_*sign_b_;
        
        A_=GREEN(nt_,ntau_,size_,sign_a_);
        B_=GREEN(nt_,ntau_,size_,sign_b_);
        
        set_G0_ref_sizeN(A_, beta_, h_, Ua_,Ea_, sign_a_);
        set_G0_ref_sizeN(B_, beta_, h_, Ub_,Eb_, sign_b_);
        
        C_=GREEN(nt_,ntau_,size_,sign_c_);
        C_ref_=GREEN(nt_,ntau_,size_,sign_c_);
        
        err_=0.0;
        
        for(int ic1=0;ic1<size_;ic1++){
            for(int ic2=0;ic2<size_;ic2++){
                get_bubble2_ref_sizeN(C_ref_, ic1, ic2,  Ua_, Ea_, sign_a_, ic1, ic2, Ub_, Eb_, sign_b_, ic1, ic2,beta_, h_);
            }
        }
        
        for(tstp_=-1;tstp_<=nt_;tstp_++){
            
            for(int ic1=0;ic1<size_;ic1++){
                for(int ic2=0;ic2<size_;ic2++){
                    cntr::Bubble2(tstp_,C_,ic1,ic2,A_,A_,ic1,ic2,B_,B_,ic1,ic2);
                }
            }
            
            err_+=cntr::distance_norm2(tstp_,C_,C_ref_);
            
        }
        REQUIRE(err_<eps_);
        
    }
    
}
