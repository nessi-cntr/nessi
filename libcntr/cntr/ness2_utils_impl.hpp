#include "ness2.hpp"
#include <fstream>


namespace ness2
{

  
   /// @private
  double  distance_norm2(fft_array &A,fft_array &B,fft_domain type){
    if (A.Nft_ != B.Nft_ || A.size1_ != B.size1_) {
      throw std::invalid_argument("ness2:distance_norm2 mismatched array sizes");
    }
    fftw_complex* A_data = (type == fft_domain::time) ? A.time_ : A.freq_;
    fftw_complex* B_data = (type == fft_domain::time) ? B.time_ : B.freq_;
    double x,sum=0.0;
    int total_size = A.Nft_ * A.element_size_;
    for (int i = 0; i < total_size; ++i) {
      x=A_data[i][0]-B_data[i][0];
      sum+=x*x;
      x=A_data[i][1]-B_data[i][1];
      sum+=x*x;
    }
    return sqrt(sum);
  }
  

  /// @private
  double  distance_norm2(herm_matrix_ness &A,herm_matrix_ness &B,fft_domain type){
    double x,sum=0.0;
    x=distance_norm2(A.ret_,B.ret_,type);
    sum+=x*x;
    x=distance_norm2(A.les_,B.les_,type);
    sum+=x*x;
    return sqrt(sum);
  }    

void Bubble1_ness(herm_matrix_ness &C,int c1,int c2,herm_matrix_ness &A,int a1,int a2,herm_matrix_ness &B,int b1,int b2)
{
  // TODO ASSERTS ETC ...    
    for(int t = 0; t < C.Nft_/2; t++){
      cdmatrix ales_t,aret_t,agtr_t;
      A.get_ret(t,aret_t,fft_domain::time); 
      A.get_les(t,ales_t,fft_domain::time); 
      agtr_t=aret_t+ales_t; 
      cdmatrix bles_t,bret_t,bgtr_t;
      B.get_ret(t,bret_t,fft_domain::time); 
      B.get_les(t,bles_t,fft_domain::time); 
      bgtr_t=bret_t+bles_t; 
      cdmatrix cles_t,cret_t;
      C.get_ret(t,cret_t,fft_domain::time); 
      C.get_les(t,cles_t,fft_domain::time); 
      cles_t(c1,c2) = NESS_II * ales_t(a1,a2) * (-bgtr_t.adjoint()(b2,b1));
      cdouble les12=cles_t(c1,c2);
      cret_t(c1,c2)= NESS_II * aret_t(a1,a2) * (-bles_t.adjoint()(b2,b1)) + NESS_II * ales_t(a1,a2) * bret_t.adjoint()(b2,b1);
      C.set_ret(t,cret_t,fft_domain::time); 
      C.set_les(t,cles_t,fft_domain::time);
      // for negative times: we set Cles_{c1,c2}(-t) = [- Cles(t)*cc]_{c1,c2} = - Cles(t)_{c2,c1}^*
      // such that C will be hermitian if all entries (c1,c2) are set
      if (t>0){
	C.get_les(-t,cles_t,fft_domain::time);
	cles_t(c2,c1) = -std::conj(les12);
	C.set_les(-t,cles_t,fft_domain::time);	  
      }
    }
    
}

void Bubble1_ness(herm_matrix_ness &C, herm_matrix_ness &A, herm_matrix_ness &B)
{
    return Bubble1_ness(C,0,0,A,0,0,B,0,0);
}


  void Bubble2_ness(herm_matrix_ness &C,int c1,int c2,herm_matrix_ness &A,int a1,int a2,herm_matrix_ness &B,int b1,int b2)
  {
    // TODO ASSERTS ETC ...    
    for(int t = 0; t < C.Nft_/2; t++){
      cdmatrix ales_t,aret_t,agtr_t;
      A.get_ret(t,aret_t,fft_domain::time); 
      A.get_les(t,ales_t,fft_domain::time); 
      agtr_t=aret_t+ales_t; 
      cdmatrix bles_t,bret_t,bgtr_t;
      B.get_ret(t,bret_t,fft_domain::time); 
      B.get_les(t,bles_t,fft_domain::time); 
      bgtr_t=bret_t+bles_t; 
      cdmatrix cles_t,cret_t;
      C.get_ret(t,cret_t,fft_domain::time); 
      C.get_les(t,cles_t,fft_domain::time);
      
      cles_t(c1,c2) = NESS_II * ales_t(a1,a2) * bles_t(b1,b2);
      cdouble les12=cles_t(c1,c2);
      cret_t(c1,c2)= NESS_II * aret_t(a1,a2) * bret_t(b1,b2) + NESS_II * ales_t(a1,a2) * bret_t(b1,b2)+ NESS_II * aret_t(a1,a2) * bles_t(b1,b2);
      C.set_ret(t,cret_t,fft_domain::time); 
      C.set_les(t,cles_t,fft_domain::time);
      // for negative times: we set Cles_{c1,c2}(-t) = [- Cles(t)*cc]_{c1,c2} = - Cles(t)_{c2,c1}^*
      // such that C will be hermitian if all entries (c1,c2) are set
      if (t>0){
	C.get_les(-t,cles_t,fft_domain::time);
	cles_t(c2,c1) = -std::conj(les12);
	C.set_les(-t,cles_t,fft_domain::time);	  
      }      
    }
  }
  
  
  ///@private
  void Bubble2_ness(herm_matrix_ness &C, herm_matrix_ness &A, herm_matrix_ness &B)
  {
    return Bubble2_ness(C,0,0,A,0,0,B,0,0);
  }
  
  
  void convolution_density_matrix(cdmatrix &result,int bosefermi,herm_matrix_ness &A, herm_matrix_ness &B,double h){
    result.resize(A.size1_, A.size1_);    
    int Nft = A.Nft_;
    int size = A.size1_;    
    if (B.Nft_!=Nft && B.size1_!=size) {
      throw std::runtime_error("convolution_density_matrix: input mismatch");
    }
    fft_array C(Nft,size);
    A.integral_transform_to_freq(h,FFT_TRAPEZ);
    B.integral_transform_to_freq(h,FFT_TRAPEZ);
    for (int w = -Nft/2; w < Nft/2; ++w) {
      cdmatrix A_les,A_ret,B_les,B_ret,B_adv,C_les;
      A.get_ret(w,A_ret,fft_domain::freq);
      A.get_les(w,A_les,fft_domain::freq);
      B.get_ret(w,B_ret,fft_domain::freq);
      B.get_les(w,B_les,fft_domain::freq);
      B_adv = B_ret.adjoint();      
      // C^< = A^< * B^A + A^R * B^<
      C_les = A_les * B_adv + A_ret * B_les;
      C.set_element(w,C_les,fft_domain::freq);
    }
    C.fft_to_time();
    grid_info grid(Nft,h);
    C.get_element(0,result,fft_domain::time);
    result *= std::complex<double>(0, bosefermi*grid.dw_*0.5/M_PI);
  }

  void density_matrix(cdmatrix &result, int bosefermi,herm_matrix_ness &A) {
    // Ensure result has correct size
    result.resize(A.size1_, A.size1_);    
    // Get the lesser component at t = 0 from the time-domain data
    A.get_les(0, result, fft_domain::time);    
    // Multiply by i
    result *= std::complex<double>(0.0, bosefermi*1.0);
  }

  fft_array downsample(fft_array &in, int factor){
    int Nft=in.Nft_;
    int Nft1=in.Nft_/factor;
    if (factor<1 || Nft%(2*factor)!=0) {
      throw std::invalid_argument("fft_array::downsample: wrong factor");
    }
    fft_array tmp(Nft1,in.size1_);
    cdmatrix M;
    for(int t=0;t<Nft1;t++){
      in.get_element(t*factor,M, fft_domain::time);
      tmp.set_element(t,M,fft_domain::time);
    }    
    for(int w=0;w<Nft1/2;w++){
      in.get_element(w,M, fft_domain::freq);
      tmp.set_element(w,M,fft_domain::freq);
    }
    for(int w=1;w<=Nft1/2;w++){
      in.get_element(Nft-w,M, fft_domain::freq);
      tmp.set_element(Nft1-w,M,fft_domain::freq);
    }
    return tmp;
  }

  herm_matrix_ness downsample(herm_matrix_ness &in, int factor){
    int Nft=in.Nft_;
    int Nft1=in.Nft_/factor;
    if (factor<1 || Nft%(2*factor)!=0) {
      throw std::invalid_argument("herm_matrix_ness::downsample: wrong factor");
    }
    herm_matrix_ness tmp(Nft1,in.size1_);
    tmp.set_zero(fft_domain::freq);
    tmp.set_zero(fft_domain::time);    
    cdmatrix M;
    for(int t=0;t<Nft1;t++){
      in.get_les(t*factor,M, fft_domain::time);
      tmp.set_les(t,M,fft_domain::time);
      in.get_ret(t*factor,M, fft_domain::time);
      tmp.set_ret(t,M,fft_domain::time);
    }        
    for(int w=0;w<Nft1/2;w++){
      in.get_ret(w,M, fft_domain::freq);
      tmp.set_ret(w,M,fft_domain::freq);
      in.get_les(w,M, fft_domain::freq);
      tmp.set_les(w,M,fft_domain::freq);
    }
    for(int w=1;w<=Nft1/2;w++){
      in.get_ret(Nft-w,M, fft_domain::freq);
      tmp.set_ret(Nft1-w,M,fft_domain::freq);
      in.get_les(Nft-w,M, fft_domain::freq);
      tmp.set_les(Nft1-w,M,fft_domain::freq);
    }
    return tmp;    
  }

  herm_matrix_ness upsample(double h_in,herm_matrix_ness &in, int factor){
    int Nft=in.Nft_;
    int Nft1=in.Nft_*factor;
    if (factor<1) {
      throw std::invalid_argument("herm_matrix_ness::upsample: wrong factor");
    }
    herm_matrix_ness tmp(Nft1,in.size1_);
    herm_matrix_ness tmp_in=in;
    tmp_in.integral_transform_to_freq(h_in,FFT_TRAPEZ);    
    double h=1.0;
    tmp.set_zero(fft_domain::freq);
    cdmatrix M;
    for(int w=0;w<Nft/2;w++){
      tmp_in.get_ret(w,M, fft_domain::freq);
      tmp.set_ret(w,M,fft_domain::freq);
      tmp_in.get_les(w,M, fft_domain::freq);
      tmp.set_les(w,M,fft_domain::freq);
    }
    for(int w=1;w<=Nft/2;w++){
      tmp_in.get_ret(Nft-w,M, fft_domain::freq);
      tmp.set_ret(Nft1-w,M,fft_domain::freq);
      tmp_in.get_les(Nft-w,M, fft_domain::freq);
      tmp.set_les(Nft1-w,M,fft_domain::freq);
    }    
    tmp.integral_transform_to_time(h_in/factor,FFT_TRAPEZ);
    return tmp;    
  }

  herm_matrix_ness resample(double h_in,herm_matrix_ness &in, int factor){
    if(factor<0) return upsample(h_in,in,-factor);
    else return downsample(in,factor);
  }

  
};//end of namespace
