#ifndef CNTR_CCLC_TSTP_IMPL
#define CNTR_CCLC_TSTP_IMPL

#include "cntr_cyclic_tstp_decl.hpp"

namespace cntr{
/* #######################################################################################
#
#   MEMBERS OF CYCLIC_TIMESTEP
#
########################################################################################*/
template <typename T> cyclic_timestep<T>::cyclic_timestep(){
   data_=0;
   ntau_=0;
   tstp_=0;
   size1_=0;
   size2_=0;
   element_size_=0;
   total_size_=0;
}
template <typename T> cyclic_timestep<T>::~cyclic_timestep(){
   delete [] data_;
}
template <typename T> cyclic_timestep<T>::cyclic_timestep(int tstp,int ntau,int size1){
   int len=((tstp+1)*2+(ntau+1))*size1*size1;
   assert(size1>=0 && tstp>=-1 && ntau>=0);
   if(len==0) data_=0;
   else{
      data_ = new cplx [len];
	  memset(data_, 0, sizeof(cplx)*len);
   }
   size1_=size1;
   size2_=size1;
   element_size_=size1*size1_;
   tstp_=tstp;
   ntau_=ntau;
   total_size_=len;
}
template <typename T> cyclic_timestep<T>::cyclic_timestep(const cyclic_timestep &g){
   tstp_=g.tstp_;
   ntau_=g.ntau_;
   size1_=g.size1_;
   size2_=g.size2_;
   element_size_=size1_*size2_;
   total_size_=g.total_size_;
   if(total_size_>0){
      data_ = new cplx [total_size_];
	  memset(data_, 0, sizeof(cplx)*total_size_);
   }else{
	   data_=0;
   }
}
template <typename T> cyclic_timestep<T> & cyclic_timestep<T>::operator=(const cyclic_timestep &g){
	if(this==&g) return *this;
	if( total_size_!=g.total_size_){
	   delete [] data_;
       total_size_=g.total_size_;
	   if(total_size_>0) data_ = new cplx [total_size_];
       else data_=0;
	}
	tstp_=g.tstp_;
    ntau_=g.ntau_;
    size1_=g.size1_;
    size2_=g.size2_;
    element_size_=size1_*size2_;
	if(total_size_>0) memcpy(data_, g.data_, sizeof(cplx)*total_size_);
    return *this;
}
template <typename T> void cyclic_timestep<T>::resize(int tstp,int ntau,int size1){
   int len=((tstp+1)*2+(ntau+1))*size1*size1;
   assert(ntau>=0 && tstp>=-1 && size1>=0);
   delete [] data_;
   if(len==0) data_=0;
   else{
      data_ = new cplx [len];
	  memset(data_, 0, sizeof(cplx)*len);
   }
   size1_=size1;
   size2_=size1;
   element_size_=size1_*size1_;
   tstp_=tstp;
   ntau_=ntau;
   total_size_=len;
}
template <typename T> void cyclic_timestep<T>::clear(void){
	if(total_size_>0) memset(data_, 0, sizeof(cplx)*total_size_);
}

} // namespace cntr

#endif  // CNTR_CCLC_TSTP_IMPL
