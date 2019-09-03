#ifndef __HDF5_INTERFACE_RPA_H_
#define __HDF5_INTERFACE_RPA_H_ 



#include "cntr.hpp"
#include "hdf5_interface.hpp"


// ********************************************************************

// Write all the K dependent observables and correlators

// ********************************************************************
template<typename T>
void write_rpa_to_hdf5_file(std::string filename, T & rpa,int kt,int print_k) {
  


  hid_t file_id = open_hdf5_file(filename);
  hid_t group_id, groupG_id,groupchi_id,groupS_id;

  // -- Impurty parameters
  if(rpa.tid_==rpa.tid_root_){
    group_id = create_group(file_id, "parm");
    store_int_attribute_to_hid(group_id, "nt", rpa.nt_); 
    store_int_attribute_to_hid(group_id, "ntau", rpa.ntau_);
  
    store_double_attribute_to_hid(group_id, "beta", rpa.beta_);
    store_double_attribute_to_hid(group_id, "h", rpa.h_);
    store_double_attribute_to_hid(group_id, "mu", rpa.latt_.mu_);

    store_real_data_to_hid(group_id, "U", rpa.latt_.U_.data(), rpa.latt_.U_.size());  
    store_real_data_to_hid(group_id, "delta", rpa.latt_.delta_.data(), rpa.latt_.delta_.size());  
    store_cplx_data_to_hid(group_id, "A", rpa.latt_.A_.data(), rpa.latt_.A_.size());  
  
    close_group(group_id); // End parameters
  }
  // -- Observables-writen
  {
    int nt = rpa.nt_;
    int nk = rpa.nk_;
    hsize_t rank=3;
    typedef boost::multi_array<double, 3> array_type;
    array_type rk(boost::extents[nk][4][nt+2]),sh(boost::extents[nk][4][nt+2]),sf(boost::extents[nk][4][nt+2]);
    array_type hkeff(boost::extents[nk][4][nt+2]),V(boost::extents[nk][4][nt+2]);

    std::vector<hsize_t> size;size.resize(rank);size[0]=nk;size[1]=4;size[2]=nt+2;
    
    for(int k=0;k<rpa.nk_;k++){
      int k1=rpa.latt_.representative_kk(k);
      Eigen::MatrixXcd tmp;
      for(int i=0;i<4;i++){
        for(int tstp=-1;tstp<=nt;tstp++){
          rpa.gather_kk_observables(tstp,kt);
          rk[k][i][tstp+1]=rpa.rk_[k](i,i).real();
          sh[k][i][tstp+1]=rpa.sh_[k](i,i).real();
          sf[k][i][tstp+1]=rpa.sf_[k](i,i).real();

          rpa.latt_.hk(tmp,tstp,rpa.latt_.kpoints_rbz_[k]);
          tmp += rpa.sf_[k]+rpa.sh_[k];
          hkeff[k][i][tstp+1]=tmp(i,i).real();

          rpa.vertex_[k1].get_value(tstp,tmp);
          V[k][i][tstp+1]=tmp(i,i).real();
        }
      }
    }
    group_id = create_group(file_id, "obsk");
    store_real_data_multi_to_hid(group_id, "rk", rk.data(),size,rank);
    store_real_data_multi_to_hid(group_id, "sh", sh.data(),size,rank);
    store_real_data_multi_to_hid(group_id, "sf", sh.data(),size,rank);
    store_real_data_multi_to_hid(group_id, "hkeff", hkeff.data(),size,rank);
    store_real_data_multi_to_hid(group_id, "V", V.data(),size,rank);
    close_group(group_id);
  }

  // -- Single particle Green's functions
  if(rpa.tid_==rpa.tid_root_){
    group_id = create_group(file_id, "Gloc");
    store_herm_greens_function(group_id, rpa.Gloc_);
    close_group(group_id);
  }

  if(print_k==1){
    group_id = create_group(file_id, "Gk");
    store_int_attribute_to_hid(group_id, std::string("nk"),rpa.nk_);

    for(int k=0;k<rpa.nk_;k++){
      if(rpa.tid_map_[k]==rpa.tid_){
        std::ostringstream nameG,namechi,nameS;
        
        nameG << "Gk" << k;
        groupG_id=create_group(group_id,nameG.str().c_str());

        store_herm_greens_function(groupG_id, rpa.kk_functions_[k].G_);
        close_group(groupG_id);

        namechi << "chik" << k;
        groupchi_id=create_group(group_id,namechi.str().c_str());
        store_herm_greens_function(groupchi_id, rpa.kk_functions_[k].chi_);
        close_group(groupchi_id);

        nameS << "sigmak" << k;
        groupS_id=create_group(group_id,nameS.str().c_str());
        store_herm_greens_function(groupS_id, rpa.kk_functions_[k].Sigma_);
        close_group(groupS_id);

      }
    }
    close_group(group_id);
  }

  close_hdf5_file(file_id);
}


template<typename T>
void write_rpa_slice_to_hdf5_file(std::string filename, T & rpa,int kt,int print_k,int outputfrequency) {

  hid_t file_id = open_hdf5_file(filename);
  hid_t group_id, groupG_id,groupchi_id,groupS_id;

  // -- Impurty parameters
  if(rpa.tid_==rpa.tid_root_){
    group_id = create_group(file_id, "parm");
    store_int_attribute_to_hid(group_id, "nt", rpa.nt_); 
    store_int_attribute_to_hid(group_id, "ntau", rpa.ntau_);
  
    store_double_attribute_to_hid(group_id, "beta", rpa.beta_);
    store_double_attribute_to_hid(group_id, "h", rpa.h_);
    store_double_attribute_to_hid(group_id, "mu", rpa.latt_.mu_);

    store_real_data_to_hid(group_id, "U", rpa.latt_.U_.data(), rpa.latt_.U_.size());  
    store_real_data_to_hid(group_id, "delta", rpa.latt_.delta_.data(), rpa.latt_.delta_.size());  
    store_cplx_data_to_hid(group_id, "A", rpa.latt_.A_.data(), rpa.latt_.A_.size());  
  
    close_group(group_id); // End parameters
  }
  // -- Observables
  {
    int nt = rpa.nt_;
    int nk = rpa.nk_;
    hsize_t rank=3;
    typedef boost::multi_array<double, 3> array_type;
    array_type rk(boost::extents[nk][4][nt+2]),sh(boost::extents[nk][4][nt+2]),sf(boost::extents[nk][4][nt+2]);
    array_type hkeff(boost::extents[nk][4][nt+2]),V(boost::extents[nk][4][nt+2]);

    std::vector<hsize_t> size;size.resize(rank);size[0]=nk;size[1]=4;size[2]=nt+2;
    
    for(int k=0;k<rpa.nk_;k++){
      int k1=rpa.latt_.representative_kk(k);
      Eigen::MatrixXcd tmp;
      for(int i=0;i<4;i++){
        for(int tstp=-1;tstp<=nt;tstp++){
          rpa.gather_kk_observables(tstp,kt);
          rk[k][i][tstp+1]=rpa.rk_[k](i,i).real();
          sh[k][i][tstp+1]=rpa.sh_[k](i,i).real();
          sf[k][i][tstp+1]=rpa.sf_[k](i,i).real();

          rpa.latt_.hk(tmp,tstp,rpa.latt_.kpoints_rbz_[k]);
          tmp += rpa.sf_[k]+rpa.sh_[k];
          hkeff[k][i][tstp+1]=tmp(i,i).real();

          rpa.vertex_[k1].get_value(tstp,tmp);
          V[k][i][tstp+1]=tmp(i,i).real();
        }
      }
    }
    group_id = create_group(file_id, "obsk");
    store_real_data_multi_to_hid(group_id, "rk", rk.data(),size,rank);
    store_real_data_multi_to_hid(group_id, "sh", sh.data(),size,rank);
    store_real_data_multi_to_hid(group_id, "sf", sh.data(),size,rank);
    store_real_data_multi_to_hid(group_id, "hkeff", hkeff.data(),size,rank);
    store_real_data_multi_to_hid(group_id, "V", V.data(),size,rank);
    close_group(group_id);
  }

  // -- Vectors of time-slices
  {
    int nt = rpa.nt_;
    int ntau = rpa.ntau_;
    int size=rpa.Gloc_.size1();
    if(rpa.tid_==rpa.tid_root_){
      group_id = create_group(file_id, "Gloc");
      for(int tstp=-1;tstp<=nt;tstp++){
        if(tstp==-1 || tstp%outputfrequency==0){
          Eigen::MatrixXcd tmp,tmp1;
          hid_t group_id_tstp;
          std::ostringstream name_tstp;
          cntr::herm_matrix_timestep<double> Gtmp(tstp,ntau,size);
          rpa.Gloc_.get_timestep(tstp,Gtmp);
          rpa.Gloc_.get_mat(0,tmp);
          Gtmp.get_mat(0,tmp1);
          std::cout << "tmp " << tstp << " " << tmp << " " << std::endl;
          std::cout <<  tmp1 << std::endl;
          name_tstp << tstp;
          group_id_tstp=create_group(group_id,name_tstp.str().c_str());
          store_herm_timestep_greens_function(group_id_tstp,Gtmp);
        }
      }
      close_group(group_id);
    }
  }

  if(print_k==1){
    int nt = rpa.nt_;
    int ntau = rpa.ntau_;
    int size=rpa.Gloc_.size1();
    group_id = create_group(file_id, "Gk");
    store_int_attribute_to_hid(group_id, std::string("nk"),rpa.nk_);

    for(int k=0;k<rpa.nk_;k++){
      if(rpa.tid_map_[k]==rpa.tid_){
        std::ostringstream nameGk,nameChik,nameSk;
        hid_t groupGk_id,groupChik_id,groupSk_id;
        nameGk << "Gk"<<k;
        nameChik << "Chik"<<k;
        nameSk << "Sk"<<k;

        groupGk_id=create_group(group_id,nameGk.str().c_str());
        groupChik_id=create_group(group_id,nameChik.str().c_str());
        groupSk_id=create_group(group_id,nameSk.str().c_str());


        for(int tstp=-1;tstp<=nt;tstp++){
          if(tstp==-1 || tstp%outputfrequency==0){
            std::cout << "tstp k " << tstp << " " << k << " " << rpa.tid_ << std::endl;
            hid_t groupG_tstp_id,groupChi_tstp_id,groupS_tstp_id;
            std::ostringstream nameG,nameChi,nameS;

            // Gk
            cntr::herm_matrix_timestep<double> Gtmp(tstp,ntau,size);
            rpa.kk_functions_[k].G_.get_timestep(tstp,Gtmp);
            nameG << "t" << tstp;
            groupG_tstp_id=create_group(groupGk_id,nameG.str().c_str());
            store_herm_timestep_greens_function(groupG_tstp_id,Gtmp);
            close_group(groupG_tstp_id);

            cntr::herm_matrix_timestep<double> Chitmp(tstp,ntau,size);
            rpa.kk_functions_[k].chi_.get_timestep(tstp,Chitmp);
            nameChi << "t" << tstp;
            groupChi_tstp_id=create_group(groupChik_id,nameChi.str().c_str());
            store_herm_timestep_greens_function(groupChi_tstp_id,Chitmp);
            close_group(groupChi_tstp_id);

            cntr::herm_matrix_timestep<double> Sigmatmp(tstp,ntau,size);
            rpa.kk_functions_[k].Sigma_.get_timestep(tstp,Sigmatmp);
            nameS << "t" << tstp;
            groupS_tstp_id=create_group(groupSk_id,nameS.str().c_str());
            store_herm_timestep_greens_function(groupG_tstp_id,Sigmatmp);
            close_group(groupS_tstp_id);
          }
        }
        close_group(groupGk_id);
        close_group(groupChik_id);
        close_group(groupSk_id);
      }
    }
    close_group(group_id);
  }

  close_hdf5_file(file_id);
}
#endif
// ********************************************************************
