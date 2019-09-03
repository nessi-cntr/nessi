#ifndef __HDF5_INTERFACE_H_
#define __HDF5_INTERFACE_H_ 

// ********************************************************************

// Minimal interface for hdf5 for storing vectors and Green's functions
// Author: Hugo Strand, hugo.strand@gmail.com (2015)

// ******************************************************************** 
#include <string>
#include <vector>
#include <assert.h>
#include <hdf5.h>
#include <complex>
#include <iostream>


#include <eigen3/Eigen/Dense>

// ********************************************************************

hid_t open_hdf5_file(std::string filename);
hid_t read_hdf5_file(std::string filename);
void close_hdf5_file(hid_t file_id);

hid_t create_group(hid_t file_id, std::string label);
hid_t open_group(hid_t file_id, std::string label);
void close_group(hid_t group_id);

void store_array_to_hid(hid_t file_id, std::string label, 
  void * data_ptr, hsize_t * shape, hsize_t len_shape, hid_t type_id);

void store_data_to_hid(hid_t file_id, std::string label, 
  void * data_ptr, size_t data_size, hid_t typie_id);

void store_attribute_to_hid(hid_t file_id, std::string label, 
  void * data_ptr, size_t data_size, hid_t type_id);

void store_cplx_array_to_hid(
  hid_t file_id, std::string label, 
  std::complex<double> * data_ptr, hsize_t * shape, hsize_t len_shape);
  
void store_cplx_data_to_hid(hid_t file_id, std::string label, 
  std::complex<double> * data_ptr, size_t data_size);

void store_real_data_to_hid(hid_t file_id, std::string label, 
  double * data_ptr, size_t data_size);

void store_real_data_multi_to_hid(hid_t file_id, std::string label, 
  double * data_ptr, std::vector<hsize_t> &data_size,hsize_t rank);

void store_int_attribute_to_hid(
  hid_t file_id, std::string label, int data);
void store_double_attribute_to_hid(
  hid_t file_id, std::string label, double data);

void read_data_to_buff(hid_t group_id, std::string name, hsize_t buff_size, void * buff);

// ********************************************************************
template<typename T>
T read_primitive_type(hid_t group_id, std::string name) {
  T buff;
  read_data_to_buff(group_id, name, sizeof(buff), (void *) &buff);
  return buff;
}

// ********************************************************************
template<typename T>
void read_primitive_type_array(hid_t group_id, std::string name,
			       hsize_t size, T * buff) {

  read_data_to_buff(group_id, name, size*sizeof(T), (void *) buff);
}

// ********************************************************************
template<typename T>
void store_matrix_type(hid_t group_id, T & mat, std::string name) {

  hsize_t len_shape = 2, shape[2] = { static_cast<hsize_t>(mat.rows()),
				      static_cast<hsize_t>(mat.cols()) };
  store_cplx_array_to_hid(group_id, name, mat.data(), shape, len_shape);
}

// ********************************************************************
template<typename T>
void store_operator_type(hid_t group_id, T & op, std::string name) {
  typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> matrix_type;
  int size = op.total_size();
  matrix_type mat(size, size);
  mat = op.to_dense();
  store_matrix_type(group_id, mat, name);
}

// ********************************************************************
template<typename T>
void store_herm_greens_function(hid_t group_id, T & g) {

  // Add attributes to file
  store_int_attribute_to_hid(group_id, std::string("ntau"), g.ntau());
  store_int_attribute_to_hid(group_id, std::string("nt"), g.nt());
  store_int_attribute_to_hid(group_id, std::string("sig"), g.sig());

  store_int_attribute_to_hid(group_id, std::string("size1"), g.size1());
  store_int_attribute_to_hid(group_id, std::string("size2"), g.size2());
  store_int_attribute_to_hid(group_id, 
    std::string("element_size"), g.element_size());

  // Keldysh components
  hsize_t len_shape=3, shape[3];
  shape[1] = g.size1(); shape[2] = g.size2();

  shape[0] = g.ntau()+1;
  store_cplx_array_to_hid(group_id, std::string("mat"),
			  g.matptr(0), shape, len_shape);

  shape[0] = (g.nt()+1)*(g.nt()+2)/2;
  store_cplx_array_to_hid(group_id, std::string("ret"),
			  g.retptr(0,0), shape, len_shape);

  store_cplx_array_to_hid(group_id, std::string("les"),
			  g.lesptr(0,0), shape, len_shape);

  shape[0] = (g.nt()+1)*(g.ntau()+1);
  store_cplx_array_to_hid(group_id, std::string("tv"),
			  g.tvptr(0,0), shape, len_shape);

}

// ********************************************************************
template<typename T>
void store_herm_timestep_greens_function(hid_t group_id, T & g) {
  int tstp=g.tstp_;
  // Add attributes to file
  store_int_attribute_to_hid(group_id, std::string("ntau"), g.ntau_);
  store_int_attribute_to_hid(group_id, std::string("tstp"), tstp);
  store_int_attribute_to_hid(group_id, std::string("sig"), g.sig_);

  store_int_attribute_to_hid(group_id, std::string("size1"), g.size1_);
  store_int_attribute_to_hid(group_id, std::string("size2"), g.size2_);
  store_int_attribute_to_hid(group_id, 
  std::string("element_size"), g.element_size_);

  // Keldysh components
  hsize_t len_shape=3, shape[3];
  shape[1] = g.size1_; shape[2] = g.size2_;

  shape[0] = g.ntau()+1;
  store_cplx_array_to_hid(group_id, std::string("mat"),
        g.matptr(0), shape, len_shape);

  shape[0] = tstp+1;
  store_cplx_array_to_hid(group_id, std::string("ret"),
        g.retptr(0), shape, len_shape);

  shape[0] = tstp+1;
  store_cplx_array_to_hid(group_id, std::string("les"),
        g.lesptr(0), shape, len_shape);

  shape[0] = g.ntau()+1;
  store_cplx_array_to_hid(group_id, std::string("tv"),
        g.tvptr(0), shape, len_shape);

}

// ********************************************************************
template<typename T>
T read_herm_matrix(hid_t group_id) {

  // -- Read dimensions

  int nt = read_primitive_type<int>(group_id, "nt");
  int ntau = read_primitive_type<int>(group_id, "ntau");
  int sig = read_primitive_type<int>(group_id, "sig");
  int size1 = read_primitive_type<int>(group_id, "size1");
  int size2 = read_primitive_type<int>(group_id, "size2");
  int element_size = read_primitive_type<int>(group_id, "element_size");

  std::cout << "nt = " << nt << std::endl;
  std::cout << "ntau = " << ntau << std::endl;
  std::cout << "size1 = " << size1 << std::endl;
  std::cout << "size2 = " << size2 << std::endl;
  std::cout << "element_size = " << element_size << std::endl;

  // -- Create instance
  
  T G(nt, ntau, size1, sig);

  // -- Read components

  typedef std::complex<double> complex;
  
  hsize_t mat_size = (ntau+1)*element_size;
  hsize_t ret_size = (nt+1)*(nt+2)/2*element_size;
  hsize_t les_size = (nt+1)*(nt+2)/2*element_size;
  hsize_t tv_size = (nt+1)*(ntau+1)*element_size;

  read_primitive_type_array(group_id, "mat", mat_size, G.matptr(0));
  read_primitive_type_array(group_id, "ret", ret_size, G.retptr(0,0));
  read_primitive_type_array(group_id, "les", les_size, G.lesptr(0,0));
  read_primitive_type_array(group_id, "tv", tv_size, G.tvptr(0,0));
    
  return G;

}

// ********************************************************************
template<typename T>
void store_full_greens_function(
  hid_t group_id, T & g) {

  // Add attributes to file
  store_int_attribute_to_hid(group_id, std::string("ntau"), g.ntau());
  store_int_attribute_to_hid(group_id, std::string("nt"), g.nt());
  store_int_attribute_to_hid(group_id, std::string("sig"), g.sig());

  store_int_attribute_to_hid(group_id, std::string("size1"), g.size1());
  store_int_attribute_to_hid(group_id, std::string("size2"), g.size2());
  store_int_attribute_to_hid(group_id, 
    std::string("element_size"), g.element_size());

  // Keldysh components
  hsize_t len_shape=3, shape[3];
  shape[1] = g.size1(); shape[2] = g.size2();

  shape[0] = g.ntau()+1;
  store_cplx_array_to_hid(group_id, std::string("mat"),
			  g.matptr(0), shape, len_shape);

  shape[0] = (g.nt()+1)*(g.nt()+2)/2;
  store_cplx_array_to_hid(group_id, std::string("les"),
			  g.lesptr(0,0), shape, len_shape);

  store_cplx_array_to_hid(group_id, std::string("gtr"),
			  g.gtrptr(0,0), shape, len_shape);

  shape[0] = (g.nt()+1)*(g.ntau()+1);
  store_cplx_array_to_hid(group_id, std::string("tv"),
			  g.tvptr(0,0), shape, len_shape);

  store_cplx_array_to_hid(group_id, std::string("vt"),
			  g.vtptr(0,0), shape, len_shape);

}

// ********************************************************************
template<typename T>
void write_gf_to_hdf5_file(std::string filename,
			   std::string groupname,
			   T & G) {

  hid_t file_id = open_hdf5_file(filename);
  hid_t group_id;

  // Single particle Green's function
  group_id = create_group(file_id, groupname);
  store_herm_greens_function(group_id, G);
  close_group(group_id);

  close_hdf5_file(file_id);
}

// ********************************************************************
template<typename T>
T read_gf_from_hdf5_file(std::string filename,
			 std::string groupname) {

  hid_t file_id = read_hdf5_file(filename);
  hid_t group_id = open_group(file_id, groupname);

  T G = read_herm_matrix<T>(group_id);
  
  close_group(group_id);
  close_hdf5_file(file_id);

  return G;
}

// ********************************************************************
template<typename T>
T read_function_real(hid_t group_id,int nt, int size1,int size2,std::string data) {
  // -- Read dimensions
  
  std::cout << "nt = " << nt << std::endl;
  std::cout << "size1 = " << size1 << std::endl;
  std::cout << "size2 = " << size2 << std::endl;

  // -- Create instance
  std::vector<double> buffer;
  
  T function(nt, size1);
  int element_size = function.element_size_;
  std::cout << "element_size = " << element_size << std::endl;
  // -- Read components
  hsize_t totsize = (nt+2)*element_size;
  buffer.resize(totsize);

  read_primitive_type_array<double>(group_id, data.c_str(), totsize, buffer.data());
  for(int i=-1;i<=nt;i++){ 
    Eigen::MatrixXcd M(size1,size2);
    if(size1 >1 || size2>1){
      std::cout << "Only implemented for size 1. Aborting" << std::endl;
      abort();
     }

    M(0,0)=std::complex<double>(buffer[i+1],0.0);
    function.set_value(i,M);
  }

  return function;

}

template<typename T>
T read_real_function_from_hdf5_file(std::string filename,std::string groupname,std::string dataname,int nt,int size1,int size2){
  hid_t file_id = read_hdf5_file(filename);
  hid_t group_id = open_group(file_id,groupname);
  T G = read_function_real<T>(group_id,nt,size1,size2,dataname);
  close_group(group_id);
  close_hdf5_file(file_id);
  return G;
}

// ********************************************************************

// Example impurity instance writer, *L1Sz case here)

// ********************************************************************
template<typename T>
void write_imp_to_hdf5_file(std::string filename, T & imp) {
  
  hid_t file_id = open_hdf5_file(filename);
  hid_t group_id, sub_group_id;

  // -- Impurty parameters
  
  group_id = create_group(file_id, "parm");
  store_int_attribute_to_hid(group_id, "nt", imp.nt_); 
  store_int_attribute_to_hid(group_id, "ntau", imp.ntau_);
  
  store_double_attribute_to_hid(group_id, "beta", imp.beta_);
  store_double_attribute_to_hid(group_id, "h", imp.h_);
  store_double_attribute_to_hid(group_id, "mu", imp.mu_);

  store_real_data_to_hid(group_id, "U", imp.U_.data(), imp.U_.size());  
  store_real_data_to_hid(group_id, "epsup", imp.epsup_.data(), imp.epsup_.size());  
  store_real_data_to_hid(group_id, "epsdo", imp.epsdo_.data(), imp.epsdo_.size());  
  
  close_group(group_id); // End parameters
  
  // -- Impurity observables
  
  std::vector<double> vtau(imp.ntau_+2);
  for(int tstp=0; tstp <= imp.ntau_; tstp++) {
    vtau[tstp] = tstp * imp.beta_ / imp.ntau_;
  }

  int nt = imp.nt_;
  std::vector<double> vt(nt+2), vdocc(nt+2), vnup(nt+2), vndo(nt+2);
  //std::vector<double> vt(nt+1), vdocc(nt+1), vnup(nt+1), vndo(nt+1);

  for(int tstp=-1;tstp<=nt;tstp++){
    double docc, nup, ndo;
    imp.local_obs(tstp,nup,ndo,docc);
    vt[tstp+1] = tstp * imp.h_;
    vnup[tstp+1]= nup;
    vndo[tstp+1] = ndo;
    vdocc[tstp+1]= docc;
  }

  group_id = create_group(file_id, "obs");
  store_real_data_to_hid(group_id, "t", vt.data(), vt.size());  
  store_real_data_to_hid(group_id, "tau", vtau.data(), vtau.size());  
  store_real_data_to_hid(group_id, "nup", vnup.data(), vnup.size());  
  store_real_data_to_hid(group_id, "ndo", vndo.data(), vndo.size());  
  store_real_data_to_hid(group_id, "docc", vdocc.data(), vdocc.size());  
  close_group(group_id); // End observables

  // -- Single particle Green's functions
  
  group_id = create_group(file_id, "Gloc_up_");
  store_herm_greens_function(group_id, imp.Gloc_up_);
  close_group(group_id);

  group_id = create_group(file_id, "Gloc_do_");
  store_herm_greens_function(group_id, imp.Gloc_do_);
  close_group(group_id);

  close_hdf5_file(file_id);
}


template<typename T>
void write_lattice_to_hdf5_file_gw(std::string filename,T & selfconsistency,int dt,int nt,int taskid,int mpi_imp){
  std::ostringstream name;
  if(taskid == mpi_imp){
    // Local component
    name.str("");
    name << filename << "_local";
    hid_t file_id = open_hdf5_file(name.str().c_str());
    hid_t sub_group_id = create_group(file_id, "parm");
    store_int_attribute_to_hid(sub_group_id, "nt", selfconsistency.nt_); 
    store_int_attribute_to_hid(sub_group_id, "ntau", selfconsistency.ntau_);
    store_double_attribute_to_hid(sub_group_id, "beta", selfconsistency.beta_);
    store_double_attribute_to_hid(sub_group_id, "h", selfconsistency.h_);
    close_group(sub_group_id);
    hid_t group_id = create_group(file_id, "Gloc_");
    store_herm_greens_function(group_id, selfconsistency.Gloc_);
    close_group(group_id);
    int ntau=selfconsistency.ntau_;
    double beta=selfconsistency.beta_;
    std::vector<double> vtau(ntau+2);
    for(int tstp=0; tstp <= ntau; tstp++) {
      vtau[tstp] = tstp * beta / ntau;
    }
    group_id = create_group(file_id, "obs");
    // store_real_data_to_hid(group_id, "t", vt.data(), vt.size());  
    // store_real_data_to_hid(group_id, "tau", vtau.data(), vtau.size());
    close_group(group_id);
    close_hdf5_file(file_id);
  }
  for(int q=0;q<selfconsistency.nklocal_;q++){

      int kfull=selfconsistency.kindex_local_[q];
      name.str("");
      name << filename << "_k" << kfull;
      hid_t file_id = open_hdf5_file(name.str().c_str());
      hid_t sub_group_id = create_group(file_id, "parm");
      store_int_attribute_to_hid(sub_group_id, "nt", selfconsistency.nt_); 
      store_int_attribute_to_hid(sub_group_id, "ntau", selfconsistency.ntau_);
      store_double_attribute_to_hid(sub_group_id, "beta", selfconsistency.beta_);
      store_double_attribute_to_hid(sub_group_id, "kweight", selfconsistency.kweight_local_[q]);
      store_double_attribute_to_hid(sub_group_id, "h", selfconsistency.h_);
      store_real_data_to_hid(sub_group_id, "k", selfconsistency.kpoints_local_[q].data(),  selfconsistency.kpoints_local_[q].size());

      std::vector<double> epsk(nt+2),fockk(nt+2),epseff(nt+2);
      for(int tstp=-1;tstp<=nt;tstp++){
        Eigen::MatrixXcd tmpeps;
        selfconsistency.epsk_[q].get_value(tstp,tmpeps);
        epsk[tstp+1]=std::real(tmpeps(0,0));

        selfconsistency.fockk_[q].get_value(tstp,tmpeps);
        fockk[tstp+1]=std::real(tmpeps(0,0));
        
        // selfconsistency.epseff_[q].get_value(tstp,tmpeps);
        // epseff[tstp+1]=std::real(tmpeps(0,0));
      }
      store_real_data_to_hid(sub_group_id, "epsk", epsk.data(),epsk.size());
      store_real_data_to_hid(sub_group_id, "fockk_", fockk.data(),fockk.size());
      // store_real_data_to_hid(sub_group_id, "epseff", epseff.data(),epseff.size());
      close_group(sub_group_id); // End parameters


  //     // -- Green's functions
      selfconsistency.Gk_[q].write_to_hdf5_slices(file_id,"G",dt);
      selfconsistency.Wk_[q].write_to_hdf5_slices(file_id,"W",dt);
      // selfconsistency.Sigmak_nl_[q].write_to_hdf5_slices(file_id,"Sigma",dt);
      // selfconsistency.Pik_nl_[q].write_to_hdf5_slices(file_id,"Pi",dt);
      close_hdf5_file(file_id);
  }

}

template<typename T>
void write_lattice_to_hdf5_file_gw_afm(std::string filename,T & selfconsistency,int dt,int nt,int taskid,int mpi_imp){
  std::ostringstream name;
  if(taskid == mpi_imp){
    // Local component
    name.str("");
    name << filename << "_local";
    hid_t file_id = open_hdf5_file(name.str().c_str());
    hid_t sub_group_id = create_group(file_id, "parm");
    store_int_attribute_to_hid(sub_group_id, "nt", selfconsistency.nt_); 
    store_int_attribute_to_hid(sub_group_id, "ntau", selfconsistency.ntau_);
    store_double_attribute_to_hid(sub_group_id, "beta", selfconsistency.beta_);
    store_double_attribute_to_hid(sub_group_id, "h", selfconsistency.h_);
    close_group(sub_group_id);
    
    hid_t group_id = create_group(file_id, "Gloc");
    store_herm_greens_function(group_id, selfconsistency.Gloc_);
    close_group(group_id);
    
    int ntau=selfconsistency.ntau_;
    double beta=selfconsistency.beta_;
    std::vector<double> vtau(ntau+2);
    for(int tstp=0; tstp <= ntau; tstp++) {
      vtau[tstp] = tstp * beta / ntau;
    }
    group_id = create_group(file_id, "obs");
    // store_real_data_to_hid(group_id, "t", vt.data(), vt.size());  
    // store_real_data_to_hid(group_id, "tau", vtau.data(), vtau.size());
    close_group(group_id);
    close_hdf5_file(file_id);
  }
  for(int q=0;q<selfconsistency.nklocal_;q++){

      int kfull=selfconsistency.kindex_local_[q];
      name.str("");
      name << filename << "_k" << kfull;
      hid_t file_id = open_hdf5_file(name.str().c_str());
      hid_t sub_group_id = create_group(file_id, "parm");
      store_int_attribute_to_hid(sub_group_id, "nt", selfconsistency.nt_); 
      store_int_attribute_to_hid(sub_group_id, "ntau", selfconsistency.ntau_);
      store_double_attribute_to_hid(sub_group_id, "beta", selfconsistency.beta_);
      store_double_attribute_to_hid(sub_group_id, "kweight", selfconsistency.kweight_local_[q]);
      store_double_attribute_to_hid(sub_group_id, "h", selfconsistency.h_);
      store_real_data_to_hid(sub_group_id, "k", selfconsistency.kpoints_local_[q].data(),  selfconsistency.kpoints_local_[q].size());

      selfconsistency.epsk_[q].write_to_hdf5(file_id,"eps");
      selfconsistency.hartree_.write_to_hdf5(file_id,"hartree");
      selfconsistency.fockk_[q].write_to_hdf5(file_id,"fock");
      selfconsistency.epseff_[q].write_to_hdf5(file_id,"eta");
      close_group(sub_group_id); // End parameters


  //     // -- Green's functions
      selfconsistency.Gk_[q].write_to_hdf5_slices(file_id,"G",dt);
      selfconsistency.Wk_[q].write_to_hdf5_slices(file_id,"W",dt);
      close_hdf5_file(file_id);
  }

}


template<typename T>
void write_lattice_to_hdf5_file_dp(std::string filename,T & selfconsistency,int dt,int nt,int taskid,int mpi_imp){
  std::ostringstream name;
  if(taskid == mpi_imp){
    // Local component
    name.str("");
    name << filename << "_local";
    hid_t file_id = open_hdf5_file(name.str().c_str());
    hid_t sub_group_id = create_group(file_id, "parm");
    store_int_attribute_to_hid(sub_group_id, "nt", selfconsistency.nt_); 
    store_int_attribute_to_hid(sub_group_id, "ntau", selfconsistency.ntau_);
    store_double_attribute_to_hid(sub_group_id, "beta", selfconsistency.beta_);
    store_double_attribute_to_hid(sub_group_id, "h", selfconsistency.h_);
    close_group(sub_group_id);
    hid_t group_id = create_group(file_id, "Gloc_");
    store_herm_greens_function(group_id, selfconsistency.Gloc_);
    close_group(group_id);
    int ntau=selfconsistency.ntau_;
    double beta=selfconsistency.beta_;
    std::vector<double> vtau(ntau+2);
    for(int tstp=0; tstp <= ntau; tstp++) {
      vtau[tstp] = tstp * beta / ntau;
    }
    group_id = create_group(file_id, "obs");
    // store_real_data_to_hid(group_id, "t", vt.data(), vt.size());  
    // store_real_data_to_hid(group_id, "tau", vtau.data(), vtau.size());
    close_group(group_id);
    close_hdf5_file(file_id);
  }
  for(int q=0;q<selfconsistency.nklocal_;q++){

      int kfull=selfconsistency.kindex_local_[q];
      name.str("");
      name << filename << "_k" << kfull;
      hid_t file_id = open_hdf5_file(name.str().c_str());
      hid_t sub_group_id = create_group(file_id, "parm");
      store_int_attribute_to_hid(sub_group_id, "nt", selfconsistency.nt_); 
      store_int_attribute_to_hid(sub_group_id, "ntau", selfconsistency.ntau_);
      store_double_attribute_to_hid(sub_group_id, "beta", selfconsistency.beta_);
      store_double_attribute_to_hid(sub_group_id, "kweight", selfconsistency.kweight_local_[q]);
      store_double_attribute_to_hid(sub_group_id, "h", selfconsistency.h_);
      store_real_data_to_hid(sub_group_id, "k", selfconsistency.kpoints_local_[q].data(),  selfconsistency.kpoints_local_[q].size());

      close_group(sub_group_id); // End parameters
      selfconsistency.epsk_[q].write_to_hdf5(file_id,"eps");
      selfconsistency.hartree_.write_to_hdf5(file_id,"hartree");
      selfconsistency.fockk_[q].write_to_hdf5(file_id,"fock");
      selfconsistency.seff_.write_to_hdf5(file_id,"seff");
      selfconsistency.eta_[q].write_to_hdf5(file_id,"eta");
      selfconsistency.Gk_[q].write_to_hdf5_slices(file_id,"G",dt);
      selfconsistency.Wk_[q].write_to_hdf5_slices(file_id,"W",dt);
      selfconsistency.Sigmak_nl_[q].write_to_hdf5_slices(file_id,"Sigma",dt);
      selfconsistency.Pik_nl_[q].write_to_hdf5_slices(file_id,"Pi",dt);
      close_hdf5_file(file_id);
  }
}


template<typename T>
void write_lattice_to_hdf5_file_dpfree(std::string filename,T & selfconsistency,int dt,int nt,int taskid,int mpi_imp){
  std::ostringstream name;
  if(taskid == mpi_imp){
    // Local component
    name.str("");
    name << filename << "_local";
    hid_t file_id = open_hdf5_file(name.str().c_str());
    hid_t sub_group_id = create_group(file_id, "parm");
    store_int_attribute_to_hid(sub_group_id, "nt", selfconsistency.nt_); 
    store_int_attribute_to_hid(sub_group_id, "ntau", selfconsistency.ntau_);
    store_double_attribute_to_hid(sub_group_id, "beta", selfconsistency.beta_);
    store_double_attribute_to_hid(sub_group_id, "h", selfconsistency.h_);
    close_group(sub_group_id);
    hid_t group_id = create_group(file_id, "Gloc_");
    store_herm_greens_function(group_id, selfconsistency.Gloc_);
    close_group(group_id);
    int ntau=selfconsistency.ntau_;
    double beta=selfconsistency.beta_;
    std::vector<double> vtau(ntau+2);
    for(int tstp=0; tstp <= ntau; tstp++) {
      vtau[tstp] = tstp * beta / ntau;
    }
    group_id = create_group(file_id, "obs");
    // store_real_data_to_hid(group_id, "t", vt.data(), vt.size());  
    // store_real_data_to_hid(group_id, "tau", vtau.data(), vtau.size());
    close_group(group_id);
    close_hdf5_file(file_id);
  }
  for(int q=0;q<selfconsistency.nklocal_;q++){

      int kfull=selfconsistency.kindex_local_[q];
      name.str("");
      name << filename << "_k" << kfull;
      hid_t file_id = open_hdf5_file(name.str().c_str());
      hid_t sub_group_id = create_group(file_id, "parm");
      store_int_attribute_to_hid(sub_group_id, "nt", selfconsistency.nt_); 
      store_int_attribute_to_hid(sub_group_id, "ntau", selfconsistency.ntau_);
      store_double_attribute_to_hid(sub_group_id, "beta", selfconsistency.beta_);
      store_double_attribute_to_hid(sub_group_id, "kweight", selfconsistency.kweight_local_[q]);
      store_double_attribute_to_hid(sub_group_id, "h", selfconsistency.h_);
      store_real_data_to_hid(sub_group_id, "k", selfconsistency.kpoints_local_[q].data(),  selfconsistency.kpoints_local_[q].size());

      close_group(sub_group_id); // End parameters
      selfconsistency.epsk_[q].write_to_hdf5(file_id,"eps");
      selfconsistency.hartree_.write_to_hdf5(file_id,"hartree");
      selfconsistency.fockk_[q].write_to_hdf5(file_id,"fock");
      selfconsistency.Gk_[q].write_to_hdf5_slices(file_id,"G",dt);
      selfconsistency.Wk_[q].write_to_hdf5_slices(file_id,"W",dt);
      selfconsistency.Sigmak_nl_[q].write_to_hdf5_slices(file_id,"Sigma",dt);
      selfconsistency.DeltaP_[q].write_to_hdf5_slices(file_id,"DeltaP",dt);
      selfconsistency.Pik_nl_[q].write_to_hdf5_slices(file_id,"Pi",dt);
      close_hdf5_file(file_id);
  }
}



#endif
// ********************************************************************
