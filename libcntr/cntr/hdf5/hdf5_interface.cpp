// ********************************************************************

// Author: Hugo Strand, hugo.strand@gmail.com (2015)

// ********************************************************************

#include <complex>
#include <iostream>

// ********************************************************************
#include "hdf5_interface.hpp"
// ********************************************************************
hid_t open_hdf5_file(std::string filename) {
  hid_t file_id = H5Fcreate(filename.c_str(), 
    H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  return file_id;
}

// ********************************************************************
hid_t read_hdf5_file(std::string filename) {
  hid_t file_id = H5Fopen(filename.c_str(),
    H5F_ACC_RDONLY, H5P_DEFAULT);
  return file_id; 
}

// ********************************************************************
void close_hdf5_file(hid_t file_id) {
  herr_t status;
  status = H5Fclose(file_id);
}

// ********************************************************************
hid_t create_group(hid_t file_id, std::string label) {
  hid_t group_id = H5Gcreate(file_id, label.c_str(), 
    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  return group_id;
}

// ********************************************************************
hid_t open_group(hid_t file_id, std::string label) {
  hid_t group_id = H5Gopen(file_id, label.c_str(), H5P_DEFAULT);
  return group_id;
}

// ********************************************************************
void close_group(hid_t group_id) {
  // Error Handling?
  herr_t status;
  status = H5Gclose(group_id);
}

// ********************************************************************
void store_array_to_hid(hid_t file_id, std::string label, 
  void * data_ptr, hsize_t * shape, hsize_t len_shape, hid_t type_id) {

  //hsize_t len_shape=1, shape[1]; shape[0] = data_size;
  hid_t space_id = H5Screate_simple(len_shape, shape, NULL);  

  hid_t data_id = H5Dcreate(file_id, label.c_str(),
    type_id, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // Error handling?
  herr_t status;
  status = H5Dwrite(data_id, type_id,
    H5S_ALL, H5S_ALL, H5P_DEFAULT, data_ptr);

  status = H5Dclose(data_id);
  status = H5Sclose(space_id);

}
void store_array_multi_to_hid(hid_t file_id, std::string label, 
  double  *data_ptr, std::vector<hsize_t> &shape, hsize_t len_shape, hid_t type_id) {
  //hsize_t len_shape=1, shape[1]; shape[0] = data_size;
  hid_t space_id = H5Screate_simple(len_shape, shape.data(), NULL);  
  hid_t data_id = H5Dcreate(file_id, label.c_str(),
    type_id, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  // Error handling?
  herr_t status;
  status = H5Dwrite(data_id, type_id,
    H5S_ALL, H5S_ALL, H5P_DEFAULT, data_ptr);
  status = H5Dclose(data_id);
  status = H5Sclose(space_id);

}

// ********************************************************************
void store_data_to_hid(hid_t file_id, std::string label, 
  void * data_ptr, size_t data_size, hid_t type_id) {

  hsize_t len_shape=1, shape[1]; shape[0] = data_size;
  store_array_to_hid(file_id, label, data_ptr, shape, len_shape, type_id);

}

void store_data_multi_to_hid(hid_t file_id, std::string label, 
  double  *data_ptr, std::vector<hsize_t> &data_size,hsize_t rank, hid_t type_id) {
  store_array_multi_to_hid(file_id, label, data_ptr, data_size, rank, type_id);
}

// ********************************************************************
void store_attribute_to_hid(hid_t file_id, std::string label, 
  void * data_ptr, size_t data_size, hid_t type_id) {

  hsize_t len_shape=1, shape[1]; shape[0] = data_size;
  hid_t space_id = H5Screate_simple(len_shape, shape, NULL);  

  hid_t data_id = H5Dcreate(file_id, label.c_str(),
    type_id, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // Error handling?
  herr_t status;
  status = H5Dwrite(data_id, type_id,
    H5S_ALL, H5S_ALL, H5P_DEFAULT, data_ptr);

  status = H5Dclose(data_id);
  status = H5Sclose(space_id);
}

// ********************************************************************
void store_cplx_array_to_hid(
  hid_t file_id, std::string label, 
  std::complex<double> * data_ptr, hsize_t * shape, hsize_t len_shape) {

  // -- COMPLEX COMPOUND HDF5 TYPE
  hid_t H5T_COMPLEX = H5Tcreate(H5T_COMPOUND, sizeof(double)*2);
  H5Tinsert(H5T_COMPLEX, "r", 0, H5T_NATIVE_DOUBLE);
  H5Tinsert(H5T_COMPLEX, "i", sizeof(double), H5T_NATIVE_DOUBLE);
  
  store_array_to_hid(file_id, label, (void *) data_ptr,
		     shape, len_shape, H5T_COMPLEX);
}

// ********************************************************************
void store_cplx_data_to_hid(
  hid_t file_id, std::string label, 
  std::complex<double> * data_ptr, size_t data_size) {

  hsize_t len_shape=1, shape[1]; shape[0] = data_size;
  store_cplx_array_to_hid(file_id, label, data_ptr, shape, len_shape);

}

// ********************************************************************
void store_real_data_to_hid(
  hid_t file_id, std::string label, 
  double * data_ptr, size_t data_size) {
  
  store_data_to_hid(file_id, label, (void *) data_ptr,
		      data_size, H5T_NATIVE_DOUBLE);
}

void store_real_data_multi_to_hid(hid_t file_id, std::string label, 
  double * data_ptr, std::vector<hsize_t> &data_size,hsize_t rank){
  store_data_multi_to_hid(file_id, label,data_ptr,
          data_size,rank,H5T_NATIVE_DOUBLE);
}
// ********************************************************************
void store_int_attribute_to_hid(
  hid_t file_id, std::string label, int data) {

  int * data_ptr;
  size_t data_size=1;
  data_ptr = &data;

  store_attribute_to_hid(file_id, label, (void *) data_ptr, 
			 data_size, H5T_NATIVE_INT);
}

// ********************************************************************
void store_double_attribute_to_hid(
  hid_t file_id, std::string label, double data) {

  double * data_ptr;
  size_t data_size=1;
  data_ptr = &data;

  store_attribute_to_hid(file_id, label, (void *) data_ptr, 
			 data_size, H5T_NATIVE_DOUBLE);
}

// ********************************************************************
void read_data_to_buff(hid_t group_id, std::string name, hsize_t buff_size, void * buff) {

  hid_t data_id = H5Dopen(group_id, name.c_str(), H5P_DEFAULT);
  hid_t type_id = H5Dget_type(data_id);
  hsize_t data_size = H5Dget_storage_size(data_id);

  assert(data_size == buff_size && "Mismatch in data and buffer size");

  H5Dread(data_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, buff);

  H5Tclose(type_id);
  H5Dclose(data_id);
}

// ********************************************************************
