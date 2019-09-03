#ifndef __HDF5_INTERFACE_CNTR_H_
#define __HDF5_INTERFACE_CNTR_H_ 

/*//////////////////////////////////////////////////////////////////////////////////////////

hdf5 output routines for green functions, time-slices, cntr::functions etc. 

//////////////////////////////////////////////////////////////////////////////////////////*/

#include "cntr/cntr.hpp"
#include "cntr/hdf5/hdf5_interface.hpp"


/////////////////////////////////////////////////////////////////////////////////
// WRITE G(tav,trel) at equally spaced time slices as 3d data array
// [...]/les/tav/  trel i j Gles(tav+trel,tav-trel)_{ij} 
// [...]/gtr/tav/  trel i j Ggtr(tav+trel,tav-trel)_{ij}
// where tav=0,dt,2*dt,...;  trel=0...min(tav,nt-tav)
template<class GG> 
void write_gf_tavtrel_to_hdf5_group(hid_t group_id,GG &g,int dt) {
	typedef std::complex<double> complex;
	char name[100];
	int nt=g.nt();
	int size1=g.size1();
	int size2=g.size2();
	int element_size=g.element_size();
	// Add attributes to file
	store_int_attribute_to_hid(group_id, std::string("nt"), nt);
	store_int_attribute_to_hid(group_id, std::string("size1"), size1);
	store_int_attribute_to_hid(group_id, std::string("size2"), size2);
	store_int_attribute_to_hid(group_id, std::string("element_size"), g.element_size());
	if(nt<0) return;
	complex *ggtr=new complex [(nt+1)*g.element_size()]; // temp storage
	complex *gles=new complex [(nt+1)*g.element_size()]; // temp storage
	hid_t group_id_les=create_group(group_id,std::string("les"));
	hid_t group_id_gtr=create_group(group_id,std::string("gtr"));
	hsize_t len_shape=3, shape[3];
	shape[1] = size1; 
	shape[2] = size2;
	for(int tav=0;tav<=nt;tav+=dt){
		int len=(tav<=nt-tav ? tav : nt-tav);
		std::sprintf(name,"%d",tav);
		shape[0] = len+1;
		// read data: trel2 is trel/2
		for(int trel2=0;trel2<=len;trel2++){
			cdmatrix ggtrtt,glestt;
			g.get_les(tav+trel2,tav-trel2,glestt);
			g.get_gtr(tav+trel2,tav-trel2,ggtrtt);
			for(int i=0;i<element_size;i++){
				gles[trel2*element_size+i]=glestt(i/size2,i%size2);
				ggtr[trel2*element_size+i]=ggtrtt(i/size2,i%size2);
			}
		}
		store_cplx_array_to_hid(group_id_les,std::string(name),gles,shape,len_shape);
		store_cplx_array_to_hid(group_id_gtr,std::string(name),ggtr,shape,len_shape);
	}
	delete [] ggtr;
	delete [] gles;
}
template<typename T>
void write_gf_tavtrel_to_hdf5_file(hid_t file_id,std::string groupname,T & g,int dt) {
	hid_t group_id;
	group_id = create_group(file_id, groupname);
	write_gf_tavtrel_to_hdf5_group(group_id,g,dt);
	close_group(group_id);
}
template<typename T>
void write_gf_tavtrel_to_hdf5_file(std::string filename,std::string groupname,T & g,int dt) {
  hid_t file_id = open_hdf5_file(filename);
  write_gf_tavtrel_to_hdf5_file(file_id,groupname,g,dt);
  close_hdf5_file(file_id);
}


// READING/WRITING CNTR::FUNCTION
template<typename T>
void write_to_hdf5_group(hid_t group_id, cntr::function<T> &f) {
  int nt=f.nt();
  // Add attributes to file
  store_int_attribute_to_hid(group_id, std::string("nt"), nt);
  store_int_attribute_to_hid(group_id, std::string("size1"), f.size1_);
  store_int_attribute_to_hid(group_id, std::string("size2"), f.size2_);
  store_int_attribute_to_hid(group_id, std::string("element_size"), f.element_size_);
  // write data
  hsize_t len_shape=3, shape[3];
  shape[1] = f.size1_; 
  shape[2] = f.size2_;
  shape[0] = nt+2;
  store_cplx_array_to_hid(group_id, std::string("data"),f.ptr(0), shape, len_shape);
}
template<typename T>
void write_to_hdf5_file(hid_t &file_id,std::string groupname,cntr::function<T> &f) {
  hid_t group_id;
  group_id = create_group(file_id, groupname);
  write_to_hdf5_group(group_id,f);
  close_group(group_id);
}
template<typename T>
void write_to_hdf5_file(std::string filename,std::string groupname,cntr::function<T> &f) {
  hid_t file_id = open_hdf5_file(filename);
  write_to_hdf5_file(file_id,groupname,f);  
  close_hdf5_file(file_id);
}
template<typename T> 
void read_from_hdf5_group(hid_t group_id,cntr::function<T> &f) {
  // -- Read dimensions
  int nt = read_primitive_type<int>(group_id, "nt");
  int size1 = read_primitive_type<int>(group_id, "size1");
  int size2 = read_primitive_type<int>(group_id, "size2");
  int element_size = read_primitive_type<int>(group_id, "element_size");
  std::cout << "nt = " << nt << std::endl;
  std::cout << "size1 = " << size1 << std::endl;
  std::cout << "size2 = " << size2 << std::endl;
  std::cout << "element_size = " << element_size << std::endl;
  // RESIZE F
  f=cntr::function<T>(nt,size1);
  hsize_t data_size = (nt+2)*element_size;
  read_primitive_type_array(group_id, "data", data_size, f.ptr(0));
}
template<typename T> 
void read_from_hdf5_file(std::string filename,std::string groupname,cntr::function<T> &f) {
  hid_t file_id = read_hdf5_file(filename);
  hid_t group_id = open_group(file_id, groupname);
  read_from_hdf5_group(group_id,f);
  close_group(group_id);
  close_hdf5_file(file_id);
}


#endif
// ********************************************************************
