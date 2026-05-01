#ifndef NESS2_HERM_MATRIX_NESS_IMPL_HPP
#define NESS2_HERM_MATRIX_NESS_IMPL_HPP

#include "ness2_herm_matrix_ness_decl.hpp"

namespace ness2 {

  herm_matrix_ness::herm_matrix_ness(int Nft, int size, unsigned FFTW_FLAG)  : Nft_(Nft), size1_(size), ret_(Nft,size,FFTW_FLAG),les_(Nft,size,FFTW_FLAG) {}
  
  herm_matrix_ness::~herm_matrix_ness() {}
  
  herm_matrix_ness::herm_matrix_ness(const herm_matrix_ness& other)
    : Nft_(other.Nft_), size1_(other.size1_), ret_(other.ret_), les_(other.les_) {}

  herm_matrix_ness& herm_matrix_ness::operator=(const herm_matrix_ness& other) {
    if (this == &other) return *this;
    Nft_ = other.Nft_;
    size1_ = other.size1_;
    ret_ = other.ret_;
    les_ = other.les_;
    return *this;
  }
  
  void herm_matrix_ness::fft_to_freq() {
    ret_.fft_to_freq();
    les_.fft_to_freq();
  }
  
  void herm_matrix_ness::fft_to_time() {
    ret_.fft_to_time();
    les_.fft_to_time();
  }
  
  void herm_matrix_ness::incr(const herm_matrix_ness& B, cplx alpha, fft_domain type) {
    ret_.incr(B.ret_, alpha, type);
    les_.incr(B.les_, alpha, type);
  }
  
  void herm_matrix_ness::smul(cplx alpha, fft_domain type) {
    ret_.smul(alpha, type);
    les_.smul(alpha, type);
  }
  
  void herm_matrix_ness::set_matrixelement(int i1, int i2, const herm_matrix_ness& B, int j1, int j2, fft_domain type) {
    ret_.set_matrixelement(i1, i2, B.ret_, j1, j2, type);
    les_.set_matrixelement(i1, i2, B.les_, j1, j2, type);
  }
  
  void herm_matrix_ness::set_zero(fft_domain type) {
    ret_.set_zero(type);
    les_.set_zero(type);
  }
  
  // Explicit instantiations for cdmatrix
  template void herm_matrix_ness::left_multiply<cdmatrix>(const cdmatrix& M, fft_domain type);
  template void herm_matrix_ness::right_multiply<cdmatrix>(const cdmatrix& M, fft_domain type);
  template void herm_matrix_ness::left_multiply_hermconj<cdmatrix>(const cdmatrix& M, fft_domain type);
  template void herm_matrix_ness::right_multiply_hermconj<cdmatrix>(const cdmatrix& M, fft_domain type);
  template void herm_matrix_ness::set_ret<cdmatrix>(int t, const cdmatrix& M, fft_domain type);
  template void herm_matrix_ness::get_ret<cdmatrix>(int t, cdmatrix& M, fft_domain type) const;
  template void herm_matrix_ness::set_les<cdmatrix>(int t, const cdmatrix& M, fft_domain type);
  template void herm_matrix_ness::get_les<cdmatrix>(int t, cdmatrix& M, fft_domain type) const;



  void herm_matrix_ness::print_to_file(const std::string& filename, int precision) const {
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
      throw std::runtime_error("Failed to open output file: " + filename);
    }
    
    outfile << "# " << Nft_ << " " << size1_ << "\n";
    
    auto dump_component = [&](const fft_array& comp, const std::string& label, fft_domain type) {
      for (int t = 0; t < Nft_; ++t) {
	outfile << label << ": " << t;
	
	cdmatrix M;
	comp.get_element(t, M, type);
	for (int i = 0; i < ret_.size1_; ++i) {
	  for (int j = 0; j < ret_.size1_; ++j) {
	    outfile << " " << std::setprecision(precision) << M(i, j).real()
		    << " " << std::setprecision(precision) << M(i, j).imag();
	  }
	}
	outfile << "\n";
      }
    };
    
    dump_component(ret_, "ret_t", fft_domain::time);
    dump_component(les_, "les_t", fft_domain::time);
    dump_component(ret_, "ret_w", fft_domain::freq);
    dump_component(les_, "les_w", fft_domain::freq);
    
    outfile.close();
  }
  
  void herm_matrix_ness::read_from_file(const std::string& filename, unsigned FFTW_FLAG) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
      throw std::runtime_error("Failed to open input file: " + filename);
    }
    
    std::string line;
    int precision_dummy;
    std::getline(infile, line);
    std::istringstream header(line.substr(1)); // skip '#'
    int file_size1,file_nft;
    header >> file_nft >> file_size1;
    
    // Resize if needed
    if (file_nft != Nft_ || file_size1 != size1_) {
      Nft_ = file_nft;
      ret_.resize(Nft_, file_size1, FFTW_FLAG);
      les_.resize(Nft_, file_size1, FFTW_FLAG);
    }
    
    auto load_component = [&](fft_array& comp, const std::string& label, fft_domain type) {
      for (int t = 0; t < Nft_; ++t) {
	if (!std::getline(infile, line)) {
	  throw std::runtime_error("Unexpected end of file when reading " + label);

	}
	std::istringstream linestream(line);
	std::string tag;
	int t_check;
	linestream >> tag >> t_check;
	if (tag != label + ":" || t_check != t) {
	  throw std::runtime_error("Mismatched line tag or time index in " + label);
	}
	
	cdmatrix M(file_size1, file_size1);
	for (int i = 0; i < size1_; ++i) {
	  for (int j = 0; j < size1_; ++j) {
	    double re, im;
	    if (!(linestream >> re >> im)) {
	      throw std::runtime_error("Failed to read complex value in " + label);
	    }
	    M(i, j) = cplx(re, im);
	  }
	}
	comp.set_element(t, M, type);
      }
    };
    
    load_component(ret_, "ret_t", fft_domain::time);
    load_component(les_, "les_t", fft_domain::time);
    load_component(ret_, "ret_w", fft_domain::freq);
    load_component(les_, "les_w", fft_domain::freq);
    
    infile.close();
  }
    
#if CNTR_USE_HDF5 == 1

  void herm_matrix_ness::write_to_hdf5(hid_t group_id) const {
    store_int_attribute_to_hid(group_id, "Nft", Nft_);
    store_int_attribute_to_hid(group_id, "size1", size1_);    
    // Create subgroups for retarded and lesser components
    hid_t ret_group = create_group(group_id, "ret");
    ret_.write_to_hdf5(ret_group);
    close_group(ret_group);
    
    hid_t les_group = create_group(group_id, "les");
    les_.write_to_hdf5(les_group);
    close_group(les_group);
  }
  
  void herm_matrix_ness::read_from_hdf5(hid_t group_id, unsigned FFTW_FLAG) {
    int size1 = read_primitive_type<int>(group_id, "size1");
    int Nft = read_primitive_type<int>(group_id, "Nft");
    
    if ( Nft != Nft_ || size1!=size1_) {
      // Resize both components consistently
      Nft_ = Nft;
      size1_=size1;
      ret_.resize(Nft_, size1, FFTW_FLAG);
      les_.resize(Nft_, size1, FFTW_FLAG);
    }
    
    hid_t ret_group = open_group(group_id, "ret");
    ret_.read_from_hdf5(ret_group, FFTW_FLAG);
    close_group(ret_group);
    
    hid_t les_group = open_group(group_id, "les");
    les_.read_from_hdf5(les_group, FFTW_FLAG);
    close_group(les_group);
  }
  
  void herm_matrix_ness::write_to_hdf5(hid_t group_id, const char* groupname) const {
    hid_t sub_group_id = create_group(group_id, groupname);
    this->write_to_hdf5(sub_group_id);
    close_group(sub_group_id);
  }



  void herm_matrix_ness::write_to_hdf5(const char* filename, const char* groupname) const {
    // 1) Open existing file for read/write; if it doesn't exist, create it
    hid_t file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
    if (file_id < 0) {
      file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      if (file_id < 0) throw std::runtime_error("write_to_hdf5: cannot open/create file");
    }
    this->write_to_hdf5(file_id, groupname);
    close_hdf5_file(file_id);
  }
  
  
  //void herm_matrix_ness::write_to_hdf5(const char* filename, const char* groupname)  const{
  //  hid_t file_id = open_hdf5_file(filename);
  //  this->write_to_hdf5(file_id, groupname);
  //  close_hdf5_file(file_id);
  // }

  void herm_matrix_ness::read_from_hdf5(hid_t group_id, const char* groupname, unsigned FFTW_FLAG) {
    hid_t sub_group_id = open_group(group_id, groupname);
    this->read_from_hdf5(sub_group_id, FFTW_FLAG);
    close_group(sub_group_id);
  }
  
  void herm_matrix_ness::read_from_hdf5(const char* filename, const char* groupname, unsigned FFTW_FLAG) {
    hid_t file_id = read_hdf5_file(filename);
    this->read_from_hdf5(file_id, groupname, FFTW_FLAG);
    close_hdf5_file(file_id);
  }
  
#endif

  void herm_matrix_ness::integral_transform_to_freq(double h,fft_integral_method method) {
    cdmatrix M;
    int nhalf=Nft_/2;
    if (Nft_ < 4) {
      throw std::runtime_error("integral_transform_to_time:cubic correction requires Nft >= 4");
    }        
    grid_info grid(Nft_,h);
    // transformation of the retarded:
    ret_.fft_to_freq();
    if(method==FFT_CUBIC){
      // actual intregral boundaries
      double a = 0.0;
      double b = (nhalf - 1) * h;
      // Prepare endpoint matrices: f_0, f_1, f_2, f_3, f_{nt-1}, f_{nt-2}, f_{nt-3}, f_{nt-4}
      std::vector<cdmatrix> endpoints(8);
      for (int m = 0; m <= 3; ++m) {
	get_ret(m, endpoints[m], fft_domain::time);                        // f_m
	get_ret(nhalf - 1 - m, endpoints[7 - m], fft_domain::time);        // f_{nt-1-m} ...
	// using only boundaries from t=0;
      }
      // Loop over all frequencies
      for (int k = 0; k < Nft_; ++k) {
	double omega = grid.freq_at(k);
	cdmatrix endcor(size1_,size1_);
	double corfac;
	for (int i = 0; i < size1_; ++i) {
	  for (int j = 0; j < size1_; ++j) {
	    cdouble endcor_ij;
	    std::complex<double> endpoints_ij[8];
	    for (int m = 0; m < 8; ++m) {
	      endpoints_ij[m] = endpoints[m](i, j);
	    }
	    fourier::complex_dftcor_cubic(omega, h, a, b, endpoints_ij, &endcor_ij, &corfac);	   
	    endcor(i,j)=endcor_ij;
	  }
	}
	get_ret(k, M, fft_domain::freq);
	M*=corfac; // same for all (i,j)
	M+=endcor*cdouble(cos(omega*a),sin(omega*a));
	set_ret(k, M, fft_domain::freq);	
      }
    }else{
      // method==FFT_TRAPEZ:
      cdmatrix G0;
      get_ret(0,G0,fft_domain::time);
      for (int k = 0; k < Nft_; ++k) {	
	get_ret(k, M, fft_domain::freq);
	M -= 0.5*G0; // boundary correction from t=0
	set_ret(k, M, fft_domain::freq);
      }
    }
    ret_.smul(h,fft_domain::freq);
    M=cdmatrix::Zero(size1_,size1_);
    // Note setting G=0 at nyquist point violares the ibvertibity to_freq->to_time
    set_ret(Nft_/2, M, fft_domain::freq); // Nyquist point
    // transformation of the lesser:
    les_.fft_to_freq();
    if(method==FFT_CUBIC){
      // actual intregral boundaries
      double a = -(nhalf - 1) * h;
      double b = (nhalf - 1) * h;
      // Prepare endpoint matrices: 
      //std::vector<cdmatrix> endpoints(8);
      // omitting boundary corrections, because functions decayed
      //for (int m = 0; m <= 3; ++m) {
      //	get_les(-(nhalf-1)+m, endpoints[m], fft_domain::time);       // f_m
      //	get_les(+(nhalf-1)-m, endpoints[7 - m], fft_domain::time);   // f_{nt-1-m}
      //}
      // Loop over all frequencies
      for (int k = 0; k < Nft_; ++k) {
	double omega = grid.freq_at(k);
	//cdmatrix endcor(size1_,size1_);
	double corfac;
	// Prepare per-element endpoint array for cubic correction
	// omitting boundary corrections, because functions decayed at boundary 
	std::complex<double> endpoints_ij[8];
	cdouble endcor_ij;
	fourier::complex_dftcor_cubic(omega, h, a, b, endpoints_ij, &endcor_ij, &corfac);
	// use this if you want B.C.:
	//for (int i = 0; i < size1_; ++i) {
	//  for (int j = 0; j < size1_; ++j) {
	//    std::complex<double> endpoints_ij[8];
	//    cdouble endcor_ij;
	//    for (int m = 0; m < 8; ++m) {
	//      endpoints_ij[m] = endpoints[m](i, j);
	//    }	  
	//    fourier::complex_dftcor_cubic(omega, h, a, b, endpoints_ij, &endcor_ij, &corfac);	   
	//    endcor(i, j) = endcor_ij;
	//  }	  
	//}
	get_les(k, M, fft_domain::freq);
	M *= corfac;
	//M += endcor*cdouble(cos(omega*a),sin(omega*a));
	set_les(k, M, fft_domain::freq);	
      }
    }
    les_.smul(h,fft_domain::freq);
    M=cdmatrix::Zero(size1_,size1_);
    set_les(Nft_/2, M, fft_domain::freq);
  }
  
  void herm_matrix_ness::integral_transform_to_time(double h,fft_integral_method method) {
    cdmatrix M,M1;
    int Nhalf=Nft_/2;
    grid_info grid(Nft_,h);
    double dw=grid.dw_;
    if (Nft_ < 4) {
      throw std::runtime_error("integral_transform_to_time:cubic correction requires nt >= 4");
    }    
    // transformation of the retarded:
    // A(w)  = -1/Pi Im G(w)
    // Gret(t) =  \int dw e^{-iwt} A(w) =  \int dw e^{-iwt} 1/Pi i*Im G(w) =  \int dw e^{-iwt} 1/2Pi (G(w)-G(w)^*)
    // construct a temporary array which stores [Gret(w)-Gret(w)^*]/2pi = -i Im G^R / pi
    // such that G(t) =  theta(t) * int dw gtmp(w) e^(-i*w*t)
    fft_array gtmp=retarded();
    for(int n=0;n<Nft_;n++){
      get_ret(n,M,fft_domain::freq);
      M1=M-M.adjoint();
      //M1=M;
      gtmp.set_element(n,M1,fft_domain::freq);      
    }
    gtmp.fft_to_time();
    if(method==FFT_CUBIC){
      // apply boundary corrections (not really needed because functions decay ...)
      double a = -(Nhalf - 1) * dw;
      double b = (Nhalf - 1) * dw;
      // Prepare endpoint matrices ... No B.C.!
      //std::vector<cdmatrix> endpoints(8);      
      // for (int m = 0; m <= 3; ++m) {
      //	gtmp.get_element(-Nhalf+1+m, endpoints[m], fft_domain::freq);      
      //	gtmp.get_element(Nhalf - 1 - m, endpoints[7 - m], fft_domain::freq);
      //}
      // Loop over all positive times:
      for (int n = 0; n < Nhalf; ++n) {
	double t = grid.time_at(n);
	double corfac;
	cdouble endcor_ij;
	std::complex<double> endpoints_ij[8];
	// end-corrections are igniored ...
	fourier::complex_dftcor_cubic(-t, dw, a, b, endpoints_ij,&endcor_ij, &corfac);	    	
	//cdmatrix endcor(size1_,size1_);
      	// Prepare per-element endpoint array for cubic correction
	//for (int i = 0; i < size1_; ++i) {
	//  for (int j = 0; j < size1_; ++j) {
	//    cdouble endcor_ij;
	//    std::complex<double> endpoints_ij[8];
	//    for (int m = 0; m < 8; ++m) {
	//      endpoints_ij[m] = endpoints[m](i, j);
	//    }	  
	//    fourier::complex_dftcor_cubic(-t, dw, a, b, endpoints_ij,&endcor_ij, &corfac);	  
	//    endcor(i, j) = endcor_ij;
	//  }
	//}
	gtmp.get_element(n, M, fft_domain::time);
	M*=corfac;
	//M+=endcor*cdouble(cos(-t*a),sin(-t*a));
	gtmp.set_element(n, M, fft_domain::time);	
      }
    }
    retarded().set_zero(fft_domain::time);
    for (int n = 0; n < Nhalf; ++n) {
    //for (int n = 0; n < Nft_; ++n) {
      gtmp.get_element(n, M, fft_domain::time);
      M *= (dw*0.5/M_PI);
      set_ret(n,M,fft_domain::time);
    }
    // transformation of the lesser:
    // such that G(t) =  int dw Gles(w) e^(-i*w*t) / 2pi    
    les_.fft_to_time();
    if(method==FFT_CUBIC){
      // apply boundary corrections:
      double a = -(Nhalf - 1) * dw;
      double b = (Nhalf - 1) * dw;
      // Prepare endpoint matrices
      //std::vector<cdmatrix> endpoints(8);      
      //for (int m = 0; m <= 3; ++m) {
      //	get_les(-Nhalf+1+m, endpoints[m], fft_domain::freq);      
      //get_les(Nhalf - 1 - m, endpoints[7 - m], fft_domain::freq);
      //}
      // Loop over all positive times:
      for (int n = 0; n < Nhalf; ++n) {
	double t = grid.time_at(n);
	double corfac;
	cdouble endcor_ij;
	std::complex<double> endpoints_ij[8];
	fourier::complex_dftcor_cubic(-t, dw, a, b, endpoints_ij,&endcor_ij, &corfac);
  	//cdmatrix endcor(size1_,size1_);
      	// Prepare per-element endpoint array for cubic correction
	//	for (int i = 0; i < size1_; ++i) {
	// for (int j = 0; j < size1_; ++j) {
	//  cdouble endcor_ij;
	//  std::complex<double> endpoints_ij[8];
	//  for (int m = 0; m < 8; ++m) {
	//    endpoints_ij[m] = endpoints[m](i, j);
	//  }	  
	//  fourier::complex_dftcor_cubic(-t, dw, a, b, endpoints_ij,&endcor_ij, &corfac);  
	//   endcor(i, j) = endcor_ij;
	//}
	//}
	get_les(n, M, fft_domain::time);
	M*=corfac;
	//M(0,0)=corfac;
	
	//M+=endcor*cdouble(cos(-t*a),sin(-t*a));
	set_les(n, M, fft_domain::time);	
      }
    }
    lesser().smul(dw*0.5/M_PI,fft_domain::time);
    for (int n = 1; n < Nhalf; ++n) {
      get_les(n, M, fft_domain::time);      
      set_les(-n,-M.adjoint(),fft_domain::time);
    }
    M=cdmatrix::Zero(size1_,size1_);
    //for (int n = nt_; n <= Nft_-nt_; ++n) {
    set_les(Nft_/2, M, fft_domain::time);	
    //}
  }

  void herm_matrix_ness::force_equilibrium(int bosefermi, double beta, double mu, double h) {
    // --- Step 1: Forward FFT of retarded part to frequency domain
    this->integral_transform_to_freq(h,FFT_TRAPEZ);
    // --- Step 2: Replace lesser Green’s function G^< in frequency domain
    this->les_.set_zero(fft_domain::freq); // clear previous content
    int Nft = this->Nft_;
    int size1 = this->size1_;
    grid_info grid(Nft, h);
    cdmatrix M;
    cdmatrix Gles;        
    for (int w = 0; w < Nft; ++w) {
      double omega = grid.freq_at(w);
      this->get_ret(w, M, fft_domain::freq);      
      // Imaginary part: (M - M†) / (2i)
      cdmatrix ImG = 0.5 * (M - M.adjoint()) ;      
      // Compute distribution function
      double arg = beta * (omega - mu);
      double f;
      if (bosefermi > 0) {
	// Bosons
	if (arg > 100) f = 0.0;
	else if (arg < -100) f = -1.0;
	else f = 1.0 / (std::exp(arg) - 1.0);
      } else {
	// Fermions
	if (arg > 100) f = 0.0;
	else if (arg < -100) f = 1.0;
	else f = 1.0 / (1.0 + std::exp(arg));
      }      
      // Apply fluctuation-dissipation relation:
      // G^< = -2i * f * Im[G^R]
      Gles = -(2.0*bosefermi) * f * ImG;
      this->set_les(w, Gles, fft_domain::freq);
    }
    // --- Step 3: Inverse FFT to time domain (only lesser)
    double dw=grid.dw_;
    this->les_.fft_to_time();
    lesser().smul(dw/(2.0*M_PI),fft_domain::time);
    for (int n = 1; n < Nft_/2; ++n) {
      get_les(n, M, fft_domain::time);      
      set_les(-n,-M.adjoint(),fft_domain::time);
    }
    M=cdmatrix::Zero(size1_,size1_);
    set_les(Nft_/2, M, fft_domain::time);	
  }

} // namespace ness2





#endif // NESS2_HERM_MATRIX_NESS_IMPL_HPP
