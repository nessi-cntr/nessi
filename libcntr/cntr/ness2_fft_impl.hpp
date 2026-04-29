#ifndef FFT_NESS2_FFT_IMPL_HPP
#define FFT_NESS2_FFT_IMPL_HPP

#include "ness2_fft_decl.hpp"
#include <cstring>
#include <stdexcept>
#include <complex>

namespace ness2 {



  grid_info::grid_info(int Nft, double h)
    : Nft_(Nft), h_(h)
  {
    if (Nft_ <= 0 || h_ <= 0.0 || (Nft_/2)*2!=Nft_ )
      throw std::invalid_argument("Nft and h must be positive, and Nft must be even");
    dw_ = M_PI / ((Nft_/2) * h_);
  }
  
  void grid_info::set_h(double new_h) {
    if (new_h <= 0.0)
      throw std::invalid_argument("h must be positive");
    h_ = new_h;
    dw_ = M_PI / ((Nft_/2) * h_);
  }
  
  void grid_info::set_dw(double new_dw) {
    if (new_dw <= 0.0)
      throw std::invalid_argument("dw must be positive");
    dw_ = new_dw;
    h_ = M_PI / ((Nft_/2) * dw_);
  }
  
  double grid_info::time_at(int n) const {
    n = ((n % Nft_) + Nft_) % Nft_;
    if (n < (Nft_/2)) return n * h_;
    return (n - Nft_) * h_;  // negative time: -(2 nt - n) * h
  }
  
  double grid_info::freq_at(int k) const {
    k = ((k % Nft_) + Nft_) % Nft_;
    if (k < (Nft_/2)) return k * dw_;    
    return (k - Nft_) * dw_;      
  }
  
  std::vector<double> grid_info::time_grid() const {
    std::vector<double> tgrid(Nft_);
    for (int n = 0; n < Nft_; ++n) {
      tgrid[n] = time_at(n);
    }
    return tgrid;
  }
  
  std::vector<double> grid_info::freq_grid() const {
    std::vector<double> fgrid(Nft_);
    for (int k = 0; k < Nft_; ++k) {
      fgrid[k] = freq_at(k);
    }
    return fgrid;
  }

  bool FFT_OMP_Manager::initialized_ = false;
  int  FFT_OMP_Manager::current_threads_ = 1;
  
  void FFT_OMP_Manager::initialize() {
#ifdef NESS2_FFT_USE_OMP
    if (initialized_) { return; }
    if (!fftw_init_threads())
      throw std::runtime_error("FFTW: fftw_init_threads() failed");
    fftw_plan_with_nthreads(1);       // affects plans created *after* this call
    current_threads_ = 1;
    initialized_ = true;
#else
    // No-op for OMP or single-thread builds
    // (void)n;
    initialized_ = true;
#endif
  }
  
  void FFT_OMP_Manager::set_threads_for_new_plans(int n) {
#ifdef NESS2_FFT_USE_OMP
    fftw_plan_with_nthreads(n);       // only affects *future* plans
    current_threads_ = n;
#else
    (void)n;                          // no-op if not using fftw3_threads
#endif
  }
  
  void FFT_OMP_Manager::finalize() {
#ifdef NESS2_FFT_USE_OMP
    if (initialized_) {
      fftw_cleanup_threads();         // call only after *all* plans destroyed
      initialized_ = false;
    }
#else
    initialized_ = false;
#endif
  }
  
  


  
  
  fft_array::fft_array(int nt, int size, unsigned FFTW_FLAG)
    : Nft_(nt), size1_(size), element_size_(size * size)
  {
    int total_size = Nft_ * element_size_;
    time_ = fftw_alloc_complex(total_size);
    freq_ = fftw_alloc_complex(total_size);
    if (!time_ || !freq_) throw std::runtime_error("Failed to allocate memory.");
    // Zero-initialize memory blocks
    std::memset(time_, 0, sizeof(fftw_complex) * total_size);
    std::memset(freq_, 0, sizeof(fftw_complex) * total_size);
    
    int n[] = { Nft_ };
    plan_for_ = fftw_plan_many_dft(
        1, n, element_size_,
        time_, nullptr, element_size_, 1,
        freq_, nullptr, element_size_, 1,
        FFTW_BACKWARD, FFTW_FLAG
				   );
    plan_back_ = fftw_plan_many_dft(
				    1, n, element_size_,
				    freq_, nullptr, element_size_, 1,
				    time_, nullptr, element_size_, 1,
				    FFTW_FORWARD, FFTW_FLAG
				    );
    if (!plan_for_ || !plan_back_)
      throw std::runtime_error("Failed to create FFTW plans.");
  }
  
  //fft_array::fft_array()   : Nft_(0), size1_(0), element_size_(0), time_(0),freq_(0)  {  }

  fft_array::~fft_array() {
    if (plan_for_) fftw_destroy_plan(plan_for_);
    if (plan_back_) fftw_destroy_plan(plan_back_);
    if (time_) fftw_free(time_);
    if (freq_) fftw_free(freq_);
  }
  
  fft_array::fft_array(const fft_array& other)
    : Nft_(other.Nft_), size1_(other.size1_), element_size_(other.element_size_)
  {
    int total_size = Nft_ * element_size_;
    time_ = fftw_alloc_complex(total_size);
    freq_ = fftw_alloc_complex(total_size);
    std::memcpy(time_, other.time_, sizeof(fftw_complex) * total_size);
    std::memcpy(freq_, other.freq_, sizeof(fftw_complex) * total_size);
    
    int n[] = { Nft_ };
    plan_for_ = fftw_plan_many_dft(
				   1, n, element_size_,
				   time_, nullptr, element_size_, 1,
				   freq_, nullptr, element_size_, 1,
				   FFTW_BACKWARD, FFTW_ESTIMATE
				   );
    plan_back_ = fftw_plan_many_dft(
				    1, n, element_size_,
				    freq_, nullptr, element_size_, 1,
				    time_, nullptr, element_size_, 1,
				    FFTW_FORWARD, FFTW_ESTIMATE
				    );
  }
  
  fft_array& fft_array::operator=(const fft_array& other) {
    if (this == &other) return *this;

    int total_size = other.Nft_ * other.element_size_;

       
    // Allocate new memory
    fftw_complex* new_time = fftw_alloc_complex(total_size);
    fftw_complex* new_freq = fftw_alloc_complex(total_size);
    if (!new_time || !new_freq) {
      if (new_time) fftw_free(new_time);
      if (new_freq) fftw_free(new_freq);
      throw std::runtime_error("fft_array::operator=: failed to allocate memory.");
    }

    std::memcpy(new_time, other.time_, sizeof(fftw_complex) * total_size);
    std::memcpy(new_freq, other.freq_, sizeof(fftw_complex) * total_size);
    
    int n[] = { other.Nft_ };
    fftw_plan new_plan_for = fftw_plan_many_dft(
						1, n, other.element_size_,
						new_time, nullptr, other.element_size_, 1,
						new_freq, nullptr, other.element_size_, 1,
						FFTW_BACKWARD, FFTW_ESTIMATE
						);
    fftw_plan new_plan_back = fftw_plan_many_dft(
						 1, n, other.element_size_,
						 new_freq, nullptr, other.element_size_, 1,
						 new_time, nullptr, other.element_size_, 1,
						 FFTW_FORWARD, FFTW_ESTIMATE
						 );
    
    if (!new_plan_for || !new_plan_back) {
      if (new_plan_for) fftw_destroy_plan(new_plan_for);
        if (new_plan_back) fftw_destroy_plan(new_plan_back);
        fftw_free(new_time);
        fftw_free(new_freq);
        throw std::runtime_error("fft_array::operator=: failed to create FFTW plans.");
    }
    
    // Free old resources only after successful setup
    if (plan_for_) fftw_destroy_plan(plan_for_);
    if (plan_back_) fftw_destroy_plan(plan_back_);
    if (time_) fftw_free(time_);
    if (freq_) fftw_free(freq_);
    
    // Assign new resources
    Nft_ = other.Nft_;
    size1_ = other.size1_;
    element_size_ = other.element_size_;
    time_ = new_time;
    freq_ = new_freq;
    plan_for_ = new_plan_for;
    plan_back_ = new_plan_back;
    
    return *this;
  }


  void fft_array::resize(int new_Nft, int new_size1, unsigned FFTW_FLAG) {
    // Destroy old FFTW plans
    if (plan_for_) fftw_destroy_plan(plan_for_);
    if (plan_back_) fftw_destroy_plan(plan_back_);
    
    // Free old memory
    if (time_) fftw_free(time_);
    if (freq_) fftw_free(freq_);
    
    // Update sizes
    Nft_ = new_Nft;
    size1_ = new_size1;
    element_size_ = size1_ * size1_;
    
    int total_size = Nft_ * element_size_;
    
    // Allocate new memory
    time_ = fftw_alloc_complex(total_size);
    freq_ = fftw_alloc_complex(total_size);
    if (!time_ || !freq_) {
      throw std::runtime_error("fft_array::resize: memory allocation failed.");
    }
    
    // Create new FFTW plans
    int n[] = { Nft_ };
    plan_for_ = fftw_plan_many_dft(
				   1, n, element_size_,
				   time_, nullptr, element_size_, 1,
				   freq_, nullptr, element_size_, 1,
				   FFTW_BACKWARD, FFTW_FLAG
				   );
    plan_back_ = fftw_plan_many_dft(
				    1, n, element_size_,
				    freq_, nullptr, element_size_, 1,
				    time_, nullptr, element_size_, 1,
				    FFTW_FORWARD, FFTW_FLAG
				    );
    
    if (!plan_for_ || !plan_back_) {
      throw std::runtime_error("fft_array::resize: FFTW plan creation failed.");
    }
  }
  

  
  void fft_array::fft_to_freq() {
    fftw_execute(plan_for_);
  }
  
  void fft_array::fft_to_time() {
    fftw_execute(plan_back_);
  }
  
  void fft_array::incr(const fft_array& B, std::complex<double> alpha, fft_domain type) {
    if (Nft_ != B.Nft_ || element_size_ != B.element_size_) {
      throw std::invalid_argument("fft_array::incr: mismatched array sizes");
    }
    
    int total_size = Nft_ * element_size_;
    fftw_complex* A_data = (type == fft_domain::time) ? time_ : freq_;
    fftw_complex* B_data = (type == fft_domain::time) ? B.time_ : B.freq_;
    
    for (int i = 0; i < total_size; ++i) {
      A_data[i][0] += alpha.real() * B_data[i][0] - alpha.imag() * B_data[i][1];
      A_data[i][1] += alpha.real() * B_data[i][1] + alpha.imag() * B_data[i][0];
    }
  }
  
  void fft_array::smul(std::complex<double> alpha, fft_domain type) {
    int total_size = Nft_ * element_size_;
    fftw_complex* data = (type == fft_domain::time) ? time_ : freq_;
    
    for (int i = 0; i < total_size; ++i) {
      double re = data[i][0];
      double im = data[i][1];
      data[i][0] = alpha.real() * re - alpha.imag() * im;
      data[i][1] = alpha.real() * im + alpha.imag() * re;
    }
  }

  void fft_array::set_matrixelement(int i1, int i2, const fft_array& B, int j1, int j2, fft_domain type) {
    if (size1_ != B.size1_ || Nft_ != B.Nft_) {
      throw std::invalid_argument("fft_array::set_matrixelement: mismatched sizes");
    }
    
    int idxA = i1 * size1_ + i2;
    int idxB = j1 * size1_ + j2;
    
    if (idxA >= element_size_ || idxB >= B.element_size_) {
      throw std::out_of_range("fft_array::set_matrixelement: index out of bounds");
    }
    
    fftw_complex* A_data = (type == fft_domain::time) ? time_ : freq_;
    fftw_complex* B_data = (type == fft_domain::time) ? B.time_ : B.freq_;
    
    for (int t = 0; t < Nft_; ++t) {
      A_data[t * element_size_ + idxA][0] = B_data[t * B.element_size_ + idxB][0];
      A_data[t * element_size_ + idxA][1] = B_data[t * B.element_size_ + idxB][1];
    }
  }

  void fft_array::set_zero(fft_domain type) {
    int total_size = Nft_ * element_size_;
    fftw_complex* data = (type == fft_domain::time) ? time_ : freq_;
    std::memset(data, 0, sizeof(fftw_complex) * total_size);
  }
  

  void fft_array::set_element(int t, cplx val, fft_domain type) {
    if (size1_ != 1) {
      throw std::logic_error("set_element(cplx) called on non-scalar fft_array");
    }
    // Wrap t periodically
    t = ((t % Nft_) + Nft_) % Nft_;
    fftw_complex* data = (type == fft_domain::time) ? time_ : freq_;
    data[t * element_size_][0] = val.real();
    data[t * element_size_][1] = val.imag();
  }

  void fft_array::get_element(int t, cplx& val, fft_domain type) const {
    if (size1_ != 1) {
      throw std::logic_error("get_element(cplx) called on non-scalar fft_array");
    }
    // Wrap t periodically
    t = ((t % Nft_) + Nft_) % Nft_;
    const fftw_complex* data = (type == fft_domain::time) ? time_ : freq_;
    val = std::complex<double>(data[t * element_size_][0], data[t * element_size_][1]);
  }
   
  
  // Explicit instantiations for cdmatrix
  template void fft_array::left_multiply<cdmatrix>(const cdmatrix& M, fft_domain type);
  template void fft_array::right_multiply<cdmatrix>(const cdmatrix& M, fft_domain type);
  template void fft_array::left_multiply_hermconj<cdmatrix>(const cdmatrix& M, fft_domain type);
  template void fft_array::right_multiply_hermconj<cdmatrix>(const cdmatrix& M, fft_domain type);
  
  // Explicit instantiations for cdmatrix
  template void fft_array::set_element<cdmatrix>(int t, const cdmatrix& M, fft_domain type);
  template void fft_array::get_element<cdmatrix>(int t, cdmatrix& M, fft_domain type) const;

  void fft_array::print_to_file(const std::string& filename, int precision) const {
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
      throw std::runtime_error("Failed to open output file: " + filename);
    }
    
    outfile << "# " << Nft_ << " " << size1_  << "\n";
    
    auto dump_block = [&](const std::string& label, fft_domain type) {
      for (int t = 0; t < Nft_; ++t) {
	outfile << label << ": " << t;
	
	cdmatrix M;
	this->get_element(t, M, type);
	for (int i = 0; i < size1_; ++i) {
	  for (int j = 0; j < size1_; ++j) {
	    outfile << " " << std::setprecision(precision) << M(i, j).real()
		    << " " << std::setprecision(precision) << M(i, j).imag();
	  }
	}
	outfile << "\n";
      }
    };
    
    dump_block("time", fft_domain::time);
    dump_block("freq", fft_domain::freq);

    outfile.close();
  }
  
  void fft_array::read_from_file(const std::string& filename, unsigned FFTW_FLAG) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
      throw std::runtime_error("Failed to open input file: " + filename);
    }
    
    std::string line;
    std::getline(infile, line);
    std::istringstream header(line.substr(1)); // skip '#'
    int file_Nft, file_size1;
    header >> file_Nft >> file_size1 ;
    
    // Resize if needed
    if (file_Nft != Nft_ || file_size1 != size1_) {
      resize(file_Nft, file_size1, FFTW_FLAG);
    }

    auto load_block = [&](const std::string& label, fft_domain type) {
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
	for (int i = 0; i < file_size1; ++i) {
	  for (int j = 0; j < file_size1; ++j) {
	    double re, im;
	    if (!(linestream >> re >> im)) {
	      throw std::runtime_error("Failed to read complex value in " + label);
	    }
	    M(i, j) = cplx(re, im);
	  }
	}
	this->set_element(t, M, type);
      }
    };
    
    load_block("time", fft_domain::time);
    load_block("freq", fft_domain::freq);

    infile.close();
  }


#if CNTR_USE_HDF5 == 1

  void fft_array::write_to_hdf5(hid_t group_id)  const {
    store_int_attribute_to_hid(group_id, "Nft", Nft_);
    store_int_attribute_to_hid(group_id, "size1", size1_);
    store_int_attribute_to_hid(group_id, "element_size", element_size_);
    
    // Shape is [Nft_, size1_, size1_]
    hsize_t shape[3] = {
      static_cast<hsize_t>(Nft_),
      static_cast<hsize_t>(size1_),
      static_cast<hsize_t>(size1_)
    };
    store_cplx_array_to_hid(group_id, "time", reinterpret_cast<std::complex<double>*>(time_), shape, 3);
    store_cplx_array_to_hid(group_id, "freq", reinterpret_cast<std::complex<double>*>(freq_), shape, 3);    
  }
  
  void fft_array::read_from_hdf5(hid_t group_id, unsigned FFTW_FLAG) {
    int Nft = read_primitive_type<int>(group_id, "Nft");
    int size1 = read_primitive_type<int>(group_id, "size1");
    int element_size = read_primitive_type<int>(group_id, "element_size");
    
    if (Nft != Nft_ || size1 != size1_ || element_size != element_size_) {
      resize(Nft, size1, FFTW_FLAG);
    }

    hsize_t total_size = static_cast<hsize_t>(Nft_) * size1_ * size1_;
    read_primitive_type_array(group_id, "time", total_size, reinterpret_cast<std::complex<double>*>(time_));
    read_primitive_type_array(group_id, "freq", total_size, reinterpret_cast<std::complex<double>*>(freq_));
  }
  
#endif

  
  
} // namespace ness2



#endif

