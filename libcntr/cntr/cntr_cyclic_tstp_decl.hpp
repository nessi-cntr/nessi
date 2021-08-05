#ifndef CNTR_CCLC_TSTP_DECL
#define CNTR_CCLC_TSTP_DECL

namespace cntr{

template <typename T> class cyclic_timestep {
 public:
	typedef std::complex<T> cplx;
	/* construction, destruction */
	cyclic_timestep();
	~cyclic_timestep();
	cyclic_timestep(int tstp,int ntau,int size1=1);
	cyclic_timestep(const cyclic_timestep &g);
	cyclic_timestep & operator=(const cyclic_timestep &g);
	void resize(int nt,int ntau,int size1);
	void clear(void);
	cplx* data_;
	int tstp_;
	int ntau_;
	int size1_;
	int size2_;
	int element_size_;
	int total_size_;
};

} //namespace cntr

#endif  // CNTR_CCLC_TSTP_DECL
