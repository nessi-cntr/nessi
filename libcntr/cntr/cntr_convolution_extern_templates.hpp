#ifndef CNTR_CONVOLUTION_EXTERN_TEMPLATES_H
#define CNTR_CONVOLUTION_EXTERN_TEMPLATES_H

#include "cntr_convolution_decl.hpp"
#include "cntr_pseudo_convolution_decl.hpp"

namespace cntr {

  extern template
  void convolution_matsubara<double, herm_matrix<double> >(herm_matrix<double> &C,herm_matrix<double> &A,
									   herm_matrix<double> &B,integration::Integrator<double> &I, double beta);
  extern template
  void convolution_timestep<double>(int n,herm_matrix<double> &C,herm_matrix<double> &A,
						    herm_matrix<double> &Acc,herm_matrix<double> &B, herm_matrix<double> &Bcc,
						    integration::Integrator<double> &I, double beta,double h);
  extern template
  void  convolution_timestep<double>(int n,herm_matrix<double> &C,herm_matrix<double> &A,
						     herm_matrix<double> &B,integration::Integrator<double> &I,
						     double beta,double h);
  extern template
  void convolution<double>(herm_matrix<double> &C,herm_matrix<double> &A,herm_matrix<double> &Acc,
					   herm_matrix<double> &B, herm_matrix<double> &Bcc, integration::Integrator<double> &I,
					   double beta,double h);
  extern template
  void convolution_timestep<double>(int n,herm_matrix<double> &C,herm_matrix<double> &A,
						    herm_matrix<double> &Acc,function<double> &ft,
						    herm_matrix<double> &B,herm_matrix<double> &Bcc, integration::Integrator<double> &I,
						    double beta,double h);
  extern template
  void  convolution_timestep<double>(int n,herm_matrix<double> &C,herm_matrix<double> &A,function<double> &ft,
				     herm_matrix<double> &B, integration::Integrator<double> &I, double beta,double h);
  extern template
  void convolution<double>(herm_matrix<double> &C,herm_matrix<double> &A,herm_matrix<double> &Acc,
			   function<double> &ft, herm_matrix<double> &B,herm_matrix<double> &Bcc,
			   integration::Integrator<double> &I, double beta,double h);
  extern template
  void convolution_density_matrix<double, herm_matrix<double> >(int tstp,std::complex<double> *rho,
								herm_matrix<double> &A,herm_matrix<double> &Acc,
								function<double> &ft,herm_matrix<double> &B,herm_matrix<double> &Bcc,
								integration::Integrator<double> &I, double beta,double h);
  extern template
  void  convolution_density_matrix<double, herm_matrix<double> >(int tstp,std::complex<double> *rho,herm_matrix<double> &A,
								 function<double> &ft, herm_matrix<double> &B,integration::Integrator<double> &I,
								 double beta,double h);
  extern template
  void convolution_density_matrix<double, herm_matrix<double> >(int tstp,std::complex<double> *rho,herm_matrix<double> &A,
								herm_matrix<double> &Acc, herm_matrix<double> &B,herm_matrix<double> &Bcc,
								integration::Integrator<double> &I, double beta,double h);
  extern template
  void  convolution_density_matrix<double, herm_matrix<double> >(int tstp,std::complex<double> *rho,herm_matrix<double> &A,
								 herm_matrix<double> &B,integration::Integrator<double> &I,
								 double beta,double h);

  extern template
  void  convolution_density_matrix<double, herm_matrix<double> >(int tstp, cdmatrix &Cles, herm_matrix<double> &A, herm_matrix<double> &B,integration::Integrator<double> &I, double  beta, double h);

  extern template
  void convolution_les_timediag<double, herm_matrix<double> >(int tstp, cdmatrix &Cles, herm_matrix<double> &A, herm_matrix<double> &B,
							     integration::Integrator<double> &I, double beta, double h);
  
  // simplified new interfaces 

extern template void convolution_timestep<double>(int n,herm_matrix<double> &C,herm_matrix<double> &A,herm_matrix<double> &Acc,
  herm_matrix<double> &B, herm_matrix<double> &Bcc, double beta,double h, int SolveOrder);
extern template void convolution_timestep<double>(int n,herm_matrix<double> &C,herm_matrix<double> &A,herm_matrix<double> &B,
  double beta,double h, int SolveOrder);
extern template void convolution<double>(herm_matrix<double> &C,herm_matrix<double> &A,herm_matrix<double> &Acc,
  herm_matrix<double> &B, herm_matrix<double> &Bcc, double beta,double h, int SolveOrder);
extern template void convolution_timestep<double>(int n,herm_matrix<double> &C,herm_matrix<double> &A,
  herm_matrix<double> &Acc,function<double> &ft, herm_matrix<double> &B,herm_matrix<double> &Bcc, double beta,double h, int SolveOrder);
extern template void convolution_timestep<double>(int n,herm_matrix<double> &C,herm_matrix<double> &A,
  function<double> &ft,herm_matrix<double> &B, double beta,double h, int SolveOrder);
extern template void convolution<double>(herm_matrix<double> &C,herm_matrix<double> &A,herm_matrix<double> &Acc,
  function<double> &ft, herm_matrix<double> &B,herm_matrix<double> &Bcc, double beta,double h, int SolveOrder);

extern template void convolution_density_matrix<double, herm_matrix<double> >(int tstp,cdmatrix &rho,herm_matrix<double> &A,
  herm_matrix<double> &Acc, function<double> &ft,herm_matrix<double> &B,herm_matrix<double> &Bcc, double beta,double h, int SolveOrder);
extern template void convolution_density_matrix<double, herm_matrix<double> >(int tstp,cdmatrix &rho,herm_matrix<double> &A,
  function<double> &ft, herm_matrix<double> &B, double beta,double h, int SolveOrder);
extern template void convolution_density_matrix<double, herm_matrix<double> >(int tstp,cdmatrix &rho,herm_matrix<double> &A,
  herm_matrix<double> &Acc, herm_matrix<double> &B,herm_matrix<double> &Bcc, double beta,double h, int SolveOrder);
extern template void convolution_density_matrix<double, herm_matrix<double> >(int tstp,cdmatrix &rho,herm_matrix<double> &A,
  herm_matrix<double> &B, double beta,double h, int SolveOrder);


#if CNTR_USE_OMP==1

extern template
void convolution_timestep_omp<double>(int omp_num_threads, int tstp, herm_matrix<double> &C,
                              herm_matrix<double> &A, herm_matrix<double> &Acc, function<double> &ft,
                              herm_matrix<double> &B, herm_matrix<double> &Bcc,
                              integration::Integrator<double> &I, double beta, double h);

extern template 
void convolution_timestep_omp<double>(int omp_num_threads, int tstp, herm_matrix<double> &C,
                              herm_matrix<double> &A, function<double> &ft, herm_matrix<double> &B,
                              integration::Integrator<double> &I, double beta, double h);

extern template 
void convolution_matsubara_omp<double>(int omp_num_threads, herm_matrix<double> &C, herm_matrix<double> &A,
                               herm_matrix<double> &B, integration::Integrator<double> &I, double beta);

extern template 
void convolution_matsubara_omp<double>(int omp_num_threads, herm_matrix<double> &C, herm_matrix<double> &A,
                               function<double> &ft, herm_matrix<double> &B,
                               integration::Integrator<double> &I, double beta);

extern template 
void convolution_omp<double>(int omp_num_threads, herm_matrix<double> &C, herm_matrix<double> &A,
                     herm_matrix<double> &Acc, function<double> &ft, herm_matrix<double> &B,
                     herm_matrix<double> &Bcc, integration::Integrator<double> &I, double beta, double h);

extern template 
void convolution_timestep_omp<double>(int omp_num_threads, int tstp, herm_matrix<double> &C,
                              herm_matrix<double> &A, herm_matrix<double> &Acc, herm_matrix<double> &B,
                              herm_matrix<double> &Bcc, integration::Integrator<double> &I, double beta,
                              double h);

extern template 
void convolution_timestep_omp<double>(int omp_num_threads, int tstp, herm_matrix<double> &C,
                              herm_matrix<double> &A, herm_matrix<double> &B,
                              integration::Integrator<double> &I, double beta, double h);

extern template
void convolution_omp<double>(int omp_num_threads, herm_matrix<double> &C, herm_matrix<double> &A,
                     herm_matrix<double> &Acc, herm_matrix<double> &B, herm_matrix<double> &Bcc,
                     integration::Integrator<double> &I, double beta, double h);

// simplified new interfaces 
extern template
void convolution_timestep_omp<double>(int omp_num_threads, int tstp, herm_matrix<double> &C,
                              herm_matrix<double> &A, herm_matrix<double> &Acc, function<double> &ft,
                              herm_matrix<double> &B, herm_matrix<double> &Bcc,
                              double beta, double h, int SolveOrder);
extern template
void convolution_timestep_omp<double>(int omp_num_threads, int tstp, herm_matrix<double> &C,
                              herm_matrix<double> &A, function<double> &ft, herm_matrix<double> &B,
                              double beta, double h, int SolveOrder);

extern template
void convolution_omp<double>(int omp_num_threads, herm_matrix<double> &C, herm_matrix<double> &A,
                     herm_matrix<double> &Acc, function<double> &ft, herm_matrix<double> &B,
                     herm_matrix<double> &Bcc, double beta, double h, int SolveOrder);
extern template
void convolution_timestep_omp<double>(int omp_num_threads, int tstp, herm_matrix<double> &C,
                              herm_matrix<double> &A, herm_matrix<double> &Acc, herm_matrix<double> &B,
                              herm_matrix<double> &Bcc, double beta,
                              double h, int SolveOrder);
extern template
void convolution_timestep_omp<double>(int omp_num_threads, int tstp, herm_matrix<double> &C,
                              herm_matrix<double> &A, herm_matrix<double> &B,
                              double beta, double h, int SolveOrder);
extern template
void convolution_omp<double>(int omp_num_threads, herm_matrix<double> &C, herm_matrix<double> &A,
                     herm_matrix<double> &Acc, herm_matrix<double> &B, herm_matrix<double> &Bcc,
                     double beta, double h, int SolveOrder);


#endif

}  // namespace cntr

#endif  // CNTR_CONVOLUTION_EXTERN_TEMPLATES_H
