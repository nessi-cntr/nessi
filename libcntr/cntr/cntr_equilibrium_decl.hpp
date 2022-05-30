#ifndef CNTR_EQUILIBRIUM_DECL_H
#define CNTR_EQUILIBRIUM_DECL_H

#include "cntr_global_settings.hpp"

namespace cntr {

  template <typename T> class function;
  template <typename T> class herm_matrix;
  template <typename T> class herm_matrix_timestep;
  template <typename T> class herm_pseudo;

  /*###########################################################################################
    #
    #   COMPUTATION OF EQUILIBRIUM GREEN FUNCTIONS FROM DOS
    #   (for size1>1, G is set to be a diagobal matrix)
    #   dome predefined DOS, otherwise:
    #   dos must have an operator (omega) (the value) and functions lo() (lower
    cutoff)
    #   and hi() (upper cutoff)
    #
    ###########################################################################################*/
  template <typename T>
  T fermi(T beta, T omega);
/// @private
  template <typename T>
  T fermi_exp(T beta, T tau, T omega);

  template<typename T>
  dvector fermi(T beta,dvector &omega);
/// @private
  template<typename T>
  dvector fermi_exp(T beta,T tau,dvector &omega);

  /// @private
  template<typename T>
  cdmatrix diag_prop(T time,dvector &omega);


  // BETHE DOS, bandwidh 4
  template <typename T>
  void green_equilibrium_mat_bethe(herm_matrix<T> &G, double beta,
    int limit = 100, int nn = 20,double mu=0.0,double Eshift=0.0);
  template <typename T>
  void green_equilibrium_bethe(herm_matrix<T> &G, double beta, double h,
    int limit = 100, int nn = 20,double mu=0.0,double Eshift=0.0);
  // user-defined DOS
  template <typename T, class dos_function>
  void green_equilibrium_mat(herm_matrix<T> &G, dos_function &dos, double beta,
    int limit = 100, int nn = 20,double mu=0.0);

  /// @private
  template <typename T, class dos_function>
  void green_equilibrium(herm_matrix<T> &G, dos_function &dos, double beta,
    double h, int limit = 100, int nn = 20, double mu=0.0);

  template <typename T, class dos_function>
  void green_equilibrium(herm_matrix<T> &G, dos_function &dos, double beta,
    double h, double mu=0.0, int limit = 100, int nn = 20);

  // other "simple" Greenfunctions: [idt + mu - H(t)]^{-1} and [idt + mu -
  // H0]^{-1} etc
  // template <typename T>
  // void green_transform(herm_matrix<T> &G0, herm_matrix<T> &G1,
  //   std::complex<T> *u01_0, std::complex<T> *u01_t);
  // template <typename T>
  // void green_transform(herm_matrix<T> &G0, herm_matrix<T> &G1,
  //   std::complex<T> *u01_0);

  //template <typename T>
  //void green_from_eps(herm_matrix<T> &G, T mu, T *eps, T beta, T h);
  template <typename T>
  void green_from_H(herm_matrix<T> &G,T mu,cdmatrix &eps,T beta,T h);
  template <typename T>
  void green_from_H(int tstp, herm_matrix_timestep<T> &G,T mu,cdmatrix &eps,T beta,T h);
  template <typename T>
  void green_from_H(int tstp, herm_matrix<T> &G,T mu,cdmatrix &eps,T beta,T h);
  template <typename T>
  void green_from_H(herm_matrix<T> &G,T mu,cntr::function<T> &eps,
           T beta,T h,int SolveOrder=MAX_SOLVE_ORDER,int cf_order=4);
  template <typename T>
  void green_from_H(int tstp, herm_matrix_timestep<T> &G,T mu,cntr::function<T> &eps,
           T beta,T h,bool fixHam=false,int SolveOrder=MAX_SOLVE_ORDER,int cf_order=4);
  template <typename T>
  void green_from_H(int tstp, herm_matrix<T> &G,T mu,cntr::function<T> &eps,
           T beta,T h,bool fixHam=false,int SolveOrder=MAX_SOLVE_ORDER,int cf_order=4);

  /// @private
  template <typename T>
  void green_from_H(herm_matrix_timestep<T> &G,T mu,cdmatrix &eps,T beta,T h);
  /// @private
  template <typename T>
  void green_from_H(herm_matrix<T> &G,T mu,cntr::function<T> &eps,
			     T beta,T h,bool fixHam=false,int SolveOrder=MAX_SOLVE_ORDER,int cf_order=4);
  /// @private
  template <typename T>
  void green_from_H(herm_matrix_timestep<T> &G,T mu,cntr::function<T> &eps,
			     T beta,T h,bool fixHam=false,int SolveOrder=MAX_SOLVE_ORDER,int cf_order=4);

  // simple bosonic GF

  // BOSONIC GREENS FUNCTION:
  // G(t,t') = - ii * <TC b(t) bdag(t') >, with H = w * b^dag b
  //template <typename T>
  //void green_single_pole_bose(herm_matrix<T> &G, T *w, T beta, T h);
  // G(t,t') = - ii * <TC X(t) X(t') > *in equilibrium*
  template <typename T>
  void green_single_pole_XX_timestep(int tstp, herm_matrix_timestep<T> &D0, T w, T beta,
				     T h);
  /// @private
  template <typename T>
  void green_single_pole_XX_timestep(int tstp, herm_matrix<T> &D0, T w, T beta,
    T h);
  /// @private
  template <typename T>
  void green_single_pole_XX(herm_matrix<T> &D0, T w, T beta, T h);

  /// @private
  template <typename T>
  void green_single_pole_XX_timestep(herm_matrix_timestep<T> &D0, T w, T beta,
             T h);

} // namespace cntr

#endif  // CNTR_EQUILIBRIUM_DECL_H
