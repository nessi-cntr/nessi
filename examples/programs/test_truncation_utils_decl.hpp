#ifndef TRUNCATION_TEST_UTILS_DECL
#define TRUNCATION_TEST_UTILS_DECL

#include <sys/stat.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <chrono>

// contour library headers
#include "cntr/cntr.hpp"
#include "cntr/utils/read_inputfile.hpp"

using namespace std;

#define CPLX complex<double>
#define GREEN cntr::herm_matrix<double>
#define CFUNC cntr::function<double>
#define GTRUNC cntr::herm_matrix_moving<double>

void calculate_weights_linear_eps(int ntheta,vector<double>& wk_,vector<double>& eps_k,double max_energy,double min_energy);

void calculate_weights_sin_square(int ntheta,vector<double>& wk_,vector<double>& eps_k );

CPLX timeder(GREEN & G, int up_down,int tstp, int order,double h);

CPLX timeder(GTRUNC & G, int tstp_max,int up_down,int tstp, int order,double h);

void print_energy(GREEN & G,GREEN & Sigma,int nt,int kt,double beta,double h,string folder_out,vector <double>& Ekin, string corpus,bool & prima_volta_energy,int tstp_max,int time_interval,int tc,GTRUNC & G_t,vector <double> & n_val,vector <double> & n_val_t);


void print_n( vector<double> & n_val, vector<double> n_val_t,  vector<double> & elapsed_time, vector<double> elapsed_time_t,  double u0,double u1, int tmax,double h,int tc,int uf_prec,int nt,string folder,string corpus);

void print_nk( vector <vector<double>> & nk_val, vector <vector<double>> &nk_val_t, int neps,vector <double> & ek,double u0,double u1, int tmax,double h,int tc,int uf_prec,int nt,string folder);

void get_sigma(int n,GREEN &Sigma,CFUNC &ufunc,GREEN &G);

void get_sigma(GTRUNC &Sigma,cntr::function_moving<double> &ufunc,GTRUNC &G,GTRUNC &chi);

#endif //TRUNCATION_TEST_UTILS_DECL
