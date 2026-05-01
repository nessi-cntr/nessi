#ifndef NESS_UTILS_H
#define NESS_UTILS_H

#include "ness.hpp"
#include <fstream>


namespace ness
{
void set_les_from_ret(GF &gf, std::vector<double> fd);
void set_les_from_ret(GF_pair &gf, std::vector<double> fd);

void fgrid_from_tgrid(double & df, int & nf, double & dt, int & nt);
void tgrid_from_fgrid(double & dt, int & nt, double & df, int & nf);

void Bubble1_ness(GF &C,int c1,int c2,GF &A,int a1,int a2,GF &B,int b1,int b2);
void Bubble1_ness(GF_pair &C,int c1,int c2,GF_pair &A,int a1,int a2,GF_pair &B,int b1,int b2);
void Bubble1_ness(GF &C, GF &A, GF &B);
void Bubble1_ness(GF_pair &C, GF_pair &A, GF_pair &B);

void Bubble2_ness(GF &C,int c1,int c2,GF &A,int a1,int a2,GF &B,int b1,int b2);
void Bubble2_ness(GF_pair &C,int c1,int c2,GF_pair &A,int a1,int a2,GF_pair &B,int b1,int b2);
void Bubble2_ness(GF &C, GF &A, GF &B);
void Bubble2_ness(GF_pair &C, GF_pair &A, GF_pair &B);

void shift_fd(std::vector<double> & func, double dw, double shift);


void force_FD(GF &gf, double temp);
void force_FD(GF_pair &gf, double temp);

void force_imag(GF &gf);
void force_diag(GF &gf);


int set_bwd_from_fwd(GF & gf_bwd, GF & gf_fwd);


double GF2norm(GF &g1, GF &g2);

void outputGF(GF &G, std::string filename);

int healthCheck(GF &G);

void appendData(const std::vector<double> & data_vec, std::string fname, int prec);

void force_tFD(fft_solver & fft, GF &gf, int nfreq, double temp);


void readGF(GF &G, std::string filename);

void impose_sym(GF &gf, cdmatrix & sym);
 

};//end of namespace

#endif // NESS_UTILS_H
