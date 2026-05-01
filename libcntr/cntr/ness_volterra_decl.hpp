#ifndef NESS_VOLTERRA_H
#define NESS_VOLTERRA_H

#include "ness.hpp"
#include <vector>
#include <array>

namespace ness
{
  void dyson_ness(GF & G, GF & Sigma, const cdmatrix & H, double eta);
  void dyson_ness(GF_pair & G, GF_pair & Sigma, const cdmatrix & H, double eta);
  void vie2_ness(GF & G, const GF & F, const GF & Q);
  void vie2_ness(GF_pair & G, const GF_pair & F, const GF_pair & Q);
  int dyson(GF &invg0, GF &self, GF &g);
  int dyson_from_inv(GF &g, GF &self);
  void volterra_IBC(GF & gf, const GF & B, const GF & C, const cdmatrix & Id);
  void volterra_IABC(GF & gf, const cdmatrix & A, const GF & B, const GF & C, const cdmatrix & Id);
  void volterra_IABC(GF & gf, const GF & A, const cdmatrix & B, const GF & C, const cdmatrix & Id);
  void volterra_IABC(GF & gf, const cdmatrix & A, const GF & B, const cdmatrix & C, const cdmatrix & Id);
  void volterra_IABC(GF & gf, const GF & A, const GF & B, const GF & C, const cdmatrix & Id);
  void volterra_IABUCU(GF & gf, const cdmatrix & A, const GF & B, const GF & C, cdmatrix U, const cdmatrix & Id);
  void volterra_IALMC(GF & gf, const GF & A, const GF & L, const GF & M, const GF & C, const cdmatrix & Id);
  void volterra_IBCgtr(GF & gf, const GF & B, const GF & C, const cdmatrix & Id);


};//namespace


#endif // NESS_VOLTERRA_H
