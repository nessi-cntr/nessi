#include "catch.hpp"
#include "ness2.hpp"
#include <iostream>
#include <fstream>
#include <complex>
#include <string>


using namespace ness2;


TEST_CASE("herm_matrix_ness integral_transform_to_time|freq roundtrip", "[ness][herm_matrix_ness][integral_transform]") {
    // === Read input parameters ===
    double h= 0.01;
    int Nft = 16384;
    int size = 1;
    double tol = 1e-3;

    // ===== FFT TRAPEZ ======
    // Initialize ret in time + round trip
    herm_matrix_ness S(Nft, size);
    // === Step 1: Initialize Sret with deterministic data ===
    int nt=Nft/2-1;
    for (int t = 0; t <= nt; ++t) {
      cdmatrix Mret(size, size);
      Mret(0,0) = cdouble( sin(h*t),cos(h*t))*exp(-h*h*t*t);
      if (t == 0) {
    cdmatrix tmp = Mret - Mret.adjoint();
    Mret=tmp;
      }
      S.set_ret(t, Mret, fft_domain::time);
    }
    // === roundtrip  ===
    S.integral_transform_to_freq(h,FFT_TRAPEZ);
    herm_matrix_ness S1=S;
    S1.integral_transform_to_time(h,FFT_TRAPEZ);
    double err = distance_norm2(S.ret_,S1.ret_,fft_domain::time);
    REQUIRE(err/Nft<tol);
    err = 0.0;
    
    // intialize ret in freq + round trip
    herm_matrix_ness G(Nft, size);

    // === Step 1: Initialize Gret.freq with deterministic data ===
    for (int w = 0; w < Nft; ++w) {
      cdmatrix Mret(size, size);
      Mret(0,0) = cdouble( sin(h*w),cos(h*w))*exp(-h*h*w*w); // any functiom
      G.set_ret(w, Mret, fft_domain::freq);
    }
    // === roundtrip  ===
    G.integral_transform_to_time(h,FFT_TRAPEZ);
    herm_matrix_ness G1=G;
    G1.integral_transform_to_freq(h,FFT_TRAPEZ);
    err = distance_norm2(G.ret_,G1.ret_,fft_domain::freq);
    REQUIRE(err/Nft<tol);
    err=0.0;
    
    S.ret_.set_zero(fft_domain::time);
    S.ret_.set_zero(fft_domain::freq);
    G.ret_.set_zero(fft_domain::time);
    G.ret_.set_zero(fft_domain::freq);
    S1.ret_.set_zero(fft_domain::time);
    S1.ret_.set_zero(fft_domain::freq);
    G1.ret_.set_zero(fft_domain::time);
    G1.ret_.set_zero(fft_domain::freq);
    
    // Initialize les in time + round trip;
    // === Step 1: Initialize Sles with deterministic data ===
    for (int t = 0; t <= nt; ++t) {
      cdmatrix Mles(size, size);
      Mles(0,0) = cdouble( sin(h*t),cos(h*t))*exp(-h*h*t*t);
      if (t == 0) {
    cdmatrix tmp = Mles - Mles.adjoint();
    Mles=tmp;
      }
      S.set_les(t, Mles, fft_domain::time);
    }
    // === roundtrip  ===
    S.integral_transform_to_freq(h,FFT_TRAPEZ);
    S1=S;
    S1.integral_transform_to_time(h,FFT_TRAPEZ);
    err = distance_norm2(S.les_,S1.les_,fft_domain::time);
    REQUIRE(err/Nft<tol);
    err = 0.0;
    
    // intialize les in freq + round trip
    // === Step 1: Initialize Gles.freq with deterministic data ===
    for (int w = 0; w < Nft; ++w) {
      cdmatrix Mles(size, size);
      Mles(0,0) = cdouble( sin(h*w),cos(h*w))*exp(-h*h*w*w); // any functiom
      G.set_les(w, Mles, fft_domain::freq);
    }
    // === roundtrip  ===
    G.integral_transform_to_time(h,FFT_TRAPEZ);
    G1=G;
    G1.integral_transform_to_freq(h,FFT_TRAPEZ);
    err = distance_norm2(G.les_,G1.les_,fft_domain::freq);
    REQUIRE(err/Nft<tol);
    err = 0.0;

    
    // ===== FFT CUBIC ======
    // Initialize ret in time + round trip
    // === Step 1: Initialize Sret with deterministic data ===
    for (int t = 0; t <= nt; ++t) {
      cdmatrix Mret(size, size);
      Mret(0,0) = cdouble( sin(h*t),cos(h*t))*exp(-h*h*t*t);
      if (t == 0) {
    cdmatrix tmp = Mret - Mret.adjoint();
    Mret=tmp;
      }
      S.set_ret(t, Mret, fft_domain::time);
    }
    // === roundtrip  ===
    S.integral_transform_to_freq(h,FFT_CUBIC);
    S1=S;
    S1.integral_transform_to_time(h,FFT_CUBIC);
     err = distance_norm2(S.ret_,S1.ret_,fft_domain::time);
    REQUIRE(err/Nft<tol);
    err = 0.0;
    
    // intialize ret in freq + round trip
    // === Step 1: Initialize Gret.freq with deterministic data ===
    for (int w = 0; w < Nft; ++w) {
      cdmatrix Mret(size, size);
      Mret(0,0) = cdouble( sin(h*w),cos(h*w))*exp(-h*h*w*w); // any functiom
      G.set_ret(w, Mret, fft_domain::freq);
    }
    // === roundtrip  ===
    G.integral_transform_to_time(h,FFT_CUBIC);
    G1=G;
    G1.integral_transform_to_freq(h,FFT_CUBIC);
    err = distance_norm2(G.ret_,G1.ret_,fft_domain::freq);
    REQUIRE(err/Nft<tol);
    err=0.0;
    
    // Initialize les in time + round trip;
    // === Step 1: Initialize Sles with deterministic data ===
    for (int t = 0; t <= nt; ++t) {
      cdmatrix Mles(size, size);
      Mles(0,0) = cdouble( sin(h*t),cos(h*t))*exp(-h*h*t*t);
      if (t == 0) {
    cdmatrix tmp = Mles - Mles.adjoint();
    Mles=tmp;
      }
      S.set_les(t, Mles, fft_domain::time);
    }
    // === roundtrip  ===
    S.integral_transform_to_freq(h,FFT_CUBIC);
    S1=S;
    S1.integral_transform_to_time(h,FFT_CUBIC);
    err = distance_norm2(S.les_,S1.les_,fft_domain::time);
    REQUIRE(err/Nft<tol);
    err = 0.0;
    
    // intialize les in freq + round trip
    // === Step 1: Initialize Gles.freq with deterministic data ===
    for (int w = 0; w < Nft; ++w) {
      cdmatrix Mles(size, size);
      Mles(0,0) = cdouble( sin(h*w),cos(h*w))*exp(-h*h*w*w); // any functiom
      G.set_les(w, Mles, fft_domain::freq);
    }
    // === roundtrip  ===
    G.integral_transform_to_time(h,FFT_CUBIC);
    G1=G;
    G1.integral_transform_to_freq(h,FFT_CUBIC);
    err = distance_norm2(G.les_,G1.les_,fft_domain::freq);
    REQUIRE(err/Nft<tol);
    err = 0.0;
}
