#include "catch.hpp"
#include "ness2_fft_decl.hpp"

using namespace ness2;

TEST_CASE("fft_array matrix element set/get", "[ness][fft_array][elements]") {
    int Nft = 4;
    int size1 = 2;
    fft_array A(Nft, size1);

    cdmatrix M(size1, size1);
    M.setConstant({1.0, 0.0});

    A.set_element(0, M, fft_domain::time);

    cdmatrix result;
    A.get_element(0, result, fft_domain::time);

    REQUIRE((M - result).norm() < 1e-12);
}

TEST_CASE("fft_array scalar set/get", "[ness][fft_array][elements]") {
    int Nft = 4;
    int size1 = 1; // scalar case
    fft_array A(Nft, size1);

    cplx val = {1.23, 4.56};
    A.set_element(0, val, fft_domain::time);

    cplx result;
    A.get_element(0, result, fft_domain::time);

    REQUIRE(std::abs(val - result) < 1e-12);
}

TEST_CASE("fft_array set_matrixelement copies between arrays", "[ness][fft_array][elements]") {
    int Nft = 4;
    int size1 = 2;
    fft_array A(Nft, size1);
    fft_array B(Nft, size1);

    // Set B(0,1)(t) to 1 + i across all t
    B.set_zero(fft_domain::time);
    for (int t = 0; t < Nft; ++t) {
        cdmatrix M(size1, size1);
        M.setZero();
        M(0, 1) = cplx{1.0, 1.0};
        B.set_element(t, M, fft_domain::time);
    }

    // Copy B(0,1)(t) → A(1,0)(t)
    A.set_zero(fft_domain::time);
    A.set_matrixelement(1, 0, B, 0, 1, fft_domain::time);

    // Check A(1,0)(t) == 1 + i
    for (int t = 0; t < Nft; ++t) {
        cdmatrix M;
        A.get_element(t, M, fft_domain::time);
        REQUIRE(std::abs(M(1, 0) - cplx{1.0, 1.0}) < 1e-12);
    }
}

