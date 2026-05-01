#include "catch.hpp"
#include "ness2_fft_decl.hpp"
#include <Eigen/Dense>

using namespace ness2;

TEST_CASE("fft_array incr and smul with nontrivial data", "[ness][fft_array][ops]") {
    int Nft = 4;
    int size1 = 2;
    fft_array A(Nft, size1);
    fft_array B(Nft, size1);

    // Fill B with deterministic values: B(t)(i,j) = t + i + j + i*j
    for (int t = 0; t < Nft; ++t) {
        cdmatrix M(size1, size1);
        for (int i = 0; i < size1; ++i) {
            for (int j = 0; j < size1; ++j) {
                M(i, j) = cplx(t + i + j + i * j, 0.0);
            }
        }
        B.set_element(t, M, fft_domain::time);
    }

    A.set_zero(fft_domain::time);
    A.incr(B, {2.0, 0.0}, fft_domain::time);

    for (int t = 0; t < Nft; ++t) {
        cdmatrix A_mat, B_mat;
        A.get_element(t, A_mat, fft_domain::time);
        B.get_element(t, B_mat, fft_domain::time);

        REQUIRE((A_mat - 2.0 * B_mat).norm() < 1e-12);
    }

    A.smul({0.5, 0.0}, fft_domain::time);

    for (int t = 0; t < Nft; ++t) {
        cdmatrix A_mat, B_mat;
        A.get_element(t, A_mat, fft_domain::time);
        B.get_element(t, B_mat, fft_domain::time);

        REQUIRE((A_mat - B_mat).norm() < 1e-12);
    }
}

TEST_CASE("fft_array left multiply with nontrivial data", "[ness][fft_array][ops][left]") {
    int Nft = 4;
    int size1 = 2;
    fft_array A(Nft, size1);
    // Fill A with patterned values: A(t)(i,j) = 1 + i + j
    for (int t = 0; t < Nft; ++t) {
        cdmatrix M(size1, size1);
        for (int i = 0; i < size1; ++i) {
            for (int j = 0; j < size1; ++j) {
                M(i, j) = cplx(1 + i + j, 1.0*t);
            }
        }
        A.set_element(t, M, fft_domain::time);
    }

    // Define swap matrix M: [[0, 1], [1, 0]], M * M = I
    Eigen::MatrixXcd M(size1, size1);
    M << 0, 1,
         1, 0;

    fft_array B = A;  // copy original
    B.left_multiply(M, fft_domain::time);
    B.left_multiply(M, fft_domain::time);  // M * (M * A) = (M * M) * A = A

    for (int t = 0; t < Nft; ++t) {
        cdmatrix A_mat, B_mat;
        A.get_element(t, A_mat, fft_domain::time);
        B.get_element(t, B_mat, fft_domain::time);

        REQUIRE((A_mat - B_mat).norm() < 1e-12);
    }
    
}

TEST_CASE("fft_array right multiply with nontrivial data", "[ness][fft_array][ops][right]") {
    int Nft = 4;
    int size1 = 2;
    fft_array A(Nft, size1);

    // Fill A with patterned values: A(t)(i,j) = 1 + i + j
    for (int t = 0; t < Nft; ++t) {
        cdmatrix M(size1, size1);
        for (int i = 0; i < size1; ++i) {
            for (int j = 0; j < size1; ++j) {
                M(i, j) = cplx(1 + i + j, 0.0);
            }
        }
        A.set_element(t, M, fft_domain::time);
    }

    // Define swap matrix M: [[0, 1], [1, 0]], M * M = I
    Eigen::MatrixXcd M(size1, size1);
    M << 0, 1,
         1, 0;

    fft_array B = A;  // copy original
    B.right_multiply(M, fft_domain::time);
    B.right_multiply(M, fft_domain::time);  // (A * M) * M = A * (M * M) = A

    for (int t = 0; t < Nft; ++t) {
        cdmatrix A_mat, B_mat;
        A.get_element(t, A_mat, fft_domain::time);
        B.get_element(t, B_mat, fft_domain::time);

        REQUIRE((A_mat - B_mat).norm() < 1e-12);
    }
}

TEST_CASE("fft_array set_zero", "[ness][fft_array][ops][zero]") {
    int Nft = 4;
    int size1 = 2;
    fft_array A(Nft, size1);
    fft_array B(Nft, size1); // dummy with zeros
    
    // Fill A with patterned values: A(t)(i,j) = 1 + i + j
    for (int t = 0; t < Nft; ++t) {
        cdmatrix M(size1, size1);
        for (int i = 0; i < size1; ++i) {
            for (int j = 0; j < size1; ++j) {
                M(i, j) = cplx(1 + i + j, 0.0);
            }
        }
        A.set_element(t, M, fft_domain::time);
        A.set_element(t, M, fft_domain::time);
    }
    
    A.set_zero(fft_domain::time);
    A.set_zero(fft_domain::freq);
    
    for (int t = 0; t < Nft; ++t) {
        cdmatrix A_mat, B_mat;
        A.get_element(t, A_mat, fft_domain::time);
        A.get_element(t, B_mat, fft_domain::freq);
        REQUIRE(A_mat.norm() < 1e-12);
        REQUIRE(B_mat.norm() < 1e-12);
    }
}
