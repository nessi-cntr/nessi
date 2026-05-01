#include "catch.hpp"
#include "ness2_fft_decl.hpp"

using namespace ness2;

TEST_CASE("fft_array basic construction and resize", "[ness][fft_array][basic]") {
    int Nft = 16;
    int size1 = 2;

    fft_array A(Nft, size1);

    SECTION("Check initial sizes") {
        REQUIRE(A.Nft_ == Nft);
        REQUIRE(A.size1_ == size1);
        REQUIRE(A.element_size_ == size1 * size1);
    }

    SECTION("Resize updates dimensions correctly") {
        A.resize(32, 3);
        REQUIRE(A.Nft_ == 32);
        REQUIRE(A.size1_ == 3);
        REQUIRE(A.element_size_ == 9);
    }
}

TEST_CASE("fft_array nontrivial initialization and round-trip FFT", "[ness][fft_array][basic]") {
    int Nft = 8;
    int size1 = 2;
    fft_array A(Nft, size1);

    // Fill time-domain with deterministic nontrivial values
    for (int t = 0; t < Nft; ++t) {
        cdmatrix M(size1, size1);
        for (int i = 0; i < size1; ++i) {
            for (int j = 0; j < size1; ++j) {
                M(i, j) = cplx(t + i * 0.1, j * 0.2);
            }
        }
        A.set_element(t, M, fft_domain::time);
    }

    // Copy original time-domain for comparison
    std::vector<cdmatrix> original(Nft);
    for (int t = 0; t < Nft; ++t) {
        A.get_element(t, original[t], fft_domain::time);
    }

    // Perform forward + backward FFT
    A.fft_to_freq();
    A.fft_to_time();

    // FFTW backward transform is unnormalized → scale by 1 / Nft
    double norm_factor = 1.0 / Nft;

    for (int t = 0; t < Nft; ++t) {
        cdmatrix recovered;
        A.get_element(t, recovered, fft_domain::time);

        recovered *= norm_factor;

        REQUIRE((recovered - original[t]).norm() < 1e-10);
    }
}

TEST_CASE("fft_array text file round-trip with nontrivial data", "[ness][fft_array][basic]") {
    int Nft = 4;
    int size1 = 2;
    fft_array A(Nft, size1);

    // Fill with nontrivial values
    for (int t = 0; t < Nft; ++t) {
        cdmatrix M(size1, size1);
        for (int i = 0; i < size1; ++i) {
            for (int j = 0; j < size1; ++j) {
                M(i, j) = cplx(i + j, t * 0.1);
            }
        }
        A.set_element(t, M, fft_domain::time);
	M*=0.1;
	A.set_element(t, M, fft_domain::freq);
    }

    std::string fname = "fft_array_test.txt";
    A.print_to_file(fname);

    fft_array B(1, 1); // dummy size to trigger resize
    B.read_from_file(fname);

    for (int t = 0; t < Nft; ++t) {
        cdmatrix A_mat, B_mat;
        A.get_element(t, A_mat, fft_domain::time);
        B.get_element(t, B_mat, fft_domain::time);
        REQUIRE((A_mat - B_mat).norm() < 1e-12);
        A.get_element(t, A_mat, fft_domain::freq);
        B.get_element(t, B_mat, fft_domain::freq);
        REQUIRE((A_mat - B_mat).norm() < 1e-12);
    }
}

#if CNTR_USE_HDF5 == 1

#include "hdf5.h"

TEST_CASE("fft_array HDF5 round-trip with nontrivial data", "[ness][fft_array][hdf5]") {
    int Nft = 4;
    int size1 = 2;
    fft_array A(Nft, size1);

    // Fill with nontrivial values
    for (int t = 0; t < Nft; ++t) {
        cdmatrix M(size1, size1);
        for (int i = 0; i < size1; ++i) {
            for (int j = 0; j < size1; ++j) {
                M(i, j) = cplx(i + j + 0.5, t * 0.1);
            }
        }
        A.set_element(t, M, fft_domain::time);
        M *= 0.2;
        A.set_element(t, M, fft_domain::freq);
    }

    const std::string fname = "fft_array_test.h5";
    hid_t file_id = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    REQUIRE(file_id >= 0);

    hid_t group_id = H5Gcreate(file_id, "/fft_data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    REQUIRE(group_id >= 0);

    A.write_to_hdf5(group_id);

    H5Gclose(group_id);
    H5Fclose(file_id);

    // Read into a new object
    fft_array B(1, 1); // dummy, will be resized
    file_id = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    REQUIRE(file_id >= 0);

    group_id = H5Gopen(file_id, "/fft_data", H5P_DEFAULT);
    REQUIRE(group_id >= 0);

    B.read_from_hdf5(group_id);

    H5Gclose(group_id);
    H5Fclose(file_id);

    // Compare
    for (int t = 0; t < Nft; ++t) {
        cdmatrix A_mat, B_mat;
        A.get_element(t, A_mat, fft_domain::time);
        B.get_element(t, B_mat, fft_domain::time);
        REQUIRE((A_mat - B_mat).norm() < 1e-12);

        A.get_element(t, A_mat, fft_domain::freq);
        B.get_element(t, B_mat, fft_domain::freq);
        REQUIRE((A_mat - B_mat).norm() < 1e-12);
    }
}


#endif
