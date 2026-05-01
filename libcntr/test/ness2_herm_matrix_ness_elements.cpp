#include "catch.hpp"
#include "ness2.hpp"

using namespace ness2;

TEST_CASE("herm_matrix_ness matrix element set/get", "[ness][herm_matrix_ness][elements]") {
    int Nft = 4;
    int size1 = 2;
    herm_matrix_ness A(Nft, size1);

    cdmatrix M(size1, size1);
    M.setConstant({1.0, 0.0});

    // ret
    A.set_ret(0, M, fft_domain::time);
    cdmatrix result;
    A.get_ret(0, result, fft_domain::time);
    REQUIRE((M - result).norm() < 1e-12);
    
    // les
    A.set_les(0, M, fft_domain::time);
    A.get_les(0, result, fft_domain::time);
    REQUIRE((M - result).norm() < 1e-12);
    
}

TEST_CASE("herm_matrix_ness scalar set/get", "[ness][herm_matrix_ness][elements]") {
    int Nft = 4;
    int size1 = 1; // scalar case
    herm_matrix_ness A(Nft, size1);

    cplx val = {1.23, 4.56};
    
    // ret
    A.set_ret(0, val, fft_domain::time);
    cplx result;
    A.get_ret(0, result, fft_domain::time);
    REQUIRE(std::abs(val - result) < 1e-12);
    
    // les
    A.set_les(0, val, fft_domain::time);
    A.get_les(0, result, fft_domain::time);
    REQUIRE(std::abs(val - result) < 1e-12);
}

TEST_CASE("herm_matrix_ness construction scalar", "[ness][herm_matrix_ness][construction]") {
    int Nft = 4;
    int size1 = 1; // scalar case
    herm_matrix_ness G(Nft, size1);
    fft_array Gret(Nft, size1);
    fft_array Gles(Nft, size1);
    
    // RET
    // deterministic values time
    std::vector<std::complex<double>> M(Nft);
    for(int t=0; t<Nft; t++){
        M[t] = std::complex<double>(1.0 * t, - 2.0 * t);
        Gret.set_element(t, M[t], fft_domain::time);
    }
    for(int t=0; t<Nft; t++){
        G.set_ret(t, M[t], fft_domain::time);
    }
    REQUIRE(distance_norm2(G.ret_, Gret, fft_domain::time) < 1e-12);
    // deterministic values freq
    for(int t=0; t<Nft; t++){
        M[t] = std::complex<double>(- 1.0 * t, 2.0 * t);
        Gret.set_element(t, M[t], fft_domain::freq);
    }
    for(int t=0; t<Nft; t++){
        G.set_ret(t, M[t], fft_domain::freq);
    }
    REQUIRE(distance_norm2(G.ret_, Gret, fft_domain::freq) < 1e-12);
    
    // LES
    // deterministic values time
    for(int t=0; t<Nft; t++){
        M[t] = std::complex<double>(1.0 * t, - 2.0 * t);
        Gles.set_element(t, M[t], fft_domain::time);
    }
    for(int t=0; t<Nft; t++){
        G.set_les(t, M[t], fft_domain::time);
    }
    REQUIRE(distance_norm2(G.les_, Gles, fft_domain::time) < 1e-12);
    // deterministic values freq
    for(int t=0; t<Nft; t++){
        M[t] = std::complex<double>(- 1.0 * t, 2.0 * t);
        Gles.set_element(t, M[t], fft_domain::freq);
    }
    for(int t=0; t<Nft; t++){
        G.set_les(t, M[t], fft_domain::freq);
    }
    REQUIRE(distance_norm2(G.les_, Gles, fft_domain::freq) < 1e-12);
}

TEST_CASE("herm_matrix_ness construction matrix", "[ness][herm_matrix_ness][construction]") {
    int Nft = 4;
    int size1 = 2; // matrix case
    herm_matrix_ness G(Nft, size1);
    fft_array Gret(Nft, size1);
    fft_array Gles(Nft, size1);
    
    // RET
    // deterministic values time
    std::vector<cdmatrix> M(Nft, cdmatrix(size1,size1));
    for(int t=0; t<Nft; t++){
        for(int i=0; i<size1; i++){
            for(int j=0; j<size1; j++){
                M[t](i,j) = std::complex<double>(1.0 * t + i, - 2.0 * t + j);
            }
        }
        Gret.set_element(t, M[t], fft_domain::time);
    }
    for(int t=0; t<Nft; t++){
        G.set_ret(t, M[t], fft_domain::time);
    }
    REQUIRE(distance_norm2(G.ret_, Gret, fft_domain::time) < 1e-12);
    // deterministic values freq
    for(int t=0; t<Nft; t++){
        for(int i=0; i<size1; i++){
            for(int j=0; j<size1; j++){
                M[t](i,j) = std::complex<double>(- 1.0 * t + i, 2.0 * t + j);
            }
        }
        Gret.set_element(t, M[t], fft_domain::freq);
    }
    for(int t=0; t<Nft; t++){
        G.set_ret(t, M[t], fft_domain::freq);
    }
    REQUIRE(distance_norm2(G.ret_, Gret, fft_domain::freq) < 1e-12);
    
    // LES
    // deterministic values time
    for(int t=0; t<Nft; t++){
        for(int i=0; i<size1; i++){
            for(int j=0; j<size1; j++){
                M[t](i,j) = std::complex<double>(1.0 * t + i, - 2.0 * t + j);
            }
        }
        Gles.set_element(t, M[t], fft_domain::time);
    }
    for(int t=0; t<Nft; t++){
        G.set_les(t, M[t], fft_domain::time);
    }
    REQUIRE(distance_norm2(G.les_, Gles, fft_domain::time) < 1e-12);
    // deterministic values freq
    for(int t=0; t<Nft; t++){
        for(int i=0; i<size1; i++){
            for(int j=0; j<size1; j++){
                M[t](i,j) = std::complex<double>(- 1.0 * t + i, 2.0 * t + j);
            }
        }
        Gles.set_element(t, M[t], fft_domain::freq);
    }
    for(int t=0; t<Nft; t++){
        G.set_les(t, M[t], fft_domain::freq);
    }
    REQUIRE(distance_norm2(G.les_, Gles, fft_domain::freq) < 1e-12);
}

