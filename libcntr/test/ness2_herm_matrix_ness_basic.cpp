#include "catch.hpp"
#include "ness2.hpp"

using namespace ness2;


TEST_CASE("green_equilibrium_ness_fermi matches CNTR Green's function via read_from_cntr", "[ness][herm_matrix_ness][to_cntr]") {
    int Nft = 64;
    int size1 = 1;
    double beta = 10.0;
    double mu = 0.0;
    double h = 0.1;

    gauss dos(0.0, 2.0);

    // === CNTR version ===
    cntr::herm_matrix<double> G_cntr(Nft/2,10, size1, FERMION);
    green_equilibrium(G_cntr, dos, beta, h, mu, 100, 20);

    // === NESS from CNTR ===
    ness2::herm_matrix_ness G_from_cntr(Nft, size1); // 
    G_from_cntr.read_from_cntr(G_cntr,Nft/2-1);

    // === NESS computed directly ===
    ness2::herm_matrix_ness G_direct(Nft, size1);
    green_equilibrium_ness(FERMION,G_direct, dos, beta, mu, h, FFT_ADAPTIVE);

    // === Compare ret_t and les_t ===
    for (int t = 0; t < Nft; ++t) {
        cdmatrix A, B;
        G_from_cntr.get_ret(t, A, fft_domain::time);
        G_direct.get_ret(t, B, fft_domain::time);
	//std::cout << A << " " << B << std::endl; 
        REQUIRE((A - B).norm() < 1e-12);
        G_from_cntr.get_les(t, A, fft_domain::time);
        G_direct.get_les(t, B, fft_domain::time);
	///std::cout << A << " " << B << std::endl; 
        REQUIRE((A - B).norm() < 1e-12);
    }
}


TEST_CASE("Compare Bubble1 between CNTR and NESS", "[ness][herm_matrix_ness][bubble1]") {
    int Nft = 32;
    int ntau=20,nt=Nft/2-1;
    double beta = 10.0;
    double h = 0.1;
    double mu = 0.0;

    int sizeA = 2;
    int sizeB = 1;


    // === Construct CNTR Green's functions ===
    cntr::herm_matrix<double> A(nt,ntau, sizeA,FERMION);
    cntr::herm_matrix<double> B(nt,ntau, sizeB,FERMION);
    cntr::herm_matrix<double> C(nt,ntau, sizeA,FERMION);

    

    // === Fill A and B manually using set_ret / set_les ===
    for (int t = 0; t <= nt; ++t) {
      for (int t1 = 0; t1 <= t; ++t1) {
        // === A: sizeA x sizeA ===
        cdmatrix MretA(sizeA, sizeA);
        cdmatrix MlesA(sizeA, sizeA);
        for (int i = 0; i < sizeA; ++i) {
	  for (int j = 0; j < sizeA; ++j) {
	    cdouble val_ret(0.1 * t + 0.01 * t1*t1 + 0.01 * i, -0.02 * j);
	    cdouble val_les(0.05 * t1 + 0.02 * i+ + 0.01 * t*t, 0.03 * j+0.005*t*t1);
	    MretA(i, j) = val_ret;
	    MlesA(i, j) = val_les;
	  }
        }
	if(t==t1){
	  cdmatrix tmp=MlesA-MlesA.adjoint();
	  MlesA=tmp;
	  tmp=MretA-MretA.adjoint();
	  MretA=tmp; 
	}
	A.set_ret(t,t1, MretA);  // only for t ≥ t1
	A.set_les(t1,t, MlesA);  // only for t ≤ t1
	
        // === B: sizeB x sizeB ===
        cdmatrix MretB(sizeB, sizeB);
        cdmatrix MlesB(sizeB, sizeB);
        MretB(0, 0) = cdouble(0.2 * t + 0.05 * t1, -0.1 * t1+0.03);
        MlesB(0, 0) = cdouble(0.1 * t1 + 0.03 * t, 0.1 * t+0.02);
	if(t==t1){
	  MretB=MretB-MretB.adjoint();
	  MlesB=MlesB-MlesB.adjoint();
	}
	B.set_ret(t,t1, MretB);
	B.set_les(t1,t, MlesB);
      }
    }
    // === Direct CNTR → NESS conversion (no file I/O) ===
    herm_matrix_ness Aness(Nft, sizeA);
    herm_matrix_ness Bness(Nft, sizeB);
    herm_matrix_ness Cness(Nft, sizeA);  // same size as C
    herm_matrix_ness Cness2(Nft, sizeA);

    
    // === Run CNTR Bubble1 ===
    for (int i = 0; i < sizeA; ++i) {
        for (int j = 0; j < sizeA; ++j) {
	  cntr::Bubble1(nt,C, i, j, A, i, j, B, 0, 0);
        }
    }    
    Aness.read_from_cntr(A, nt);
    Bness.read_from_cntr(B, nt);
    Cness.read_from_cntr(C, nt);
        // === Recompute Bubble1 in NESS ===
    for (int i = 0; i < sizeA; ++i) {
        for (int j = 0; j < sizeA; ++j) {
            Bubble1_ness(Cness2, i, j, Aness, i, j, Bness, 0, 0);
        }
    }
    // === Compare results ===
    double d = distance_norm2(Cness, Cness2, fft_domain::time);
    REQUIRE(d < 1e-10); // Tight tolerance

    // === Run CNTR Bubble1 ===
    for (int i = 0; i < sizeA; ++i) {
        for (int j = 0; j < sizeA; ++j) {
	  cntr::Bubble1(nt,C, i, j, B, 0, 0, A, i, j);
        }
    }    
    Aness.read_from_cntr(A, nt);
    Bness.read_from_cntr(B, nt);
    Cness.read_from_cntr(C, nt);
        // === Recompute Bubble1 in NESS ===
    for (int i = 0; i < sizeA; ++i) {
        for (int j = 0; j < sizeA; ++j) {
            Bubble1_ness(Cness2, i, j, Bness, 0, 0, Aness, i, j);
        }
    }
    // === Compare results ===
    d = distance_norm2(Cness, Cness2, fft_domain::time);
    REQUIRE(d < 1e-10); // Tight tolerance

    
    // === Run CNTR Bubble2 ===
    for (int i = 0; i < sizeA; ++i) {
        for (int j = 0; j < sizeA; ++j) {
	  cntr::Bubble2(nt,C, i, j, A, i, j, B, 0, 0);
        }
    }    
    Aness.read_from_cntr(A, nt);
    Bness.read_from_cntr(B, nt);
    Cness.read_from_cntr(C, nt);
        // === Recompute Bubble1 in NESS ===
    for (int i = 0; i < sizeA; ++i) {
        for (int j = 0; j < sizeA; ++j) {
            Bubble2_ness(Cness2, i, j, Aness, i, j, Bness, 0, 0);
        }
    }
    // === Compare results ===
    d = distance_norm2(Cness, Cness2, fft_domain::time);
    REQUIRE(d < 1e-10); // Tight tolerance

    // === Run CNTR Bubble2 ===
    for (int i = 0; i < sizeA; ++i) {
        for (int j = 0; j < sizeA; ++j) {
	  cntr::Bubble2(nt,C, i, j, B, 0, 0, A, i, j);
        }
    }    
    Aness.read_from_cntr(A, nt);
    Bness.read_from_cntr(B, nt);
    Cness.read_from_cntr(C, nt);
    // === Recompute Bubble1 in NESS ===
    for (int i = 0; i < sizeA; ++i) {
        for (int j = 0; j < sizeA; ++j) {
            Bubble2_ness(Cness2, i, j, Bness, 0, 0, Aness, i, j);
        }
    }
    // === Compare results ===
    d = distance_norm2(Cness, Cness2, fft_domain::time);
    REQUIRE(d < 1e-10); // Tight tolerance

}


TEST_CASE("NESS <--> CNTR round-trip with distance norm", "[ness][herm_matrix_ness][ness2cntr]") {
    int Nft = 32;
    int size = 2;
    int nt = Nft / 2 - 1;
    int ntau=20;
    int tc = nt;

    // === Step 1: Initialize Gness with deterministic data ===
    ness2::herm_matrix_ness Gness(Nft, size);
    for (int t = 0; t <= nt; ++t) {
        cdmatrix Mret(size, size), Mles(size, size);
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                Mret(i, j) = cdouble(0.01 * t + 0.001 * i, -0.002 * j+0.1*t);
                Mles(i, j) = cdouble(0.002 * t + 0.001 * i, 0.001 * j+0.1*t*t);
            }
        }
        if (t == 0) {
            cdmatrix tmp = Mret - Mret.adjoint();
	    Mret-=tmp;
            tmp = Mles - Mles.adjoint();
	    Mles=tmp;
        }
        Gness.set_ret(t, Mret, fft_domain::time);
        Gness.set_les(t, Mles, fft_domain::time);
        if (t > 0) {
            Gness.set_les(-t, -Mles.adjoint(), fft_domain::time);
        }
    }

    // === Step 2: Convert Gness to Gcntr and Gmoving ===
    cntr::herm_matrix<double> Gcntr(nt, ntau, size, FERMION);
    ness2::ness2cntr(Gcntr, Gness);

    cntr::herm_matrix_moving<double> Gmoving(tc, size, FERMION);
    ness2::ness2cntr(Gmoving, Gness);

    // === Step 3: Generate GfromGcntr from full Gcntr using backward mapping ===
    cntr::herm_matrix_moving<double> GfromGcntr(tc, size, FERMION);
    GfromGcntr.set_from_G_backward(Gcntr, Gcntr, tc); // Acc = Gcntr because hermitian

    // === Step 4: Compare Gmoving vs GfromGcntr using distance_norm2 on last step ===
    double d = distance_norm2(Gmoving, 0, GfromGcntr, 0);
    //Gmoving.print_to_file("A1.txt");
    //GfromGcntr.print_to_file("A2.txt");    
    REQUIRE(d < 1e-12);

    // === Step 5: Back-convert Gcntr -> Gness_from_cntr ===
    ness2::herm_matrix_ness Gness_from_cntr(Nft, size);
    ness2::cntr2ness(Gness_from_cntr, Gcntr, nt);

    // === Step 6: Back-convert Gmoving → Gness_from_moving ===
    ness2::herm_matrix_ness Gness_from_moving(Nft, size);
    ness2::cntr2ness(Gness_from_moving, Gmoving, 0);

    // === Step 7: Compare both back-converted GFs to original ===
    double d1 = ness2::distance_norm2(Gness, Gness_from_cntr, fft_domain::time);
    double d2 = ness2::distance_norm2(Gness, Gness_from_moving, fft_domain::time);

    REQUIRE(d1 < 1e-12);
    REQUIRE(d2 < 1e-12);
}

TEST_CASE("force_equilibrium matches green_equilibrium_ness", "[ness][herm_matrix_ness][force_equilibrium]") {
    int Nft = 64;
    int size1 = 1;
    double beta = 10.0;
    double mu = 0.0;
    double h = 0.1;

    gauss dos(0.0, 2.0);

    // reference G from dos
    ness2::herm_matrix_ness G_ref(Nft, size1);
    green_equilibrium_ness(FERMION, G_ref, dos, beta, mu, h, FFT_ADAPTIVE);

    // set only ret from reference G, and then force les
    ness2::herm_matrix_ness G_eq(Nft, size1);
    for(int t=0;t<Nft;t++){
        cdmatrix m(size1,size1);
        G_ref.get_ret(t,m,fft_domain::time);
        G_eq.set_ret(t,m,fft_domain::freq);
    }
    G_eq.force_equilibrium(FERMION, beta, mu, h);

    REQUIRE((ness2::distance_norm2(G_eq,G_ref,fft_domain::freq)/Nft)<1e-10);
}

#if CNTR_USE_HDF5 == 1

#include "hdf5.h"

TEST_CASE("herm_matrix_ness HDF5 round-trip with nontrivial data, group_id", "[ness][herm_matrix_ness][hdf5]") {
    int Nft = 4;
    int size1 = 2;
    herm_matrix_ness A(Nft, size1);

    // Fill with nontrivial values
    for (int t = 0; t < Nft; ++t) {
        cdmatrix M(size1, size1);
        for (int i = 0; i < size1; ++i) {
            for (int j = 0; j < size1; ++j) {
                M(i, j) = cplx(i + j + 0.5, t * 0.1);
            }
        }
        A.set_ret(t, M, fft_domain::time);
        A.set_les(t, M, fft_domain::time);
        M *= 0.2;
        A.set_ret(t, M, fft_domain::freq);
        A.set_les(t, M, fft_domain::freq);
    }

    const std::string fname = "herm_matrix_ness_test.h5";
    hid_t file_id = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    REQUIRE(file_id >= 0);

    hid_t group_id = H5Gcreate(file_id, "/herm_matrix_ness_data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    REQUIRE(group_id >= 0);

    A.write_to_hdf5(group_id);

    H5Gclose(group_id);
    H5Fclose(file_id);

    // Read into a new object
    herm_matrix_ness B(1, 1); // dummy, will be resized
    file_id = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    REQUIRE(file_id >= 0);

    group_id = H5Gopen(file_id, "/herm_matrix_ness_data", H5P_DEFAULT);
    REQUIRE(group_id >= 0);

    B.read_from_hdf5(group_id);

    H5Gclose(group_id);
    H5Fclose(file_id);

    // Compare
    for (int t = 0; t < Nft; ++t) {
        cdmatrix A_ret, B_ret, A_les, B_les;
        A.get_ret(t, A_ret, fft_domain::time);
        B.get_ret(t, B_ret, fft_domain::time);
        A.get_les(t, A_les, fft_domain::time);
        B.get_les(t, B_les, fft_domain::time);

        REQUIRE((A_ret - B_ret).norm() < 1e-12);
        REQUIRE((A_les - B_les).norm() < 1e-12);

        A.get_ret(t, A_ret, fft_domain::freq);
        B.get_ret(t, B_ret, fft_domain::freq);
        A.get_les(t, A_les, fft_domain::freq);
        B.get_les(t, B_les, fft_domain::freq);

        REQUIRE((A_ret - B_ret).norm() < 1e-12);
        REQUIRE((A_les - B_les).norm() < 1e-12);
    }
}

TEST_CASE("herm_matrix_ness HDF5 round-trip with nontrivial data, file_id", "[ness][herm_matrix_ness][hdf5]") {
    int Nft = 4;
    int size1 = 2;
    herm_matrix_ness A(Nft, size1);

    // Fill with nontrivial values
    for (int t = 0; t < Nft; ++t) {
        cdmatrix M(size1, size1);
        for (int i = 0; i < size1; ++i) {
            for (int j = 0; j < size1; ++j) {
                M(i, j) = cplx(i + j + 0.5, t * 0.1);
            }
        }
        A.set_ret(t, M, fft_domain::time);
        A.set_les(t, M, fft_domain::time);
        M *= 0.2;
        A.set_ret(t, M, fft_domain::freq);
        A.set_les(t, M, fft_domain::freq);
    }

    const std::string fname = "herm_matrix_ness_test2.h5";
    const std::string gname = "/herm_matrix_ness_data2";
    hid_t file_id = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    REQUIRE(file_id >= 0);

    A.write_to_hdf5(file_id, gname.c_str());

    H5Fclose(file_id);

    // Read into a new object
    herm_matrix_ness B(1, 1); // dummy, will be resized
    file_id = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    REQUIRE(file_id >= 0);

    B.read_from_hdf5(file_id, gname.c_str());

    H5Fclose(file_id);

    // Compare
    for (int t = 0; t < Nft; ++t) {
        cdmatrix A_ret, B_ret, A_les, B_les;
        A.get_ret(t, A_ret, fft_domain::time);
        B.get_ret(t, B_ret, fft_domain::time);
        A.get_les(t, A_les, fft_domain::time);
        B.get_les(t, B_les, fft_domain::time);

        REQUIRE((A_ret - B_ret).norm() < 1e-12);
        REQUIRE((A_les - B_les).norm() < 1e-12);

        A.get_ret(t, A_ret, fft_domain::freq);
        B.get_ret(t, B_ret, fft_domain::freq);
        A.get_les(t, A_les, fft_domain::freq);
        B.get_les(t, B_les, fft_domain::freq);

        REQUIRE((A_ret - B_ret).norm() < 1e-12);
        REQUIRE((A_les - B_les).norm() < 1e-12);
    }
}

TEST_CASE("herm_matrix_ness HDF5 round-trip with nontrivial data, filename", "[ness][herm_matrix_ness][hdf5]") {
    int Nft = 4;
    int size1 = 2;
    herm_matrix_ness A(Nft, size1);

    // Fill with nontrivial values
    for (int t = 0; t < Nft; ++t) {
        cdmatrix M(size1, size1);
        for (int i = 0; i < size1; ++i) {
            for (int j = 0; j < size1; ++j) {
                M(i, j) = cplx(i + j + 0.5, t * 0.1);
            }
        }
        A.set_ret(t, M, fft_domain::time);
        A.set_les(t, M, fft_domain::time);
        M *= 0.2;
        A.set_ret(t, M, fft_domain::freq);
        A.set_les(t, M, fft_domain::freq);
    }

    const std::string fname = "herm_matrix_ness_test3.h5";
    const std::string gname = "/herm_matrix_ness_data3";
    hid_t file_id = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    REQUIRE(file_id >= 0);

    A.write_to_hdf5(fname.c_str(), gname.c_str());

    H5Fclose(file_id);
    // Read into a new object
    herm_matrix_ness B(1, 1); // dummy, will be resized
    file_id = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    REQUIRE(file_id >= 0);
    
    B.read_from_hdf5(fname.c_str(), gname.c_str());

    H5Fclose(file_id);
    // Compare
    for (int t = 0; t < Nft; ++t) {
        cdmatrix A_ret, B_ret, A_les, B_les;
        A.get_ret(t, A_ret, fft_domain::time);
        B.get_ret(t, B_ret, fft_domain::time);
        A.get_les(t, A_les, fft_domain::time);
        B.get_les(t, B_les, fft_domain::time);

        REQUIRE((A_ret - B_ret).norm() < 1e-12);
        REQUIRE((A_les - B_les).norm() < 1e-12);

        A.get_ret(t, A_ret, fft_domain::freq);
        B.get_ret(t, B_ret, fft_domain::freq);
        A.get_les(t, A_les, fft_domain::freq);
        B.get_les(t, B_les, fft_domain::freq);

        REQUIRE((A_ret - B_ret).norm() < 1e-12);
        REQUIRE((A_les - B_les).norm() < 1e-12);
    }
}


#endif
