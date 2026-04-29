#include <iostream>
#include <fstream>
#include <complex>
#include <string>

#include "cntr/ness2.hpp"
#include "cntr/utils/read_inputfile.hpp"

using namespace ness2;

class box_dos {
public:
    double hi_, lo_, wc_, nu_;
    box_dos() : wc_(10.0), nu_(3.0), hi_(15.0), lo_(-15.0) {}
    double operator()(double w) {
        return 1.0 / ((1 + std::exp(nu_ * (w - wc_))) * (1 + std::exp(-nu_ * (wc_ + w))))/M_PI;
    }
};

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: ness_test_siam.x inputfile.txt outputfile.h5\n";
        return 1;
    }

    // === Read input parameters ===
    std::string outfile = argv[2];

    double U, beta, V, epsd, h, eta;
    int Nft;
    find_param(argv[1], "__U=", U);
    find_param(argv[1], "__beta=", beta);
    find_param(argv[1], "__V=", V);
    find_param(argv[1], "__epsd=", epsd);
    find_param(argv[1], "__Nft=", Nft);
    find_param(argv[1], "__h=", h);
    find_param(argv[1], "__eta=", eta);

    int size = 1;

    // === Allocate Green's functions ===
    herm_matrix_ness G0(Nft, size);
    herm_matrix_ness G(Nft, size);
    herm_matrix_ness SL(Nft, size);
    herm_matrix_ness SR(Nft, size);
    herm_matrix_ness Sbath(Nft, size);
    herm_matrix_ness SU(Nft, size);
    herm_matrix_ness S(Nft, size);
    herm_matrix_ness W(Nft, size);  // auxiliary

    // === Set up bath self-energies ===
    box_dos dos;
    green_equilibrium_ness(FERMION,SL, dos, beta, +0.5 * V, h, FFT_TRAPEZ);
    green_equilibrium_ness(FERMION,SR, dos, beta, -0.5 * V, h, FFT_TRAPEZ);
    Sbath = SL;
    Sbath.incr(SR,1.0, fft_domain::time);

    // === Solve noninteracting Dyson equation ===
    cdmatrix eps_matrix(1, 1);
    eps_matrix(0, 0) = epsd;
    double mu = 0.0;
    dyson(G0, mu, eps_matrix, Sbath, h, FFT_TRAPEZ, BATH_GAUSS, eta, beta, mu);

    // === Compute 2nd order diagram ===
    cdmatrix U_matrix(1, 1);
    U_matrix(0, 0) = U;
    Bubble1_ness(W, G0, G0);  // W = ii * G * G
    W.smul(U*U,fft_domain::time);
    //W.left_multiply(U_matrix, fft_domain::time);
    //W.right_multiply(U_matrix, fft_domain::time);
    Bubble2_ness(SU, G0, W);  // SU = -ii * G * W
    SU.smul(-1.0, fft_domain::time);

    // === Solve interacting Dyson ===
    S = Sbath;
    S.incr(SU,1.0, fft_domain::time);
    dyson(G, mu, eps_matrix, S, h, FFT_TRAPEZ, BATH_GAUSS, eta, beta, mu);

    // make a round-trip test, to compare
    //herm_matrix_ness G1=G;
    //G1.integral_transform_to_freq(h,FFT_TRAPEZ);
    //double d= distance_norm2(G.ret_,G1.ret_,fft_domain::freq);
    //std::cout << "Roundtrip G ->  transform_to_time  ->  transform_to_freq = S1 :  |G-G1|=  " << d << std::endl;
    //G.print_to_file("G.out");
    //G1.print_to_file("G1.out");

    
    

    
    
    // === Observables ===
    cdmatrix SUG,SLG,GSL,rho;
    // ii*[SL*G]^<(t=0):
    ness2::convolution_density_matrix(SLG,FERMION,SL,G,h); 
    ness2::convolution_density_matrix(SLG,FERMION,G,SL,h); 
    // ii*[SU*G]^<(t=0):
    ness2::convolution_density_matrix(SUG,FERMION,SU,G,h); 
    ness2::density_matrix(rho,FERMION,G); //   = ii*G^<(t=0)
    double Current = 4.0*SLG.trace().imag(); 
    double Eint = (SUG).trace().real(); // interaction energy
    double dens = 2*rho.trace().real(); // <n_up + n_do >
    // check:
    double Current1 = -4.0*GSL.trace().imag();
    std::cout << "test hermiticity pof current: " << Current1-Current << std::endl;

    // === Fourier transform for spectra ===
    G.integral_transform_to_freq(h, FFT_TRAPEZ);
    G0.integral_transform_to_freq(h, FFT_TRAPEZ);
    SL.integral_transform_to_freq(h, FFT_TRAPEZ);
    SR.integral_transform_to_freq(h, FFT_TRAPEZ);
    SU.integral_transform_to_freq(h, FFT_TRAPEZ);

    // === Write results to HDF5 ===
    hid_t file_id = open_hdf5_file(outfile.c_str());

    store_int_attribute_to_hid(file_id, "Nft", Nft);
    store_double_attribute_to_hid(file_id, "dens", dens);
    store_double_attribute_to_hid(file_id, "Eint", Eint);
    store_double_attribute_to_hid(file_id, "Current", Current);

    G0.write_to_hdf5(file_id, "G0");
    G.write_to_hdf5(file_id, "G");
    SL.write_to_hdf5(file_id, "SL");
    SR.write_to_hdf5(file_id, "SR");
    SU.write_to_hdf5(file_id, "SU");
    S.write_to_hdf5(file_id, "S");

    close_hdf5_file(file_id);

    std::cout << "Simulation complete. Data written to " << outfile << "\n";
    return 0;
}


