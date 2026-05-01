//#include "ness_global_settings.hpp"
//#include "ness_cpp.hpp"
#include "ness.hpp"

#include <fstream>


namespace ness
{

/** \brief <b> Set the lesser component of a `GF` from the retarded one given a distribution function \f$ F(\omega)\f$
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 * > Set the lesser component of a `GF` from its retarded component and a given real distribution function \f$ F(\omega)\f$
 * > using the Fluctuation-Dissipation Theorem in frequency space:
 * > \f$ G^<(\omega) = F(\omega)G^A(\omega) - G^R(\omega)F(\omega)\f$.
 * > This is for one-orbital GFs.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param gf
 *  [GF] Frequency GF object to be set
 * @param fd
 *  [std::vector<double>] Real distribution function
 */
void set_les_from_ret(GF &gf, std::vector<double> fd)
{
    assert(gf.ngrid_ == fd.size());
    assert(gf.size1_ == 1);
    assert(gf.size2_ == 1);
    assert(gf.gf_type_ == freq_gf);
    for(int w = 0; w < gf.ngrid_; w++)
    {
        gf.Lesser[w] = fd[w] * gf.Retarded[w].adjoint() - gf.Retarded[w] * fd[w];
    }
}

/** \brief <b> Set the lesser component of a `GF_pair` from the retarded one given a distribution function \f$ F(\omega)\f$
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 * > Set the lesser component of a `GF_pair` from its retarded component and a given real distribution function \f$ F(\omega)\f$
 * > using the Fluctuation-Dissipation Theorem in frequency space:
 * > \f$ G^<(\omega) = F(\omega)G^A(\omega) - G^R(\omega)F(\omega)\f$.
 * > This is for one-orbital GFs.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param gf
 *  [GF_pair] GF pair object to be set
 * @param fd
 *  [std::vector<double>] Real distribution function.
 */
void set_les_from_ret(GF_pair &gf, std::vector<double> fd)
{
    assert(gf.freq.ngrid_ == fd.size());
    assert(gf.freq.size1_ == 1);
    assert(gf.freq.size2_ == 1);
    gf.to_freq();
    for(int w = 0; w < gf.freq.ngrid_; w++)
    {
        gf.freq.Lesser[w] = fd[w] * gf.freq.Retarded[w].adjoint() - gf.freq.Retarded[w] * fd[w];
    }
    gf.to_time();
}

/** \brief <b> Calculate  the FFT frequency grid \f$ d\omega, N_{\omega} \f$ from a given time grid \f$ h, N_t \f$</b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 * > Calculate  the FFT frequency grid \f$ d\omega, N_{\omega}\f$ from a given time grid \f$dt, N_t\f$,
 * > so that it matches the FFT requirements:
 * > \f$ N_{\omega} = \frac{2 (N_t -1)} {3}  \f$,
 * >  \f$ d\omega = \frac{\pi }{(N_t-1)  dt} \f$.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param df
 *  [double] frequency grid spacing
 * @param nf
 *  [int] Number of frequency grid points
 * @param dt
 *  [double] Given time grid spacing
 * @param nt
 *  [int] Given number of time grid points
 */
void fgrid_from_tgrid(double & df, int & nf, double & dt, int & nt)
{
    nf = 2 * (nt -1) / 3;
    df = M_PI / ((nt-1) * dt);
}

/** \brief <b> Calculate  the FFT time grid \f$dt, N_t\f$ from a given frequency grid \f$ d\omega, N_{\omega} \f$</b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 * > Calculate  the FFT time grid \f$h, N_t\f$ from a given frequency grid \f$ d\omega, N_{\omega} \f$,
 * > so that it matches the FFT requirements:
 * > \f$ N_t = \frac{3 N_{\omega}}{2} + 1  \f$,
 * >  \f$ dt = \frac{2 \pi} { 3  N_{\omega}  d\omega} \f$.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param df
 *  [double] frequency grid spacing
 * @param nf
 *  [int] Number of frequency grid points
 * @param dt
 *  [double] Given time grid spacing h
 * @param nt
 *  [int] Given number of time grid points
 */
void tgrid_from_fgrid(double & dt, int & nt, double & df, int & nf)
{
    nt = 3 * nf / 2 + 1;
    dt = 2 * M_PI / (3 * nf *df);
}

/** \brief <b> Calculate the bubble diagram of first type: \f$ C_{c_1,c_2} (t)= i A_{a_1,a_2} (t) B_{b_2,b_1} (-t) \f$ </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Evaluate `GF` Green's function \f$ C \f$ given by the bubble diagram equation
 * > \f$ C_{c_1,c_2} (t)= i A_{a_1,a_2} (t) B_{b_2,b_1} (-t) \f$.
 * > The evaluation is done for all time points and for the lesser ad retarded component.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param C
 *  [GF] Time Green's function to store result in
 * @param c1
 *  First orbital index of C
 * @param c2
 *  Second orbital index of C
 * @param A
 *  [GF] Time Green's function
 * @param a1
 *  First orbital index of A
 * @param a2
 *  Second orbital index of A
 * @param B
 *  [GF] Time Green's function
 * @param b1
 *  First orbital index of B
 * @param b2
 *  Second orbital index of B
 */
void Bubble1_ness(GF &C,int c1,int c2,GF &A,int a1,int a2,GF &B,int b1,int b2)
{
    assert(C.ngrid_ == A.ngrid_);
    assert(C.ngrid_ == B.ngrid_);
    assert(C.gf_type_ == time_gf);
    assert(B.gf_type_ == time_gf);
    assert(A.gf_type_ == time_gf);
    assert(a1 <= A.size1_);
    assert(a2 <= A.size1_);
    assert(b1 <= B.size1_);
    assert(b2 <= B.size1_);
    assert(c1 <= C.size1_);
    assert(c2 <= C.size1_);
    
    for(auto w = 0; w < C.ngrid_; w++){
        C.Lesser[w](c1, c2) = NESS_II * A.Lesser[w](a1,a2) * (-B.Greater(w).adjoint()(b1,b2));
        C.Retarded[w](c1, c2) = NESS_II * A.Retarded[w](a1,a2) * (-B.Lesser[w].adjoint()(b1,b2)) + NESS_II * A.Lesser[w](a1,a2) * B.Retarded[w].adjoint()(b1,b2);
    }
}

/** \brief <b> Calculate the bubble diagram of first type: \f$ C_{c_1,c_2} (t)= i A_{a_1,a_2} (t) B_{b_2,b_1} (-t) \f$ </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Evaluate `GF_pair` Green's function \f$ C \f$ given by the bubble diagram equation
 * > \f$ C_{c_1,c_2} (t)= i A_{a_1,a_2} (t) B_{b_2,b_1} (-t) \f$.
 * > The evaluation is done for all time points and for the lesser ad retarded component.
 * > This is the version of the function which uses the time member of Green's functions pairs
 * > to calculate the bubble.
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param C
 * [GF_pair] Time Green's function pair, store result in time member
 * @param c1
 *  First orbital index of C
 * @param c2
 *  Second orbital index of C
 * @param A
 *   [GF_pair] Green's function pair
 * @param a1
 *  First orbital index of A
 * @param a2
 *  Second orbital index of A
 * @param B
 *  [GF_pair] Green's function pair
 * @param b1
 *  First orbital index of B
 * @param b2
 *  Second orbital index of B
 */
void Bubble1_ness(GF_pair &C,int c1,int c2,GF_pair &A,int a1,int a2,GF_pair &B,int b1,int b2)
{
    assert(C.time.ngrid_ == A.time.ngrid_);
    assert(C.time.ngrid_ == B.time.ngrid_);
    assert(a1 <= A.time.size1_);
    assert(a2 <= A.time.size1_);
    assert(b1 <= B.time.size1_);
    assert(b2 <= B.time.size1_);
    assert(c1 <= C.time.size1_);
    assert(c2 <= C.time.size1_);
    
    for(auto w = 0; w < C.time.ngrid_; w++){
        C.time.Lesser[w](c1, c2) = NESS_II * A.time.Lesser[w](a1,a2) * (-B.time.Greater(w).adjoint()(b1,b2));
        C.time.Retarded[w](c1, c2) = NESS_II * A.time.Retarded[w](a1,a2) * (-B.time.Lesser[w].adjoint()(b1,b2)) + NESS_II * A.time.Lesser[w](a1,a2) * B.time.Retarded[w].adjoint()(b1,b2);
    }
}

///@private
void Bubble1_ness(GF &C, GF &A, GF &B)
{
    return Bubble1_ness(C,0,0,A,0,0,B,0,0);
}

///@private
void Bubble1_ness(GF_pair &C, GF_pair &A, GF_pair &B)
{
    return Bubble1_ness(C,0,0,A,0,0,B,0,0);
}

/** \brief <b> Calculate the bubble diagram of second type:  \f$ C_{c_1,c_2} (t)= i A_{a_1,a_2} (t) B_{b_1,b_2} (t) \f$ </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Evaluate `GF` Green's function \f$ C\f$ given by the bubble diagram equation
 * > \f$ C_{c_1,c_2} (t)= i A_{a_1,a_2} (t) B_{b_1,b_2} (t) \f$.
 * > The evaluation is done for all time points and for the lesser ad retarded component.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param C
 *  [GF] Time Green's function to store result in
 * @param c1
 *  First orbital index of C
 * @param c2
 *  Second orbital index of C
 * @param A
 *   [GF] Time Green's function
 * @param a1
 *  First orbital index of A
 * @param a2
 *  Second orbital index of A
 * @param B
 *  [GF]  Time Green's function
 * @param b1
 *  First orbital index of B
 * @param b2
 *  Second orbital index of B
 */
void Bubble2_ness(GF &C,int c1,int c2,GF &A,int a1,int a2,GF &B,int b1,int b2)
{
    assert(C.ngrid_ == A.ngrid_);
    assert(C.ngrid_ == B.ngrid_);
    assert(C.gf_type_ == time_gf);
    assert(B.gf_type_ == time_gf);
    assert(A.gf_type_ == time_gf);
    assert(a1 <= A.size1_);
    assert(a2 <= A.size1_);
    assert(b1 <= B.size1_);
    assert(b2 <= B.size1_);
    assert(c1 <= C.size1_);
    assert(c2 <= C.size1_);
    
    for(auto w = 0; w < C.ngrid_; w++){
        C.Lesser[w](c1, c2) = NESS_II * A.Lesser[w](a1,a2) * B.Lesser[w](b1,b2);
        C.Retarded[w](c1, c2) = NESS_II * A.Retarded[w](a1,a2) * B.Retarded[w](b1,b2) + NESS_II * A.Lesser[w](a1,a2) * B.Retarded[w](b1,b2) + NESS_II * A.Retarded[w](a1,a2) * B.Lesser[w](b1,b2);
    }
}


/** \brief <b> Calculate the bubble diagram of second type:  \f$ C_{c_1,c_2} (t)= i A_{a_1,a_2} (t) B_{b_1,b_2} (t) \f$ </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Evaluate `GF_pair` Green's function \f$ C\f$ given by the bubble diagram equation
 * > \f$ C_{c_1,c_2} (t)= i A_{a_1,a_2} (t) B_{b_1,b_2} (t) \f$.
 * > The evaluation is done for all time points and for the lesser ad retarded component.
 * > This is the version of the function which uses the time member of Green's functions pairs
 * > to calculate the bubble.
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param C
 *  [GF_pair] Time Green's function pair, store result in time member
 * @param c1
 *  First orbital index of C
 * @param c2
 *  Second orbital index of C
 * @param A
 *  [GF_pair] Green's function pair
 * @param a1
 * First orbital index of A
 * @param a2
 *  Second orbital index of A
 * @param B
 *  [GF_pair] Green's function pair
 * @param b1
 *  First orbital index of B
 * @param b2
 *  Second orbital index of B
 */
void Bubble2_ness(GF_pair &C,int c1,int c2,GF_pair &A,int a1,int a2,GF_pair &B,int b1,int b2)
{
    assert(C.time.ngrid_ == A.time.ngrid_);
    assert(C.time.ngrid_ == B.time.ngrid_);
    assert(a1 <= A.time.size1_);
    assert(a2 <= A.time.size1_);
    assert(b1 <= B.time.size1_);
    assert(b2 <= B.time.size1_);
    assert(c1 <= C.time.size1_);
    assert(c2 <= C.time.size1_);
    
    for(auto w = 0; w < C.time.ngrid_; w++){
        C.time.Lesser[w](c1, c2) = NESS_II * A.time.Lesser[w](a1,a2) * B.time.Lesser[w](b1,b2);
        C.time.Retarded[w](c1, c2) = NESS_II * A.time.Retarded[w](a1,a2) * B.time.Retarded[w](b1,b2) + NESS_II * A.time.Lesser[w](a1,a2) * B.time.Retarded[w](b1,b2) + NESS_II * A.time.Retarded[w](a1,a2) * B.time.Lesser[w](b1,b2);
    }
}

///@private
void Bubble2_ness(GF &C, GF &A, GF &B)
{
    return Bubble2_ness(C,0,0,A,0,0,B,0,0);
}

///@private
void Bubble2_ness(GF_pair &C, GF_pair &A, GF_pair &B)
{
    return Bubble2_ness(C,0,0,A,0,0,B,0,0);
}

///@private
/** \brief <b> Impose a symmetry onto the Green's function  </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Impose a symmetry onto the `GF` Green's function \f$ A\f$ where the symmetry is given by a matrix `sym` (for all freqeuncies/times and reatrded and lesser component). \f$ A^{R,<}(t) \to 0.5 ( A^{R,<}(t) + {\tt sym}^\dagger \cdot A^{R,<}(t) \cdot {\tt sym}) \f$
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param gf
 *   [GF] GF to symmetrize
 * @param  sym
 *  [cdmatrix] Matrix representing the symmetry
 */
void impose_sym(GF &gf, cdmatrix & sym)
{
        GF sgf(gf);     
        cdmatrix tsym = sym.adjoint();
        sgf.left_multiply(tsym);
        sgf.right_multiply(sym);
        gf.incr(sgf);
        gf.smul(0.5);
}

///@private
/** \brief <b> Shift a (distribution) function in frequency  </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Shift a (distribution) function given by a vector of doubles by a specified frequency.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param func
 *  [std::vector<double>] (Distribution) function on frequency domain
 * @param dw
 *  [double] Frequency grid spacing
 * @param shift
 *  [double] Frequency shift to be applied
 */
void shift_fd(std::vector<double> & func, double dw, double shift)
{
	int Nw = func.size();
	int dshift = int(shift / dw + 0.00000001);
	if(dshift > 0)
	{
		for(int w = Nw / 2 + 1; w < Nw - dshift; w++)
		{
			func[w] = func[w + dshift];
		}
		for(int w = Nw - dshift; w < Nw; w++)
		{
			func[w] = func[w - Nw + dshift];
		}

		for(int w = 0; w < Nw / 2 + 1 - dshift; w++)
		{
			func[w] = func[w + dshift];
		}
		for(int w = Nw / 2 - dshift + 1; w < Nw / 2 + 1; w++)
		{
			func[w] = 0.0;
		}

	}
	else
	{
		int adshift = -dshift;

		for(int w = Nw / 2; w >= adshift; w--)
		{
			func[w] = func[w - adshift];
		}
		for(int w = adshift - 1; w >= 0; w--)
		{
			func[w] = func[Nw + w - adshift];
		}

		for(int w = Nw - 1; w >= Nw / 2 + 1 + adshift; w--)
		{
			func[w] = func[w - adshift];
		}
		for(int w = Nw / 2 + adshift ; w >= Nw / 2 + 1; w--)
		{
			func[w] = 1.0;
		}

	}
}

/** \brief <b> Set the lesser Green's function (in frequency) by imposing a thermal fermionic distribution  </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > This sets the lesser Green's function (`GF`) by imposing a thermal fermionic distribution of a given temperature `temp` ( \f$ \beta = 1/{\tt temp }\f$ using the Fluctuation-Dissipation theorem.
 * > The Green's function is given in frequency domain.
 * > Constructs \f$ G^<(\omega) \f$ via:
 *   \f[
 *   G^<(\omega) =  \mp 2i\, f(\omega)\, \mathrm{Im}\, G^R(\omega),
 *   \f]
 *   where \f$ f(\omega) = 1 / (e^{\beta \omega} - 1) \f$ is the Fermi-Dirac distribution,
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param gf
 *  [GF] Frequency Green's function
 * @param temp
 *  [double] Temperature of the Fermi function
 */
void force_FD(GF &gf, double temp)
{
	for(int w = 0; w < gf.ngrid_; w++)
	{
		double fd = -tanh(gf.grid_[w] / temp)/2.0 + 0.5;
		gf.Lesser[w] = -(gf.Retarded[w] - gf.Retarded[w].adjoint()) * fd;
	}
}

/** \brief <b> Set the lesser Green's function (in frequency) by imposing a thermal fermionic distribution  </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > This sets the lesser Green's function (`GF_pair`) in frequency domain by imposing a thermal fermionic distribution of a given temperature `temp` ( \f$ \beta = 1/{\tt temp }\f$  using the Fluctuation-Dissipation theorem.
 * > Constructs \f$ G^<(\omega) \f$ via:
 *   \f[
 *   G^<(\omega) =  \mp 2i\, f(\omega)\, \mathrm{Im}\, G^R(\omega),
 *   \f]
 *   where \f$ f(\omega) = 1 / (e^{\beta \omega} - 1) \f$ is the Fermi-Dirac distribution,
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param gf
 *  [GF_pair] Green's function pair
 * @param temp
 *  [double] Temperature of the Fermi function
 */
void force_FD(GF_pair &gf, double temp)
{
    gf.to_freq();
    for(int w = 0; w < gf.freq.ngrid_; w++)
    {
        double fd = -tanh(gf.freq.grid_[w] / temp)/2.0 + 0.5;
        gf.freq.Lesser[w] = -(gf.freq.Retarded[w] - gf.freq.Retarded[w].adjoint()) * fd;
    }
    gf.to_time();
}

///@private
/** \brief <b> This forces the lesser Green's function to be purely imaginary  </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > This forces the lesser Green's function `GF` (in frequenc) to be purely imaginary.
 * > Constructs \f$ G^<(\omega) \f$ via:
 *   \f[
 *   G^<(\omega) \to  0.5(G^<(\omega)-\big(G^<\big)^\dagger(\omega)).
 *   \f]
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param gf
 * [GF] Time Green's function
 */
void force_imag(GF &gf)
{
	if(gf.gf_type_ == time_gf)
		return;
	for(int w = 0; w < gf.ngrid_; w++)
	{
		gf.Lesser[w] = 0.5 * (gf.Lesser[w] - gf.Lesser[w].adjoint());
	}

}

/** \brief <b> This sets a Green's function \f$ G(t)\f$ to the backward \f$ G(-t)\f$ </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > This sets a `GF` Green's function \f$ G(t)\f$ to \f$ G(-t)\f$ in time domain according to the relations:
 * In the lesser component: \f$ G^<(-t) = -[G^>(t)]^\dagger \f$
 * In the retarded component (is zero for \f$ t<0 \f$ but store greater here effectifely) : \f$ G^R(-t) = [G^>(t)]^\dagger - [G^<(t)]^\dagger\f$
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param gf_bwd
 * [GF] Green's function for \f$ -t\f$, result
 * @param gf_fwd
 * [GF] Green's function for \f$ t\f$
 */
int set_bwd_from_fwd(GF & gf_bwd, GF & gf_fwd)
{
	assert(gf_fwd.ngrid_ == gf_bwd.ngrid_);
	assert(gf_fwd.size1_ == gf_bwd.size1_);
	assert(gf_fwd.size2_ == gf_bwd.size2_);

	for(int w = 0; w < gf_fwd.ngrid_; w++)
	{
		gf_bwd.Lesser[w] = -gf_fwd.Greater(w).adjoint();
		gf_bwd.Retarded[w] = -gf_fwd.Lesser[w].adjoint() - gf_bwd.Lesser[w]; // This is to store ggtr (gret is actually zero for t < 0)
	}
	gf_bwd.Retarded[0] *= 0.5;
	return 0;
}


/** \brief <b> Calculates a difference mesaure between two Green's functions </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Computes the L2 norm (Euclidean distance)
 * > \f$ \|A - B\|_2 = \sqrt{ \sum_i |A_i - B_i|^2 } \f$
 * > between two `GF` objects in either time or frequency domain.
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param g1
 *  [GF] First Green's function
 * @param g2
 *  [GF] Second Green's function
 */
double GF2norm(GF &g1, GF &g2)
{
	double norm = 0.0;
	assert(g1.data_size_ == g2.data_size_);

	for (int w = 0; w < g1.data_size_ ; w++)
	{
		norm += std::norm(g1.g_ret[w] - g2.g_ret[w]);
		norm += std::norm(g1.g_les[w] - g2.g_les[w]);
	}
	norm = sqrt(norm) / (double) g1.data_size_;

	return norm;
}

/** \brief <b> Print out a Green's function into a text file </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > This prints out the frequenc or time `GF` Green's function into a text file of specified name.
 * > Output order of each line: Frequency (time), Real \f$[G^R]\f$, Imag \f$[G^R]\f$, Imag \f$[G^<]\f$, Real \f$[G^<]\f$...
 * > (Repeat in row-major order for other matrix entries.)
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param G
 *  [GF] Green's function
 * @param filename
 * [std::string] Filename
*/
void outputGF(GF &G, std::string filename){
	std::ofstream outputFile;
	outputFile.open(filename.c_str());
	outputFile << "# frequency/time" << "\t" << "retarded Real" << "\t"<< "retarded imag"<< "\t"<< "lesser imag"<< "\t"<< "lesser real ..."<< std::endl;
	long N=G.ngrid_;

	outputFile << std::setprecision(50);

	for (long w = 0; w < N; w++)
	{
		int i = w * G.el_size_;
		outputFile << G.grid_[w] << "\t";
		for (int j = 0; j < G.el_size_; j++)
			outputFile << G.g_ret[i + j].real() << "\t" << G.g_ret[i + j].imag() << "\t"<< G.g_les[i + j].imag() << "\t"<< G.g_les[i + j].real() << "\t";
		outputFile << std::endl;

	}
	
	outputFile.close();
}

///@private
/** \brief <b> Health check for a Green's function </b>*/
int healthCheck(GF &G){
	for(int idx = 0; idx < G.data_size_; idx++)
	{
		if(G.g_ret[idx] != G.g_ret[idx] || G.g_les[idx] != G.g_les[idx])
			return idx;
	}
	return -1;
}

/** \brief <b> Write data append style from a vector to a textfile. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Write data append style from a vector of doubles to a textfile with a given precision.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param data_vec
 * [std::vector<double>] Vector containing the data
 * @param fname
 * [std::string] Filename
 * @param prec
 * [int] Precision of the data to be stored in the text file
 */
void appendData(const std::vector<double> & data_vec, std::string fname, int prec)
{
	std::ofstream f(fname, std::ofstream::app);	
	f.precision(prec);
	for(const auto & data : data_vec)
	{
		f << data << "\t";
	}
	f << std::endl;
	f.close();
}


/** \brief <b> Read in a Green's function from a text file </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > This reads in the `GF` Green's function from a text file of specified name.
 * > Output order of each line: Frequency (time), Real \f$[G^R]\f$, Imag \f$[G^R]\f$, Imag \f$[G^<]\f$, Real \f$[G^<]\f$...
 * > (Repeat in row-major order for other matrix entries.)
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param G
 * [GF] Green's function
 * @param filename
 * [std::string] Filename
*/
void readGF(GF &G, std::string filename){
	std::ifstream outputFile;
	outputFile.open(filename.c_str());
	std::string line;
	std::getline(outputFile, line);
	long N=G.ngrid_;


	for (long w = 0; w < N; w++)
	{
		int i = w * G.el_size_;
		outputFile >> G.grid_[w];
		for (int j = 0; j < G.el_size_; j++)
		{
			double re, im;
			outputFile >> re >> im;
			G.g_ret[i + j] = re + NESS_II * im;
			outputFile >> im >> re;
			G.g_les[i + j] = re + NESS_II * im;
		}

	}

	outputFile.close();
}

/** \brief <b> Force the Green's function to be diagonal </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Force the `GF` Green's function to be diagonal in orbital space by discarding the off-diagonal elements,
 * > for all frequencies / times.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param gf
 *  [GF] Green's function
 */
void force_diag(GF &gf)
{
	for(int w = 0; w < gf.ngrid_; w++)
	{
		cdmatrix tmp(gf.size1_, gf.size2_);
		tmp = gf.Lesser[w].diagonal().asDiagonal();
		gf.Lesser[w] = tmp;
		tmp = gf.Retarded[w].diagonal().asDiagonal();
		gf.Retarded[w] = tmp;
	}

}

};//end of namespace
