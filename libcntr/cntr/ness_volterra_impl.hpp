#include "ness.hpp"
#include <vector>
#include <array>

namespace ness
{

/** \brief  Solves the Dyson equation in frequency space 
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Solves the Dyson equation in frequency space:
 * > \f$ G(\omega) = 1/[ \omega -H  -\Sigma(\omega) - i\eta]\f$
 * > for the retarded and lesser component of the `GF` Green's function for all frequencies.
 * > Here, are given: \f$\Sigma(\omega)\f$, \f$H\f$, and \f$\eta\f$.
 *
 * \note: Check if algorithm works with setting \f$\eta=0\f$, only needed if resulting Green's functions are slowly decaying and in combination with FFT.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param G
 *  [GF] Green's function to store the solution
 * @param Sigma
 *  [GF] Self-energy
 * @param H
 *  [cdmatrix] hamiltonian and chemical potential
 * @param eta
 *  [double] regulator for stability
 */
void dyson_ness(GF & G, GF & Sigma, const cdmatrix & H, double eta)
{
        assert(G.ngrid_ == Sigma.ngrid_);
        assert(G.gf_type_ == freq_gf);
        assert(Sigma.gf_type_ == freq_gf);
        assert(G.size1_ == Sigma.size1_);
        int size = Sigma.size1_;
        cdmatrix Id = cdmatrix::Identity(size, size);
                for(auto w = 0; w < G.ngrid_; ++w)
                {
                        G.Retarded[w] = (Id*G.grid_[w] - H - NESS_II * Id * eta - Sigma.Retarded[w]).inverse();
                        G.Lesser[w] = G.Retarded[w] * Sigma.Lesser[w] * G.Retarded[w].adjoint();

                }
}

/** \brief  Solves the Dyson equation in frequency space using 'GF_pair' objects 
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Solves the Dyson equation in frequency space:
 * > \f$ G(\omega) = 1/[ \omega -H  -\Sigma(\omega) - i\eta]\f$
 * > for the retarded and lesser component of the `GF_pair` Green's function for all frequencies.
 * > Here, are given: \f$\Sigma(\omega)\f$, \f$H\f$, and \f$ \eta\f$.
 * > This version of the function uses the time members of 'GF_pair' objects.
 *
 * \note: Check if algorithm works with setting \f$\eta=0\f$, only needed if resulting Green's functions are slowly decaying and in combination with FFT.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param  G
 *  [GF_pair] Green's function pair to store the solution (time)
 * @param  Sigma
 *  [GF_pair] Self-energy, Green's function pair (time)
 * @param H
 *  [cdmatrix] hamiltonian and chemical potential
 * @param eta
 *  [double] regulator for stability
 */
void dyson_ness(GF_pair & G, GF_pair & Sigma, const cdmatrix & H, double eta)
{
        assert(G.time.ngrid_ == Sigma.time.ngrid_);
        assert(G.time.size1_ == Sigma.time.size1_);
        assert(G.freq.ngrid_ == Sigma.freq.ngrid_);
        assert(G.freq.size1_ == Sigma.freq.size1_);
        int size = Sigma.freq.size1_;
        cdmatrix Id = cdmatrix::Identity(size, size);
        Sigma.to_freq();
                for(auto w = 0; w < G.freq.ngrid_; ++w)
                {
                    G.freq.Retarded[w].noalias() = (Id*G.freq.grid_[w] - H - NESS_II * Id * eta - Sigma.freq.Retarded[w]).inverse();
                    G.freq.Lesser[w].noalias() = G.freq.Retarded[w] * Sigma.freq.Lesser[w] * G.freq.Retarded[w].adjoint();
                }
        G.to_time();
}


/** \brief  Solve the volterra equation \f$(I + F) \ast G = Q\f$ 
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Solve the volterra equation \f$(I + F) \ast G = Q\f$ for \f$G\f$,
 * > with \f$ F \f$, \f$ Q \f$ hermitian Green's functions and \f$ I \f$ the identity.
 * > Steady state version of real time ´vie2´ function.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param  G
 *  [GF] Green's function to be solved for
 * @param  F
 *  [GF] Green's function
 * @param  Q
 *  [GF] Green's function
*/
void vie2_ness(GF & G, const GF & F, const GF & Q)
{
    assert(F.ngrid_ == Q.ngrid_);
    assert(G.ngrid_ == Q.ngrid_);
    assert(G.gf_type_ == freq_gf);
    assert(F.gf_type_ == freq_gf);
    assert(Q.gf_type_ == freq_gf);
    assert(F.size1_ == Q.size1_);
    assert(Q.size1_ == G.size1_);
    int size = F.size1_;
    cdmatrix Id = cdmatrix::Identity(size, size);
    for(auto w = 0; w < G.ngrid_; ++w)
    {
            G.Retarded[w].noalias() = (Id + F.Retarded[w]).inverse() *
                Q.Retarded[w];
            G.Lesser[w].noalias() = (Id + F.Retarded[w]).inverse() *
                (Q.Lesser[w] - F.Lesser[w] * G.Retarded[w].adjoint());
    }
}

/** \brief  Solve the volterra equation \f$(I + F) \ast G = Q\f$ 
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Solve the volterra equation \f$(I + F) \ast G = Q\f$ for \f$G\f$,
 * > with \f$ F \f$, \f$ Q \f$ hermitian Green's functions and \f$ I \f$ the identity.
 * > Steady state version of real time ´vie2´ function using Green's function pairs (time).
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param  G
 * [GF_pair] Green's function pair to store solution in
 * @param  F
 * [GF_pair] Green's function pair
 * @param  Q
 * [GF_pair] Green's function pair
 */
void vie2_ness(GF_pair & G, GF_pair & F, GF_pair & Q)
{
    assert(F.time.ngrid_ == Q.time.ngrid_);
    assert(G.time.ngrid_ == Q.time.ngrid_);
    assert(F.time.size1_ == Q.time.size1_);
    assert(Q.time.size1_ == G.time.size1_);
    int size = F.time.size1_;
    cdmatrix Id = cdmatrix::Identity(size, size);
    F.to_freq();
    Q.to_freq();
    for(auto w = 0; w < G.freq.ngrid_; ++w)
    {
            G.freq.Retarded[w].noalias() = (Id + F.freq.Retarded[w]).inverse() *
                Q.freq.Retarded[w];
            G.freq.Lesser[w].noalias() = (Id + F.freq.Retarded[w]).inverse() *
                (Q.freq.Lesser[w] - F.freq.Lesser[w] * G.freq.Retarded[w].adjoint());
    }
    G.to_time();
}


/** \brief  Solves the Dyson equation in frequency space given the free propagator 
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Solves the Dyson equation in frequency space:
 * > \f$ G(\omega) = 1/[ G_0(\omega)^{-1}  -\Sigma(\omega)]\f$
 * > for the retarded and lesser component of the Green's function for all frequencies.
 * > Here, are given: \f$\Sigma(\omega)\f$, \f$G_0(\omega)^{-1}\f$.
 *
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param invg0
 *  [GF] Inverse free propagator
 * @param self
 *  [GF] Self-energy
 * @param g
 *  [GF] Green's function to store the solution
 */
int dyson(GF &invg0, GF &self, GF &g)
{
    assert(self.data_size_ == invg0.data_size_);
    cplx temp[self.el_size_];
    cplx tempc[self.el_size_];
    //gret
    for (int w = 0; w < self.ngrid_; w++)
    {
        g.Retarded[w] = invg0.Retarded[w] - self.Retarded[w];
    }
    //glss
    for (int w = 0; w < g.ngrid_; w++)
    {
        g.Lesser[w] = g.Retarded[w] * self.Lesser[w] * g.Retarded[w].adjoint();
    }
    return 0;
}


/** \brief  Solves the Dyson equation self-consistently 
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Solves the Dyson equation in frequency space:
 * > \f$ G(\omega) = 1/[ G_{old}(\omega)^{-1}  -\Sigma(\omega)]\f$
 * > for the retarded and lesser component of the Green's function for all frequencies.
 * > Here, are given: \f$\Sigma(\omega)\f$, \f$G_{old}(\omega)^{-1}\f$.
 *
 * \note: Here, the interacting propagator from the previous step is getting dressed with \f$ \Sigma \f$.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param g
 *  [GF] Green's function to store the solution and providing the previous step
 * @param self
 *  [GF] Self-energy
 */
int dyson_from_inv(GF &g, GF &self)
{
    assert(self.data_size_ == g.data_size_);
    assert(self.size1_ == g.size1_);
    assert(self.size2_ == g.size2_);

    cplx temp[self.el_size_];
    cplx tempc[self.el_size_];
    for (int w = 0; w < self.ngrid_; w++)
    {
        //gret
        g.Retarded[w] = (g.Retarded[w] - self.Retarded[w]).inverse();
        g.Lesser[w] = g.Retarded[w] * self.Lesser[w] * g.Retarded[w].adjoint();

    }
    return 0;
}

/** \brief  Solve the volterra equation \f$(I + B) \ast {\tt gf} = C\f$ 
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Solve the volterra equation \f$(I + B) \ast {\tt gf} = C\f$ for `gf` in frequency domain,
 * > with \f$B,C\f$ Green's functions and \f$ I\f$ the identity.
 *
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param  gf
 *  [GF] Green's function to be solved for
 * @param  B
 *  [GF] Green's function
 * @param  C
 *  [GF] Green's function
 * @param  Id
 *  [cdmatrix] Identity matrix of specified size
 */
void volterra_IBC(GF & gf, const GF & B, const GF & C, const cdmatrix & Id)
{
	assert(B.ngrid_ == C.ngrid_);
	assert(gf.ngrid_ == C.ngrid_);
	assert(gf.gf_type_ == freq_gf);
	assert(B.gf_type_ == freq_gf);
	assert(C.gf_type_ == freq_gf);
	assert(B.size1_ == C.size1_);
	assert(C.size1_ == gf.size1_);

	for(auto w = 0; w < gf.ngrid_; ++w)
	{
			gf.Retarded[w].noalias() = (Id + B.Retarded[w]).inverse() * 
				C.Retarded[w];
			gf.Lesser[w].noalias() = (Id + B.Retarded[w]).inverse() * 
				(C.Lesser[w] - B.Lesser[w] * gf.Retarded[w].adjoint());
	}
}

/** \brief  Solve the volterra equation \f$(I + AB) \ast {\tt gf} = C\f$ where \f$ A\f$ is a constant 
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Solve the volterra equation \f$(I + AB) \ast {\tt gf} = C\f$  for `gf` in frequency domain, where \f$ A\f$ is a constant
 * > with \f$ B,C\f$ Green's functions, \f$ A\f$ a matrix and \f$ I\f$ the identity.
 *
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param  gf
 *  [GF] Green's function to be solved for
 * @param  A
 *  [cdmatrix] Constant matrix
 * @param  B
 *  [GF] Green's function
 * @param  C
 *  [GF] Green's function
 * @param  Id
 *  [cdmatrix] Identity matrix of specified size
 */
void volterra_IABC(GF & gf, const cdmatrix & A, const GF & B, const GF & C, const cdmatrix & Id)
{
	assert(B.ngrid_ == C.ngrid_);
	assert(gf.ngrid_ == C.ngrid_);
	assert(gf.gf_type_ == freq_gf);
	assert(B.gf_type_ == freq_gf);
	assert(C.gf_type_ == freq_gf);
	assert(A.rows() == C.size1_);
	assert(B.size1_ == C.size1_);
	assert(C.size1_ == gf.size1_);

	for(auto w = 0; w < gf.ngrid_; ++w)
	{
			gf.Retarded[w].noalias() = (Id + A * B.Retarded[w]).inverse() * 
				C.Retarded[w];
			gf.Lesser[w].noalias() = (Id + A * B.Retarded[w]).inverse() * 
				(C.Lesser[w] - A * B.Lesser[w] * gf.Retarded[w].adjoint());
	}
}

/** \brief  Solve the volterra equation \f$(I + AB) \ast {\tt gf} = C\f$ where \f$ B\f$ is a constant 
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Solve the volterra equation \f$(I + AB) \ast {\tt gf} = C\f$  for `gf` in frequency domain, where \f$ B\f$  is a constant
 * > with \f$ A,C\f$ Green's functions, \f$ B\f$ a matrix and \f$ I\f$ the identity.
 *
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param  gf
 *  [GF] Green's function to be solved for
 * @param  A
 *  [GF] Green's function
 * @param  b
 *  [cdmatrix] Constant matrix
 * @param  B
 *  [GF] Green's function
 * @param  C
 *  [GF] Green's function
 * @param  Id
 *  [cdmatrix] Identity matrix of specified size
 */
void volterra_IABC(GF & gf, const GF & A, const cdmatrix & B, const GF & C, const cdmatrix & Id)
{
	assert(A.ngrid_ == C.ngrid_);
	assert(gf.ngrid_ == C.ngrid_);
	assert(gf.gf_type_ == freq_gf);
	assert(A.gf_type_ == freq_gf);
	assert(C.gf_type_ == freq_gf);
	assert(A.size1_ == C.size1_);
	assert(B.rows() == C.size1_);
	assert(C.size1_ == gf.size1_);

	for(auto w = 0; w < gf.ngrid_; ++w)
	{
			gf.Retarded[w].noalias() = (Id + A.Retarded[w] * B).inverse() * 
				C.Retarded[w];
			gf.Lesser[w].noalias() = (Id + A.Retarded[w] * B).inverse() * 
				(C.Lesser[w] - A.Lesser[w] * B * gf.Retarded[w].adjoint());
	}
}

/** \brief  Solve the volterra equation \f$(I + AB) \ast {\tt gf} = C\f$ where \f$ A, C\f$ are constants
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Solve the volterra equation \f$(I + AB) \ast {\tt gf} = C\f$  for `gf` in frequency domain, where \f$ A, C\f$ are constants
 * > with  \f$ B\f$ a Green's function,  \f$ A, C\f$ matrices and  \f$ I\f$ the identity.
 *
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param  gf
 *  [GF] Green's function to be solved for
 * @param  A
 *  [cdmatrix] Constant matrix
 * @param  B
 *  [GF] Green's function
 * @param  C
 *  [cdmatrix] Constant matrix
 * @param  Id
 *  [cdmatrix] Identity matrix of specified size
 */
void volterra_IABC(GF & gf, const cdmatrix & A, const GF & B, const cdmatrix & C, const cdmatrix & Id)
{
	assert(gf.ngrid_ == B.ngrid_);
	assert(gf.gf_type_ == freq_gf);
	assert(B.gf_type_ == freq_gf);
	assert(A.rows() == C.rows());
	assert(B.size1_ == C.rows());
	assert(B.size1_ == gf.size1_);

	for(auto w = 0; w < gf.ngrid_; ++w)
	{
			gf.Retarded[w].noalias() = (Id + A * B.Retarded[w]).inverse() * 
				C;
			gf.Lesser[w].noalias() = (Id + A * B.Retarded[w]).inverse() * 
				(- A * B.Lesser[w] * gf.Retarded[w].adjoint());
	}
}


/** \brief  Solve the volterra equation \f$(I + AB) \ast  = C\f$ 
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Solve the volterra equation \f$(I + AB) \ast {\tt gf} = C\f$  for `gf` in frequency domain,
 * > with \f$  A,B,C\f$ Green's functions and \f$ I \f$ the identity.
 *
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param  gf
 *  [GF] Green's function to be solved for
 * @param  A
 *  [GF] Green's function
 * @param  B
 *  [GF] Green's function
 * @param  C
 *  [GF] Green's function
 * @param  Id
 *  [cdmatrix] Identity matrix of specified size
 */
void volterra_IABC(GF & gf, const GF & A, const GF & B, const GF & C, const cdmatrix & Id)
{
	assert(A.ngrid_ == B.ngrid_);
	assert(B.ngrid_ == C.ngrid_);
	assert(gf.ngrid_ == C.ngrid_);
	assert(gf.gf_type_ == freq_gf);
	assert(A.gf_type_ == freq_gf);
	assert(B.gf_type_ == freq_gf);
	assert(C.gf_type_ == freq_gf);
	assert(A.size1_ == B.size1_);
	assert(B.size1_ == C.size1_);
	assert(C.size1_ == gf.size1_);

	for(auto w = 0; w < gf.ngrid_; ++w)
	{
			gf.Retarded[w].noalias() = (Id + A.Retarded[w] * B.Retarded[w]).inverse() * 
				C.Retarded[w];
			gf.Lesser[w].noalias() = (Id + A.Retarded[w] * B.Retarded[w]).inverse() * 
				(C.Lesser[w] - (A.Retarded[w] * B.Lesser[w]
				 + A.Lesser[w] * B.Retarded[w].adjoint()) *  gf.Retarded[w].adjoint());
	}

}

/** \brief  Solve the volterra equation \f$(I + AB) \ast {\tt gf} = UCU\f$ where \f$ A\f$ and \f$ U\f$ are constants 
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Solve the volterra equation \f$ (I + AB) \ast {\tt gf} = UCU \f$  for `gf` in frequency domain, where \f$ A \f$ and \f$ U\f$ are constants,
 * > with \f$ B,C\f$ Green's functions, \f$ A,U\f$ matrices and \f$ I\f$ the identity.
 *
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param  gf
 *  [GF] Green's function to be solved for
 * @param  A
 *  [cdmatrix] Green's function
 * @param  B
 *  [GF] Green's function
 * @param  C
 *  [GF] Green's function
 * @param  U
 *  [cdmatrix] Green's function
 * @param  Id
 *  [cdmatrix] Identity matrix of specified size
 */
void volterra_IABUCU(GF & gf, const cdmatrix & A, const GF & B, const GF & C, cdmatrix U, const cdmatrix & Id)
{
	assert(B.ngrid_ == C.ngrid_);
	assert(gf.ngrid_ == C.ngrid_);
	assert(gf.gf_type_ == freq_gf);
	assert(B.gf_type_ == freq_gf);
	assert(C.gf_type_ == freq_gf);
	assert(A.rows() == C.size1_);
	assert(U.rows() == C.size1_);
	assert(B.size1_ == C.size1_);
	assert(C.size1_ == gf.size1_);

	for(auto w = 0; w < gf.ngrid_; ++w)
	{
			gf.Retarded[w].noalias() = (Id + A * B.Retarded[w]).inverse() * 
				U * C.Retarded[w] * U;
			gf.Lesser[w].noalias() = (Id + A * B.Retarded[w]).inverse() * 
				(U * C.Lesser[w] * U - A * B.Lesser[w] * gf.Retarded[w].adjoint());
	}
}



/** \brief  Solve the volterra equation \f$(I + A(L+M)) \ast {\tt gf} = C\f$  
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Solve the volterra equation  \f$(I + A(L+M)) \ast {\tt gf} = C\f$ for `gf` in frequency domain,
 * > with \f$ A,L,M,C\f$  Green's functions and \f$ I\f$ the identity.
 *
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param  gf
 *  [GF] Green's function to be solved for
 * @param  A
 *  [GF] Green's function
 * @param  L
 *  [GF] Green's function
 * @param  M
 *  [GF] Green's function
 * @param  C
 *  [GF] Green's function
 * @param  Id
 *  [cdmatrix] Identity matrix of specified size
 */
void volterra_IALMC(GF & gf, const GF & A, const GF & L, const GF & M, const GF & C, const cdmatrix & Id)
{
	assert(A.ngrid_ == L.ngrid_);
	assert(L.ngrid_ == M.ngrid_);
	assert(M.ngrid_ == C.ngrid_);
	assert(gf.ngrid_ == C.ngrid_);
	assert(gf.gf_type_ == freq_gf);
	assert(A.gf_type_ == freq_gf);
	assert(L.gf_type_ == freq_gf);
	assert(M.gf_type_ == freq_gf);
	assert(C.gf_type_ == freq_gf);
	assert(A.size1_ == L.size1_);
	assert(L.size1_ == M.size1_);
	assert(L.size1_ == C.size1_);
	assert(C.size1_ == gf.size1_);

	for(auto w = 0; w < gf.ngrid_; ++w)
	{
			gf.Retarded[w].noalias() = (Id + A.Retarded[w] * (L.Retarded[w] + M.Retarded[w])).inverse() * 
				C.Retarded[w];
			gf.Lesser[w].noalias() = (Id + A.Retarded[w] * (L.Retarded[w] + M.Retarded[w])).inverse() * 
				(C.Lesser[w] - (A.Retarded[w] * (L.Lesser[w] + M.Lesser[w])
				 + A.Lesser[w] * (L.Retarded[w] + M.Retarded[w]).adjoint()) *  gf.Retarded[w].adjoint());
	}

}

/** \brief  Alternative way of solving the volterra equation \f$(I + B) \ast {\tt gf} = C\f$ 
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Solve the volterra equation \f$(I + B) \ast {\tt gf} = C\f$ for `gf`in frequency domain,
 * > with \f$ B,C\f$ Green's functions and \f$ I\f$ the identity.
 * > Calculate the retarded component by first calculating the greater component.
 *
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param  gf
 *  [GF] Green's function to be solved for
 * @param  B
 *  [GF] Green's function
 * @param  C
 *  [GF] Green's function
 * @param  Id
 *  [cdmatrix] Identity matrix of specified size
 */
void volterra_IBCgtr(GF & gf, const GF & B, const GF & C, const cdmatrix & Id)
{
	assert(B.ngrid_ == C.ngrid_);
	assert(gf.ngrid_ == C.ngrid_);
	assert(gf.gf_type_ == freq_gf);
	assert(B.gf_type_ == freq_gf);
	assert(C.gf_type_ == freq_gf);
	assert(B.size1_ == C.size1_);
	assert(C.size1_ == gf.size1_);

	for(int w = 0; w < gf.ngrid_; ++w) 
	{			
		gf.Lesser[w].noalias() = (Id + B.Retarded[w]).inverse() * 
				(C.Lesser[w] - B.Lesser[w] * gf.Retarded[w].adjoint());
		gf.Retarded[w].noalias() = (Id + B.Retarded[w]).inverse() * 
				(C.Greater(w) - B.Greater(w) * gf.Retarded[w].adjoint());

		gf.Retarded[w] -= gf.Lesser[w];
		gf.Retarded[w] = (gf.Retarded[w] - gf.Retarded[w].adjoint()) * 0.5;
	}
}


};//namespace


