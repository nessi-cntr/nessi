#ifndef CNTR_GETSET_IMPL_H
#define CNTR_GETSET_IMPL_H

#include "eigen_typedef.h"
#include "cntr_herm_matrix_decl.hpp"
#include "cntr_herm_matrix_timestep_decl.hpp"
#include "cntr_herm_matrix_timestep_view_decl.hpp"

namespace cntr {

//---------------------------------------------------------------------
//-------                  herm_matrix                          -------
//---------------------------------------------------------------------

// versions assuming hermitian symmetry
template <typename T>
void get_les(const int n, const int j, herm_matrix<T> &G, std::complex<T> &G_les){
	assert(n <= G.nt() && j <= G.nt());

	G.get_les(n, j, G_les);
}

template <typename T>  
void get_les(const int n, const int j, herm_matrix<T> &G, cdmatrix &G_les){
	assert(n <= G.nt() && j <= G.nt());

	G.get_les(n, j, G_les);
}

template <typename T>
void get_ret(const int n, const int j, herm_matrix<T> &G, std::complex<T> &G_ret){
	assert(n <= G.nt() && j <= G.nt());

	G.get_ret(n, j, G_ret);
}

template <typename T> template
void get_ret(const int n, const int j, herm_matrix<T> &G, cdmatrix &G_ret){
	assert(n <= G.nt() && j <= G.nt());

	G.get_ret(n, j, G_ret);
}

template <typename T>
void get_gtr(const int n, const int j, herm_matrix<T> &G, std::complex<T> &G_gtr){
	assert(n <= G.nt() && j <= G.nt());

	G.get_gtr(n, j, G_gtr);
}

template <typename T>
void get_gtr(const int n, const int j, herm_matrix<T> &G, cdmatrix &G_gtr){
	assert(n <= G.nt() && j <= G.nt());

	G.get_gtr(n, j, G_gtr);
}


template <typename T>
void get_tv(const int n, const int mtau, herm_matrix<T> &G, std::complex<T> &G_tv){
	assert(n <= G.nt() && mtau <= G.ntau());

	G.get_tv(n, mtau, G_tv);
}

template <typename T>
void get_tv(const int n, const int mtau, herm_matrix<T> &G, cdmatrix &G_tv){
	assert(n <= G.nt() && mtau <= G.ntau());

	G.get_tv(n, mtau, G_tv);
}

template <typename T>
void get_vt(const int mtau, const int n, herm_matrix<T> &G, std::complex<T> &G_vt){
	assert(n <= G.nt() && mtau <= G.ntau());

	G.get_vt(mtau, n, G_vt);
}

template <typename T>
void get_vt(const int mtau, const int n, herm_matrix<T> &G, cdmatrix &G_vt){
	assert(n <= G.nt() && j <= G.ntau());

	G.get_vt(mtau, n, G_tv);
}

template <typename T>
void get_mat(const int mtau, herm_matrix<T> &G, std::complex<T> &G_mat){
	assert(mtau <= G.ntau());

	G.get_mat(mtau, G_mat);
}

template <typename T>
void get_mat(const int mtau, herm_matrix<T> &G, cdmatrix &G_mat){
	assert(mtau <= G.ntau());

	G.get_mat(mtau, G_mat);
}


// versions without hermitian symmetry
template <typename T>
void get_les(const int n, const int j, herm_matrix<T> &G, herm_matrix<T> &Gcc, std::complex<T> &G_les){
	assert(n <= G.nt() && j <= G.nt());

	if(n <= j){
		G.get_les(n, j, G_les);
	} else {
		Gcc.get_les(j, n, G_les);
		G_les = -std::conj(G_les);
	}
}

template <typename T> template <class Matrix>
void get_les(const int n, const int j, herm_matrix<T> &G, herm_matrix<T> &Gcc, cdmatrix &G_les){
	assert(n <= G.nt() && j <= G.nt());

	if(n <= j){
		G.get_les(n, j, G_les);
	} else {
		Gcc.get_les(j, n, G_les);
		G_les = -G_les.adjoint();
	}
}

template <typename T>
void get_ret(const int n, const int j, herm_matrix<T> &G, herm_matrix<T> &Gcc, std::complex<T> &G_ret){
	assert(n <= G.nt() && j <= G.nt());

	if(j <= n){
		G.get_ret(n, j, G_ret);
	} else {
		Gcc.get_ret(j, n, G_ret);
		G_ret = -std::conj(G_ret);
	}
}

template <typename T>
void get_ret(const int n, const int j, herm_matrix<T> &G, herm_matrix<T> &Gcc, cdmatrix &G_ret){
	assert(n <= G.nt() && j <= G.nt());

	if(j <= n){
		G.get_ret(n, j, G_ret);
	} else {
		Gcc.get_ret(j, n, G_ret);
		G_ret = -G_ret.adjoint();
	}
}

template <typename T>
void get_gtr(const int n, const int j, herm_matrix<T> &G, herm_matrix<T> &Gcc, std::complex<T> &G_gtr){
	assert(n <= G.nt() && j <= G.nt());
	std::complex<T> G_ret, G_les;

	get_ret(n, j, G, Gcc, G_ret);
	get_les(n, j, G, Gcc, G_les);
	G_gtr = G_ret + G_les;
}

template <typename T>
void get_gtr(const int n, const int j, herm_matrix<T> &G, herm_matrix<T> &Gcc, cdmatrix &G_gtr){
	assert(n <= G.nt() && j <= G.nt());
	cdmatrix G_ret, G_les;

	get_ret(n, j, G, Gcc, G_ret);
	get_les(n, j, G, Gcc, G_les);
	G_gtr = G_ret + G_les;
}


template <typename T>
void get_tv(const int n, const int mtau, herm_matrix<T> &G, herm_matrix<T> &Gcc, std::complex<T> &G_tv){
	assert(n <= G.nt() && mtau <= G.ntau());

	G.get_tv(n, mtau, G_tv);
}

template <typename T> template
void get_tv(const int n, const int mtau, herm_matrix<T> &G, herm_matrix<T> &Gcc, cdmatrix &G_tv){
	assert(n <= G.nt() && mtau <= G.ntau());

	G.get_tv(n, mtau, G_tv);
}

template <typename T>
void get_vt(const int mtau, const int n, herm_matrix<T> &G, herm_matrix<T> &Gcc, std::complex<T> &G_vt){
	assert(n <= G.nt() && mtau <= G.ntau());
	std::complex<T> Gcc_tv;
	int ntau = G.ntau();
	int sig = G.sig();

	Gcc.get_tv(n, ntau - mtau, Gcc_tv);
	G_vt = -sig * std::conj(Gcc_tv);
}

template <typename T>
void get_vt(const int mtau, const int n, herm_matrix<T> &G, herm_matrix<T> &Gcc, Matrix &G_vt){
	assert(n <= G.nt() && mtau <= G.ntau());
	cdmatrix Gcc_tv;
	int ntau = G.ntau();
	int sig = G.sig();

	Gcc.get_tv(n, ntau - mtau, Gcc_tv);
	G_vt = -sig * Gcc_tv.adjoint();
}

//---------------------------------------------------------------------
//-------              herm_matrix_timestep                     -------
//---------------------------------------------------------------------

// versions assuming hermitian symmetry
template <typename T>
void get_les(const int n, const int j, herm_matrix_timestep<T> &G, std::complex<T> &G_les){
	assert(j == G.tstp_ && n <= G.tstp_);

	G.get_les(n, j, G_les);
}

template <typename T>  
void get_les(const int n, const int j, herm_matrix_timestep<T> &G, cdmatrix &G_les){
	assert(n <= G.nt() && j <= G.nt());

	G.get_les(n, j, G_les);
}

template <typename T>
void get_ret(const int n, const int j, herm_matrix_timestep<T> &G, std::complex<T> &G_ret){
	assert(n <= G.nt() && j <= G.nt());

	G.get_ret(n, j, G_ret);
}

template <typename T> template
void get_ret(const int n, const int j, herm_matrix_timestep<T> &G, cdmatrix &G_ret){
	assert(n <= G.nt() && j <= G.nt());



	G.get_ret(n, j, G_ret);
}

template <typename T>
void get_gtr(const int n, const int j, herm_matrix<T> &G, std::complex<T> &G_gtr){
	assert(n <= G.nt() && j <= G.nt());

	G.get_gtr(n, j, G_gtr);
}

template <typename T>
void get_gtr(const int n, const int j, herm_matrix<T> &G, cdmatrix &G_gtr){
	assert(n <= G.nt() && j <= G.nt());

	G.get_gtr(n, j, G_gtr);
}


template <typename T>
void get_tv(const int n, const int mtau, herm_matrix<T> &G, std::complex<T> &G_tv){
	assert(n <= G.nt() && mtau <= G.ntau());

	G.get_tv(n, mtau, G_tv);
}

template <typename T>
void get_tv(const int n, const int mtau, herm_matrix<T> &G, cdmatrix &G_tv){
	assert(n <= G.nt() && mtau <= G.ntau());

	G.get_tv(n, mtau, G_tv);
}

template <typename T>
void get_vt(const int mtau, const int n, herm_matrix<T> &G, std::complex<T> &G_vt){
	assert(n <= G.nt() && mtau <= G.ntau());

	G.get_vt(mtau, n, G_vt);
}

template <typename T>
void get_vt(const int mtau, const int n, herm_matrix<T> &G, cdmatrix &G_vt){
	assert(n <= G.nt() && j <= G.ntau());

	G.get_vt(mtau, n, G_tv);
}

template <typename T>
void get_mat(const int mtau, herm_matrix<T> &G, std::complex<T> &G_mat){
	assert(mtau <= G.ntau());

	G.get_mat(mtau, G_mat);
}

template <typename T>
void get_mat(const int mtau, herm_matrix<T> &G, cdmatrix &G_mat){
	assert(mtau <= G.ntau());

	G.get_mat(mtau, G_mat);
}


} // namespace cntr


#endif  // CNTR_GETSET_IMPL_H