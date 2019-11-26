/** \brief <b> Returns the Matsubara component of a general contour function at given times.  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Returns the lesser component \f$ G^\mathrm{M}(\tau_m) \f$ at given imaginary time \f$ \tau_m\f$.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param m
* > [int] Index of time \f$ \tau_m\f$ .
* @param G_mat
* > [complex<T>] The Matsubara component (scalar GF).
* @param G
* > [herm_matrix] Contour function G
* @param Gcc
* > [herm_matrix] Hermitian conjugate \f$G^\ddagger\f$ of \f$G\f$.
*/
template <typename T>
inline void get_mat(const int m, std::complex<T> &G_mat, herm_matrix<T> &G, herm_matrix<T> &Gcc){
    assert(m <= G.ntau());
    assert(G.ntau() == Gcc.ntau());
    assert(G.size1() == Gcc.size1());
    assert(G.size2() == Gcc.size2());
    assert(G.sig() == Gcc.sig());

    G_mat = *G.matptr(m);
}

/** \brief <b> Returns the Matsubara component of a general contour function at given times.  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Returns the lesser component \f$ G^\mathrm{M}(\tau_m) \f$ at given imaginary time \f$ \tau_m\f$.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param m
* > [int] Index of time \f$ \tau_m\f$ .
* @param G_mat
* > [complex<T>] The Matsubara component (scalar GF).
* @param G
* > [herm_matrix] Contour function G
* @param Gcc
* > [herm_matrix] Hermitian conjugate \f$G^\ddagger\f$ of \f$G\f$.
*/
template <typename T>
void get_mat(const int m, cdmatrix &G_mat, herm_matrix<T> &G, herm_matrix<T> &Gcc){
	assert(m <= G.ntau());
    assert(G.ntau() == Gcc.ntau());
    assert(G.size1() == Gcc.size1());
    assert(G.size2() == Gcc.size2());
    assert(G.sig() == Gcc.sig());
	std::complex<T> *mat;
	int size1=G.size1(), size2=G.size2();

    mat = G.matptr(m);
    map_ptr2matrix(size1, size2, mat, G_mat);
}



/** \brief <b> Returns the Matsubara component of a general contour function at given times.  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Returns the lesser component \f$ G^\mathrm{M}(\tau_m) \f$ at given imaginary time \f$ \tau_m\f$.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param m
* > [int] Index of time \f$ \tau_m\f$ .
* @param G_mat
* > [complex<T>] The Matsubara component (scalar GF).
* @param G
* > [herm_matrix] Contour function G
*/
template <typename T>
inline void get_mat(const int m, std::complex<T> &G_mat, herm_matrix<T> &G){
    assert(m <= G.ntau());

    G_mat = *G.matptr(m);
}

/** \brief <b> Returns the Matsubara component of a general contour function at given times.  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Returns the lesser component \f$ G^\mathrm{M}(\tau_m) \f$ at given imaginary time \f$ \tau_m\f$.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param m
* > [int] Index of time \f$ \tau_m\f$ .
* @param G_mat
* > [complex<T>] The Matsubara component (scalar GF).
* @param G
* > [herm_matrix] Contour function G
*/
template <typename T>
void get_mat(const int m, cdmatrix &G_mat, herm_matrix<T> &G){
    assert(m <= G.ntau());
    int size1=G.size1(), size2=G.size2();
    std::complex<T> *mat;

    mat = G.matptr(m);
    map_ptr2matrix(size1, size2, mat, G_mat);
}



/** \brief <b> Returns the lesser component of a general contour function at given times.  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Returns the lesser component \f$ G^<(t_i,t_j) \f$ at given times \f$ t_i\f$
* > and \f$ t_j\f$. For \f$ i > j\f$ the hermitian conjugate \f$G^\ddagger\f$ is used.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > [int] Index of time \f$ t_i\f$ .
* @param j
* > [int] Index of time \f$ t_j\f$ .
* @param G_les
* > [complex<T>] The lesser component (scalar).
* @param G
* > [herm_matrix] Contour function G
* @param Gcc
* > [herm_matrix] Hermitian conjugate \f$G^\ddagger\f$ of \f$G\f$.
*/
template <typename T>
inline void get_les(const int i, const int j, std::complex<T> &G_les, herm_matrix<T> &G, herm_matrix<T> &Gcc){
    assert(i <= G.nt() && j <= G.nt());
    assert(i <= Gcc.nt() && j <= Gcc.nt());
    assert(G.nt() == Gcc.nt());
    assert(G.ntau() == Gcc.ntau());
    assert(G.size1() == Gcc.size1());
    assert(G.size2() == Gcc.size2());
    assert(G.sig() == Gcc.sig());
    int size1 = G.size1(), size2 = G.size2();

    if (i <= j) {
        G_les = *G.lesptr(i,j);
    } else {
        G_les = *Gcc.lesptr(j,i);
        G_les = -std::conj(G_les);
    }
}



/** \brief <b> Returns the lesser component of a general contour function at given times.  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Returns the lesser component \f$ G^<(t_i,t_j) \f$ at given times \f$ t_i\f$
* > and \f$ t_j\f$. For \f$ i > j\f$ the hermitian conjugate \f$G^\ddagger\f$ is used.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > [int] Index of time \f$ t_i\f$ .
* @param j
* > [int] Index of time \f$ t_j\f$ .
* @param G_les
* > [cdmatrix] The lesser component (returned to an eigen3 matrix).
* @param G
* > [herm_matrix] Contour function G
* @param Gcc
* > [herm_matrix] Hermitian conjugate \f$G^\ddagger\f$ of \f$G\f$.
*/
template <typename T>
void get_les(const int i, const int j, cdmatrix &G_les, herm_matrix<T> &G, herm_matrix<T> &Gcc){
	assert(i <= G.nt() && j <= G.nt());
	assert(i <= Gcc.nt() && j <= Gcc.nt());
	assert(G.nt() == Gcc.nt());
	assert(G.ntau() == Gcc.ntau());
	assert(G.size1() == Gcc.size1());
	assert(G.size2() == Gcc.size2());
	assert(G.sig() == Gcc.sig());
	int size1 = G.size1(), size2 = G.size2();
	std::complex<T> *les;

	if (i <= j) {
		les = G.lesptr(i,j);
		map_ptr2matrix(size1, size2, les, G_les);
	} else {
		les = Gcc.lesptr(j,i);
		map_ptr2matrix(size1, size2, les, G_les);
		G_les.adjointInPlace();
		G_les = -G_les;
	}
}


/** \brief <b> Returns the retarded component of a general contour function at given times.  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Returns the retarded component \f$ G^\mathrm{R}(t_i,t_j) \f$ at given times \f$ t_i\f$
* > and \f$ t_j\f$. We assume \f$G^\mathrm{R}(t_i,t_j)\f$ can be analytically continued
* > to \f$ j > i\f$, for which the hermitian conjugate \f$G^\ddagger\f$ is used.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > [int] Index of time \f$ t_i\f$ .
* @param j
* > [int] Index of time \f$ t_j\f$ .
* @param G_ret
* > [complex<T>] The retarded component (scalar GF).
* @param G
* > [herm_matrix] Contour function G
* @param Gcc
* > [herm_matrix] Hermitian conjugate \f$G^\ddagger\f$ of \f$G\f$.
*/
template <typename T>
inline void get_ret(const int i, const int j, std::complex<T> &G_ret, herm_matrix<T> &G, herm_matrix<T> &Gcc){
	assert(i <= G.nt() && j <= G.nt());
	assert(i <= Gcc.nt() && j <= Gcc.nt());
	assert(G.nt() == Gcc.nt());
	assert(G.ntau() == Gcc.ntau());
	assert(G.size1() == Gcc.size1());
	assert(G.size2() == Gcc.size2());
	assert(G.sig() == Gcc.sig());

	if (j <= i) {
		G_ret = *G.retptr(i,j);
	} else {
		G_ret = *Gcc.retptr(j,i);
		G_ret = -std::conj(G_ret);
	}
}

/** \brief <b> Returns the retarded component of a general contour function at given times.  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Returns the retarded component \f$ G^\mathrm{R}(t_i,t_j) \f$ at given times \f$ t_i\f$
* > and \f$ t_j\f$. We assume \f$G^\mathrm{R}(t_i,t_j)\f$ can be analytically continued
* > to \f$ j > i\f$, for which the hermitian conjugate \f$G^\ddagger\f$ is used.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > [int] Index of time \f$ t_i\f$ .
* @param j
* > [int] Index of time \f$ t_j\f$ .
* @param G_ret
* > [cdmatrix] The retarded component (returned to an eigen3 matrix).
* @param G
* > [herm_matrix] Contour function G
* @param Gcc
* > [herm_matrix] Hermitian conjugate \f$G^\ddagger\f$ of \f$G\f$.
*/
template <typename T>
void get_ret(const int i, const int j, cdmatrix &G_ret, herm_matrix<T> &G, herm_matrix<T> &Gcc){
	assert(i <= G.nt() && j <= G.nt());
	assert(i <= Gcc.nt() && j <= Gcc.nt());
	assert(G.nt() == Gcc.nt());
	assert(G.ntau() == Gcc.ntau());
	assert(G.size1() == Gcc.size1());
	assert(G.size2() == Gcc.size2());
	assert(G.sig() == Gcc.sig());
	int size1 = G.size1(), size2 = G.size2();
	std::complex<T> *ret;

	if (i >= j) {
		ret = G.retptr(i,j);
		map_ptr2matrix(size1, size2, ret, G_ret);
	} else {
		ret = Gcc.retptr(j,i);
		map_ptr2matrix(size1, size2, ret, G_ret);
		G_ret.adjointInPlace();
		G_ret = -G_ret;
	}
}



/** \brief <b> Returns the left-mixing component of a general contour function at given times.  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Returns the left-mixing component \f$ G^\rceil(t_i,\tau_m) \f$ at given times \f$ t_i\f$
* > and \f$ \tau_m\f$.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > [int] Index of time \f$ t_i\f$ .
* @param m
* > [int] Index of time \f$ \tau_m\f$ .
* @param G_tv
* > [complex<T>] The left-mixing component (scalar GF).
* @param G
* > [herm_matrix] Contour function G
* @param Gcc
* > [herm_matrix] Hermitian conjugate \f$G^\ddagger\f$ of \f$G\f$.
*/
template <typename T>
inline void get_tv(const int i, const int m, std::complex<T> &G_tv, herm_matrix<T> &G, herm_matrix<T> &Gcc){
	assert(i <= G.nt() && m <= G.ntau());
	assert(i <= Gcc.nt() && m <= Gcc.ntau());
	assert(G.nt() == Gcc.nt());
	assert(G.ntau() == Gcc.ntau());
	assert(G.size1() == Gcc.size1());
	assert(G.size2() == Gcc.size2());
	assert(G.sig() == Gcc.sig());

	G_tv = *G.tvptr(i,m);
}

/** \brief <b> Returns the left-mixing component of a general contour function at given times.  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Returns the left-mixing component \f$ G^\rceil(t_i,\tau_m) \f$ at given times \f$ t_i\f$
* > and \f$ \tau_m\f$.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > [int] Index of time \f$ t_i\f$ .
* @param m
* > [int] Index of time \f$ \tau_m\f$ .
* @param G_tv
* > [cdmatrix] The left-mixing component (returned to an eigen3 matrix).
* @param G
* > [herm_matrix] Contour function G
* @param Gcc
* > [herm_matrix] Hermitian conjugate \f$G^\ddagger\f$ of \f$G\f$.
*/
template <typename T>
void get_tv(const int i, const int m, cdmatrix &G_tv, herm_matrix<T> &G, herm_matrix<T> &Gcc){
	assert(i <= G.nt() && m <= G.ntau());
	assert(i <= Gcc.nt() && m <= Gcc.ntau());
	assert(G.nt() == Gcc.nt());
	assert(G.ntau() == Gcc.ntau());
	assert(G.size1() == Gcc.size1());
	assert(G.size2() == Gcc.size2());
	assert(G.sig() == Gcc.sig());
	int size1 = G.size1(), size2 = G.size2();
	std::complex<T> *tv;

	tv = G.tvptr(i,m);
	map_ptr2matrix(size1, size2, tv, G_tv);
}



/** \brief <b> Returns the right-mixing component of a general contour function at given times.  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Returns the right-mixing component \f$ G^\lceil(\tau_m,t_i) \f$ at given times \f$ t_i\f$
* > and \f$ \tau_m\f$. The hermitian conjugate \f$G^\ddagger\f$ is used.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param m
* > [int] Index of time \f$ \tau_m\f$ .
* @param i
* > [int] Index of time \f$ t_i\f$ .
* @param G_vt
* > [complex<T>] The right-mixing component (scalar GF).
* @param G
* > [herm_matrix] Contour function G
* @param Gcc
* > [herm_matrix] Hermitian conjugate \f$G^\ddagger\f$ of \f$G\f$.
*/
template <typename T>
inline void get_vt(const int m, const int i, std::complex<T> &G_vt, herm_matrix<T> &G, herm_matrix<T> &Gcc){
	assert(i <= G.nt() && m <= G.ntau());
	assert(i <= Gcc.nt() && m <= Gcc.ntau());
	assert(G.nt() == Gcc.nt());
	assert(G.ntau() == Gcc.ntau());
	assert(G.size1() == Gcc.size1());
	assert(G.size2() == Gcc.size2());
	assert(G.sig() == Gcc.sig());

	G_vt = *Gcc.tvptr(i, G.ntau() - m);
	G_vt = std::complex<T>(-G.sig(),0.0) * std::conj(G_vt);
}

/** \brief <b> Returns the right-mixing component of a general contour function at given times.  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Returns the right-mixing component \f$ G^\lceil(\tau_m,t_i) \f$ at given times \f$ t_i\f$
* > and \f$ \tau_m\f$. The hermitian conjugate \f$G^\ddagger\f$ is used.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param m
* > [int] Index of time \f$ \tau_m\f$ .
* @param i
* > [int] Index of time \f$ t_i\f$ .
* @param G_vt
* > [cdmatrix] The right-mixing component (returned to an eigen3 matrix).
* @param G
* > [herm_matrix] Contour function G
* @param Gcc
* > [herm_matrix] Hermitian conjugate \f$G^\ddagger\f$ of \f$G\f$.
*/
template <typename T>
void get_vt(const int m, const int i, cdmatrix &G_vt, herm_matrix<T> &G, herm_matrix<T> &Gcc){
	assert(i <= G.nt() && m <= G.ntau());
	assert(i <= Gcc.nt() && m <= Gcc.ntau());
	assert(G.nt() == Gcc.nt());
	assert(G.ntau() == Gcc.ntau());
	assert(G.size1() == Gcc.size1());
	assert(G.size2() == Gcc.size2());
	assert(G.sig() == Gcc.sig());
	int size1 = G.size1(), size2 = G.size2();
	std::complex<T> *vt;

	vt = Gcc.tvptr(i,G.ntau() - m);
	map_ptr2matrix(size1, size2, vt, G_vt);
	G_vt.adjointInPlace();
	G_vt = std::complex<T>(-G.sig(),0.0)* G_vt;
}




/** \brief <b> Returns the greatesr component of a general contour function at given times.  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Returns the greater component \f$ G^>(t_i,t_j) \f$ at given times \f$ t_i\f$
* > and \f$ t_j\f$ from \f$G^<(t_i,t_j)\f$ and \f$G^\mathrm{R}(t_i,t_j)\f$.
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > [int] Index of time \f$ t_i\f$ .
* @param j
* > [int] Index of time \f$ t_j\f$ .
* @param G_gtr
* > [complex<T>] The greater component (scalar GF).
* @param G
* > [herm_matrix] Contour function G
* @param Gcc
* > [herm_matrix] Hermitian conjugate \f$G^\ddagger\f$ of \f$G\f$.
*/
template <typename T>
inline void get_gtr(const int i, const int j, std::complex<T> &G_gtr, herm_matrix<T> &G, herm_matrix<T> &Gcc){
	std::complex<T> G_les, G_ret;

	get_les(i, j, G_les, G, Gcc);
	get_ret(i, j, G_ret, G, Gcc);
	G_gtr = G_les + G_ret;
}

/** \brief <b> Returns the greater component of a general contour function at given times.  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Returns the greater component \f$ G^>(t_i,t_j) \f$ at given times \f$ t_i\f$
* > and \f$ t_j\f$ from \f$G^<(t_i,t_j)\f$ and \f$G^\mathrm{R}(t_i,t_j)\f$.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > [int] Index of time \f$ t_i\f$ .
* @param j
* > [int] Index of time \f$ t_j\f$ .
* @param G_gtr
* > [cdmatrix] The greater component (returned to an eigen3 matrix).
* @param G
* > [herm_matrix] Contour function G
* @param Gcc
* > [herm_matrix] Hermitian conjugate \f$G^\ddagger\f$ of \f$G\f$.
*/
template <typename T>
void get_gtr(const int i, const int j, cdmatrix &G_gtr, herm_matrix<T> &G, herm_matrix<T> &Gcc){
	cdmatrix G_les, G_ret;

	get_les(i, j, G_les, G, Gcc);
	get_ret(i, j, G_ret, G, Gcc);
	G_gtr.noalias() = G_les + G_ret;

}







/** \brief <b> Returns the lesser component of a hermitian contour function at given times.  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Returns the lesser component \f$ G^<(t_i,t_j) \f$ at given times \f$ t_i\f$
* > and \f$ t_j\f$. Hermitian symmetry is assumed.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > [int] Index of time \f$ t_i\f$ .
* @param j
* > [int] Index of time \f$ t_j\f$ .
* @param G_les
* > [complex<T>] The lesser component (scalar).
* @param G
* > [herm_matrix] Contour function G
*/
template <typename T>
inline void get_les(const int i, const int j, std::complex<T> &G_les, herm_matrix<T> &G){
    get_les(i, j, G_les, G, G);
}


/** \brief <b> Returns the lesser component of a hermitian contour function at given times.  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Returns the lesser component \f$ G^<(t_i,t_j) \f$ at given times \f$ t_i\f$
* > and \f$ t_j\f$. Hermitian symmetry is assumed.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > [int] Index of time \f$ t_i\f$ .
* @param j
* > [int] Index of time \f$ t_j\f$ .
* @param G_les
* > [cdmatrix] The lesser component (returned to an eigen3 matrix).
* @param G
* > [herm_matrix] Contour function G
*/
template <typename T>
void get_les(const int i, const int j, cdmatrix &G_les, herm_matrix<T> &G){
	get_les(i, j, G_les, G, G);
}


/** \brief <b> Returns the retarded component of a hermitian contour function at given times.  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Returns the retarded component \f$ G^\mathrm{R}(t_i,t_j) \f$ at given times \f$ t_i\f$
* > and \f$ t_j\f$. We assume \f$G^\mathrm{R}(t_i,t_j)\f$ can be analytically continued
* > to \f$ j > i\f$, for which the hermitian conjugate \f$G^\ddagger\f$ is used.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > [int] Index of time \f$ t_i\f$ .
* @param j
* > [int] Index of time \f$ t_j\f$ .
* @param G_ret
* > [complex<T>] The retarded component (scalar GF).
* @param G
* > [herm_matrix] Contour function G
*/
template <typename T>
inline void get_ret(const int i, const int j, std::complex<T> &G_ret, herm_matrix<T> &G){
	get_ret(i, j, G_ret, G, G);
}

/** \brief <b> Returns the retarded component of a hermitian contour function at given times.  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Returns the retarded component \f$ G^\mathrm{R}(t_i,t_j) \f$ at given times \f$ t_i\f$
* > and \f$ t_j\f$. We assume \f$G^\mathrm{R}(t_i,t_j)\f$ can be analytically continued
* > to \f$ j > i\f$, for which the hermitian conjugate \f$G^\ddagger\f$ is used.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > [int] Index of time \f$ t_i\f$ .
* @param j
* > [int] Index of time \f$ t_j\f$ .
* @param G_ret
* > [cdmatrix] The retarded component (returned to an eigen3 matrix).
* @param G
* > [herm_matrix] Contour function G
*/
template <typename T>
void get_ret(const int i, const int j, cdmatrix &G_ret, herm_matrix<T> &G){
	get_ret(i, j, G_ret, G, G);
}



/** \brief <b> Returns the left-mixing component of a hermitian contour function at given times.  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Returns the left-mixing component \f$ G^\rceil(t_i,\tau_m) \f$ at given times \f$ t_i\f$
* > and \f$ \tau_m\f$.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > [int] Index of time \f$ t_i\f$ .
* @param m
* > [int] Index of time \f$ \tau_m\f$ .
* @param G_tv
* > [complex<T>] The left-mixing component (scalar GF).
* @param G
* > [herm_matrix] Contour function G
*/
template <typename T>
inline void get_tv(const int i, const int m, std::complex<T> &G_tv, herm_matrix<T> &G){
	get_tv(i, m, G_tv, G, G);
}

/** \brief <b> Returns the left-mixing component of a hermitian contour function at given times.  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Returns the left-mixing component \f$ G^\rceil(t_i,\tau_m) \f$ at given times \f$ t_i\f$
* > and \f$ \tau_m\f$.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > [int] Index of time \f$ t_i\f$ .
* @param m
* > [int] Index of time \f$ \tau_m\f$ .
* @param G_tv
* > [cdmatrix] The left-mixing component (returned to an eigen3 matrix).
* @param G
* > [herm_matrix] Contour function G
*/
template <typename T>
void get_tv(const int i, const int m, cdmatrix &G_tv, herm_matrix<T> &G){
	get_tv(i, m, G_tv, G, G);
}



/** \brief <b> Returns the right-mixing component of a hermitian contour function at given times.  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Returns the right-mixing component \f$ G^\lceil(\tau_m,t_i) \f$ at given times \f$ t_i\f$
* > and \f$ \tau_m\f$ assuming hermitian symmetry.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param m
* > [int] Index of time \f$ \tau_m\f$ .
* @param i
* > [int] Index of time \f$ t_i\f$ .
* @param G_vt
* > [complex<T>] The right-mixing component (scalar GF).
* @param G
* > [herm_matrix] Contour function G
*/
template <typename T>
inline void get_vt(const int m, const int i, std::complex<T> &G_vt, herm_matrix<T> &G){
	get_vt(m, i, G_vt, G, G);
}

/** \brief <b> Returns the right-mixing component of a hermitian contour function at given times.  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Returns the right-mixing component \f$ G^\lceil(\tau_m,t_i) \f$ at given times \f$ t_i\f$
* > and \f$ \tau_m\f$ assuming hermitian symmetry.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param m
* > [int] Index of time \f$ \tau_m\f$ .
* @param i
* > [int] Index of time \f$ t_i\f$ .
* @param G_vt
* > [cdmatrix] The right-mixing component (returned to an eigen3 matrix).
* @param G
* > [herm_matrix] Contour function G
*/
template <typename T>
void get_vt(const int m, const int i, cdmatrix &G_vt, herm_matrix<T> &G){
	get_vt(m, i, G_vt, G, G);
}




/** \brief <b> Returns the greatesr component of a hermitian contour function at given times.  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Returns the greater component \f$ G^>(t_i,t_j) \f$ at given times \f$ t_i\f$
* > and \f$ t_j\f$ from \f$G^<(t_i,t_j)\f$ and \f$G^\mathrm{R}(t_i,t_j)\f$.
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > [int] Index of time \f$ t_i\f$ .
* @param j
* > [int] Index of time \f$ t_j\f$ .
* @param G_gtr
* > [complex<T>] The greater component (scalar GF).
* @param G
* > [herm_matrix] Contour function G
*/
template <typename T>
inline void get_gtr(const int i, const int j, std::complex<T> &G_gtr, herm_matrix<T> &G){
	get_gtr(i, j, G_gtr, G, G);
}

/** \brief <b> Returns the greater component of a hermitian contour function at given times.  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Returns the greater component \f$ G^>(t_i,t_j) \f$ at given times \f$ t_i\f$
* > and \f$ t_j\f$ from \f$G^<(t_i,t_j)\f$ and \f$G^\mathrm{R}(t_i,t_j)\f$.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > [int] Index of time \f$ t_i\f$ .
* @param j
* > [int] Index of time \f$ t_j\f$ .
* @param G_gtr
* > [cdmatrix] The greater component (returned to an eigen3 matrix).
* @param G
* > [herm_matrix] Contour function G
*/
template <typename T>
void get_gtr(const int i, const int j, cdmatrix &G_gtr, herm_matrix<T> &G){
	get_gtr(i, j, G_gtr, G, G);
}
