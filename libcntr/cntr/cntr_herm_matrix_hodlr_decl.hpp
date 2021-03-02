#ifndef CNTR_HERM_MATRIX_HODLR_DECL_H
#define CNTR_HERM_MATRIX_HODLR_DECL_H

#include "cntr_global_settings.hpp"

namespace cntr {

template <typename T> class function;
template <typename T> class herm_matrix_timestep;


template <typename T>
/** \brief <b> Class `herm_matrix_hodlr` for two-time contour objects \f$ C(t,t') \f$
 * compressed with the hodlr structure and with hermitian symmetry .</b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 *  This class contains the data structures for representing bosonic (\f$\eta = +1 \f$)
 *  or fermionic (\f$\eta = -1 \f$) or two-time contour
 *  functions \f$ C(t,t') \f$ in the hierarchical off-diagonal low-rank structure.
 *  The class `herm_matrix_hodlr` stores TODO
 *   -  TODO
 *
 * The contour function \f$ C(t,t') \f$ can be of scalar type or matrix-valued.
 * Square-matrix-type contour functions are fully supported, while non-square-
 * matrix contour functions are under development.
 *
 *  If `nt = 0`, only the Matsubara component is stored.
 *
 */
class hodlr_box {
  public:
    typedef std::complex<T> cplx;
    typedef T scalar_type;

    /* construction, destruction */
    hodlr_box();
    hodlr_box(cdmatrix &M,double eps);
    hodlr_box(int rows,int cols,int epsrank,double eps);
    ~hodlr_box();
    int cols(void) const {return cols_;};
    int rows(void) const {return rows_;};
    int epsrank(void) const {return epsrank_;};
    double svdtol(void) const {return svdtol_;};
    Eigen::VectorXd singular(void) const;
    cdmatrix U(void) const;
    cdmatrix V(void) const;
    cdmatrix direct(void) const;
// // MPI UTILS
// #if CNTR_USE_MPI == 1
//     void Reduce_timestep(int tstp, int root);
//     void Bcast_timestep(int tstp, int root);
//     void Send_timestep(int tstp, int dest, int tag);
//     void Recv_timestep(int tstp, int root, int tag);
// #endif
//   private:
//     int les_offset(int t, int t1) const;
//     int ret_offset(int t, int t1) const;
//     int tv_offset(int t, int tau) const;
//     int mat_offset(int tau) const;

  private:
    /// @private
    /** \brief <b> Matrices corresponding to the left (U) and right (V) eigenvectors */
    cdmatrix U_,V_;
    /// @private
    /** \brief <b> Vector of diagonal eigenvalues </b> */
    Eigen::VectorXd S_;
    /// @private
    /** \brief <b> SVD tolerance </b> */
    double svdtol_; 
    /// @private
    /** \brief <b> Epsilon rank \f$ </b> */
    int epsrank_;
    /// @private
    /** \brief <b>  Number of rows in block.</b> */
    int rows_;
    /// @private
    /** \brief <b> Number of columns in block.</b> */
    int cols_;
    /// @private
    /** \brief <b> Direct storage - only for debuggin TODO: Remove at the end</b> */
    int M_;

};

template <typename T>
/** \brief <b> Class `herm_matrix_hodlr` for two-time contour objects \f$ C(t,t') \f$
 * compressed with the hodlr structure and with hermitian symmetry .</b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 *  This class contains the data structures for representing bosonic (\f$\eta = +1 \f$)
 *  or fermionic (\f$\eta = -1 \f$) or two-time contour
 *  functions \f$ C(t,t') \f$ in the hierarchical off-diagonal low-rank structure.
 *  The class `herm_matrix_hodlr` stores TODO
 *   -  TODO
 *
 * The contour function \f$ C(t,t') \f$ can be of scalar type or matrix-valued.
 * Square-matrix-type contour functions are fully supported, while non-square-
 * matrix contour functions are under development.
 *
 *  If `nt = 0`, only the Matsubara component is stored.
 *
 */
class herm_matrix_hodlr {
  public:
    typedef std::complex<T> cplx;
    typedef T scalar_type;

    /* construction, destruction */
    herm_matrix_hodlr();
    ~herm_matrix_hodlr();
    herm_matrix_hodlr(int nt,int nlvl,double svdtol);

    void getbkl_indices(int maxdir);

    int nlvl(void) const {return nlvl_;};
    int nt(void) const {return nt_;};
    int ntau(void) const {return ntau_;};
    int size1(void) const {return size1_;};
    int size2(void) const {return size2_;};
    int element_size(void) const {return element_size_;};
    int sig(void) const {return sig_;};
    double svdtol(void) const {return svdtol_;};
    std::vector<int> blkr1(void) const;
    std::vector<int> blkr2(void) const;
    std::vector<int> blkc1(void) const;
    std::vector<int> blkc2(void) const;
    std::vector<int> blklevel(void) const;
    std::vector<int> ranks_ret(void) const;
    std::vector<int> ranks_les(void) const;
    std::vector<int> ranks_tv(void) const;
    std::vector<cntr::hodlr_box<T> > ret(void) const;
    std::vector<cntr::hodlr_box<T> > les(void) const;
    std::vector<cntr::hodlr_box<T> > tv(void) const;
    std::vector<cdmatrix> ret_dir(void) const;
    std::vector<cdmatrix> les_dir(void) const;
    std::vector<cdmatrix> tv_dir(void) const;

// // MPI UTILS
// #if CNTR_USE_MPI == 1
//     void Reduce_timestep(int tstp, int root);
//     void Bcast_timestep(int tstp, int root);
//     void Send_timestep(int tstp, int dest, int tag);
//     void Recv_timestep(int tstp, int root, int tag);
// #endif
//   private:
//     int les_offset(int t, int t1) const;
//     int ret_offset(int t, int t1) const;
//     int tv_offset(int t, int tau) const;
//     int mat_offset(int tau) const;

  private:
    /// @private
    /** \brief <b> First rows in blocks </b> */
    std::vector<int> blkr1_;
    /// @private
    /** \brief <b> Last rows in blocks </b> */
    std::vector<int> blkr2_; 
    /// @private
    /** \brief <b> First cols in blocks </b> */
    std::vector<int> blkc1_;
    /// @private
    /** \brief <b> Last cols in blocks </b> */
    std::vector<int> blkc2_;
    /// @private
    /** \brief <b> Levels of blocks </b> */
    std::vector<int> blklevel_;
    /// @private
    /** \brief <b> Ranks of blocks for retarded component  </b> */
    std::vector<int> ranks_ret_;
    /// @private
    /** \brief <b> Ranks of blocks for les component  </b> */
    std::vector<int> ranks_les_;
    /// @private
    /** \brief <b> Ranks of blocks for tv component  </b> */
    std::vector<int> ranks_tv_;
    /// @private
    /** \brief <b> Retarded blocks  </b> */
    std::vector<cntr::hodlr_box<T> > ret_;
    /// @private
    /** \brief <b> Lesser blocks   </b> */
    std::vector<cntr::hodlr_box<T> > les_;
    /// @private
    /** \brief <b> tv blocks   </b> */
    std::vector<cntr::hodlr_box<T> > tv_;
    /// @private
    /** \brief <b> Directly stored retarded components  </b> */
    std::vector<cdmatrix> ret_dir_;
    /// @private
    /** \brief <b> Lesser blocks   </b> */
    std::vector<cdmatrix> les_dir_;
    /// @private
    /** \brief <b> tv blocks   </b> */
    std::vector<cdmatrix> tv_dir_;
    /// @private
    /** \brief <b> mat TODO   </b> */
    // std::vector<cdmatrix> mat_; TODO
    /// @private
    /** \brief <b> Number of levels in hodlr.</b> */
    int nlvl_;
    /// @private
    /** \brief <b> SVD tolerance.</b> */
    double svdtol_;
    /// @private
    /** \brief <b> Maximum number of the time steps.</b> */
    int nt_;
    /// @private
    /** \brief <b> Number of the time grids on the Matsubara axis.</b> */
    int ntau_;
    /// @private
    /** \brief <b> Number of the colums in the Matrix form.</b> */
    int size1_;
    /// @private
    /** \brief <b> Number of the rows in the Matrix form.</b> */
    int size2_;
    /// @private
    /** \brief <b> Size of the Matrix form; size1*size2. </b> */
    int element_size_;
    /// @private
    /** \brief <b> Bose = +1, Fermi =-1. </b> */
    int sig_; // Bose = +1, Fermi =-1
};




}  // namespace cntr

#endif  // CNTR_HERM_MATRIX_HODLR_DECL_H
