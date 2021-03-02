#ifndef CNTR_HERM_MATRIX_HODLR_IMPL_H
#define CNTR_HERM_MATRIX_HODLR_IMPL_H

#include "cntr_herm_matrix_decl.hpp"
//#include "cntr_exception.hpp"
#include "cntr_elements.hpp"
#include "cntr_function_decl.hpp"
#include "cntr_herm_matrix_timestep_decl.hpp"
#include "cntr_herm_matrix_timestep_view_impl.hpp"
#include <Eigen/Core>
#include <Eigen/SVD>

namespace cntr {

/* #######################################################################################
#
#   CONSTRUCTION/DESTRUCTION
#
########################################################################################*/
template <typename T>
herm_matrix_hodlr<T>::herm_matrix_hodlr() {
    blkr1_=std::vector<int>();
    blkr2_=std::vector<int>();
    blkc1_=std::vector<int>();
    blkc2_=std::vector<int>();
    blklevel_=std::vector<int>();

    ranks_ret_=std::vector<int>();
    ranks_les_=std::vector<int>();
    ranks_tv_=std::vector<int>();

    ret_=std::vector<cntr::hodlr_box<T> >();
    les_=std::vector<cntr::hodlr_box<T>>();
    tv_=std::vector<cntr::hodlr_box<T> >();

    ret_dir_=std::vector<cdmatrix>();
    les_dir_=std::vector<cdmatrix>();
    tv_dir_=std::vector<cdmatrix>();

    nlvl=0;
    svdtol_=0.0;

    ntau_ = 0;
    nt_ = 0;
    size1_ = 0;
    size2_ = 0;
    element_size_ = 0;
    sig_ = -1;
}


template <typename T>
herm_matrix_hodlr<T>::~herm_matrix_hodlr() {

}

template <typename T>
herm_matrix_hodlr<T>::herm_matrix_hodlr(int nt,int nlvl,double svdtol) {
    assert(nt >= 0 && nlvl >= 0 && svdtol > 0);
    nt_=nt;
    nlvl_=nlvl;
    svdtol_=svdtol;
    std::cout << "Levels " << int(round(pow(2,nlvl)-1)) << std::endl;

    int tmp=(int)round(pow(2,nlvl)-1);
    blkr1_=std::vector<int>(tmp);
    blkr2_=std::vector<int>(tmp);
    blkc1_=std::vector<int>(tmp);
    blkc2_=std::vector<int>(tmp);
    blklevel_=std::vector<int>(tmp);

    int maxdir=0;  //Maximum number of directly-stored entries in a row
    getbkl_indices(maxdir);    
}




template <typename T>
void herm_matrix_hodlr<T>::getbkl_indices(int maxdir) {
    int tmp=(int)round(pow(2,nlvl_));
    int l;

    std::vector<int> is(tmp*2);
    //  First row indices of all blocks
    // is(i+1)=1+i*nt/(2^{nlvl})
    for(int i=0;i<tmp;i++){
        is[i]=(int)round(1+nt_*i/double(tmp));
        if(i>=1 && i<tmp) blkr1_[i-1]=is[i];
        // std::cout << "index " << i << " " << is[i] << " " << blkr1_[i] << std::endl;
    }
    for(int i=0; i<is.size(); ++i) std::cout << is[i] << ' ';
    std::cout << std::endl;
    for(int i=0; i<blkr1_.size(); ++i) std::cout << blkr1_[i] << ' ';
    std::cout << std::endl;

    // Levels of all blocks [closest to diagonal=1, \ldots]
    std::fill(blklevel_.begin(), blklevel_.end(), 1);
    for(int l=1;l<nlvl_;l++){
        // std::cout << int(pow(2,l)) << " " << int(pow(2,nlvl_)) << " " << int(pow(2,l)) << std::endl;
        for(int j=(int)round(pow(2,l));j<(int)round(pow(2,nlvl_));j=j+(int)round(pow(2,l))){
            blklevel_[j-1]=l+1;
            // std::cout << "Index: " << j << " " << l+1 << std::endl; 
        }
    }
    for(int i=0; i<blklevel_.size(); ++i) std::cout << blklevel_[i] << ' ';
    std::cout << std::endl;
    // First column indices of all blocks
    std::vector<int> il(nlvl_,1);
    std::fill(blkc1_.begin(), blkc1_.end(), 0);
    for(int i=0;i<tmp-1;i++){
        l=blklevel_[i];
        // std::cout <<" Values " << i << " " << l << " " << blkc1_[i]<< " " << (il[l]-1)*(int)round(pow(2,l))  << std::endl;
        blkc1_[i]=is[(il[l]-1)*(int)round(pow(2,l))];
        std::cout <<" Values " << i << " " << l << " " << blkc1_[i]<< " " << il[l] << " " << (il[l-1]-1)*(int)round(pow(2,l))  << std::endl;
        il[l]=il[l]+1;        
    }
    for(int i=0; i<blkc1_.size(); ++i) std::cout << blkc1_[i] << ' ';
    std::cout << std::endl;


}

template <typename T>
std::vector<int> herm_matrix_hodlr<T>::blkr1() const{
    return blkr1_;
}

template <typename T>
std::vector<int> herm_matrix_hodlr<T>::blkr2() const{
    return blkr2_;
}

template <typename T>
std::vector<int> herm_matrix_hodlr<T>::blkc1() const{
    return blkc1_;
}

template <typename T>
std::vector<int> herm_matrix_hodlr<T>::blkc2() const{
    return blkc2_;
}

template <typename T>
std::vector<int> herm_matrix_hodlr<T>::ranks_ret() const{
    return ranks_ret_;
}

template <typename T>
std::vector<int> herm_matrix_hodlr<T>::ranks_les() const{
    return ranks_les_;
}

template <typename T>
std::vector<int> herm_matrix_hodlr<T>::ranks_tv() const{
    return ranks_tv_;
}

template <typename T>
std::vector<cntr::hodlr_box<T> > herm_matrix_hodlr<T>::ret() const{
    return ret_;
}

template <typename T>
std::vector<cntr::hodlr_box<T> > herm_matrix_hodlr<T>::les() const{
    return les_;
}

template <typename T>
std::vector<cdmatrix> herm_matrix_hodlr<T>::ret_dir() const{
    return ret_dir_;
}

template <typename T>
std::vector<cdmatrix> herm_matrix_hodlr<T>::les_dir() const{
    return les_dir_;
}

template <typename T>
std::vector<cdmatrix> herm_matrix_hodlr<T>::tv_dir() const{
    return tv_dir_;
}


template <typename T>
hodlr_box<T>::hodlr_box() {
    rows_=0;
    cols_=0;
    epsrank_=0;
    svdtol_=0;
    S_=Eigen::VectorXd();
    U_=cdmatrix();
    V_=cdmatrix();
}

template <typename T>
hodlr_box<T>::hodlr_box(int rows,int cols,int epsrank,double eps) {
    rows_=rows;
    cols_=cols;
    epsrank_=epsrank;
    svdtol_=eps;
    S_=Eigen::VectorXd(epsrank);
    U_.resize(rows_,epsrank_); // Resize for the eps-rank
    V_.resize(epsrank_,cols_);
}



template <typename T>
hodlr_box<T>::~hodlr_box() {

}


template <typename T>
Eigen::VectorXd hodlr_box<T>::singular() const{
    return S_;
}


template <typename T>
cdmatrix hodlr_box<T>::U() const{
    return U_;
}


template <typename T>
cdmatrix hodlr_box<T>::V() const{
    return V_;
}

// TODO: remove at the end
template <typename T>
cdmatrix hodlr_box<T>::direct() const{
    return M_;
}



// /** \brief <b> Decompose the matrix into SVD </b>
// *
// * <!-- ====== DOCUMENTATION ====== -->
// *
// *  \par Purpose
// * <!-- ========= -->
// *
// * ***
// *
// * <!-- ARGUMENTS
// *      ========= -->
// *
// * @param M
// * > Number of time steps
// * @param ntau
// * > Number of points on Matsubara axis
// * @param size1
// * > Matrix rank of the contour function
// * @param sig
// * > Set `sig = -1` for fermions or `sig = +1` for bosons.
// */

template <typename T>
hodlr_box<T>::hodlr_box(cdmatrix &M,double svdtol) {
    assert(svdtol >= 0 && M.rows() >= 0 && M.rows() >= 0);
    svdtol_=svdtol;
    rows_=M.rows();
    cols_=M.cols();
    Eigen::BDCSVD<MatrixXcd> svd(M,Eigen::ComputeThinU | Eigen::ComputeThinV); 
    epsrank_=(svd.singularValues().array() >= svdtol_).count(); // determined eps rank
    U_.resize(M.rows(),epsrank_); // Resize for the eps-rank
    V_.resize(epsrank_,M.cols());

    S_=svd.singularValues().head(epsrank_);
    U_=svd.matrixU().block(0,0,M.rows(),epsrank_);
    V_=svd.matrixV().block(0,0,epsrank_,M.cols());
}



} // namespace cntr

#endif  // CNTR_HERM_MATRIX_HODLR_IMPL_H
