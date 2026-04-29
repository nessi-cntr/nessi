#ifndef CNTR_HERM_PSEUDO_IMPL_H
#define CNTR_HERM_PSEUDO_IMPL_H

#include "cntr_herm_pseudo_decl.hpp"
//#include "cntr_exception.hpp"
#include "cntr_elements.hpp"
#include "cntr_herm_matrix_timestep_decl.hpp"
#include "cntr_function_decl.hpp"

namespace cntr {

/* #######################################################################################
#
#   CONSTRUCTION/DESTRUCTION
#
########################################################################################*/
/// @private
template <typename T> herm_pseudo<T>::herm_pseudo() {
    les_ = 0;
    tv_ = 0;
    ret_ = 0;
    mat_ = 0;
    ntau_ = 0;
    nt_ = 0;
    size1_ = 0;
    size2_ = 0;
    element_size_ = 0;
    sig_ = -1;
}
/// @private
template <typename T> herm_pseudo<T>::~herm_pseudo() {
    delete[] les_;
    delete[] ret_;
    delete[] tv_;
    delete[] mat_;
}
/// @private
template <typename T> herm_pseudo<T>::herm_pseudo(int nt, int ntau, int size1, int sig) {
    assert(size1 >= 0 && nt >= -1 && sig * sig == 1 && ntau >= 0);
    nt_ = nt;
    ntau_ = ntau;
    sig_ = sig;
    size1_ = size1;
    size2_ = size1;
    element_size_ = size1 * size1;
    if (size1 > 0) {
        mat_ = new cplx[(ntau_ + 1) * element_size_];
        memset(mat_, 0, sizeof(cplx) * (ntau_ + 1) * element_size_);
    } else {
        mat_ = 0;
    }
    if (nt >= 0 && size1 > 0) {
        les_ = new cplx[((nt_ + 1) * (nt_ + 2)) / 2 * element_size_];
        ret_ = new cplx[((nt_ + 1) * (nt_ + 2)) / 2 * element_size_];
        tv_ = new cplx[(nt_ + 1) * (ntau_ + 1) * element_size_];
        memset(les_, 0, sizeof(cplx) * ((nt_ + 1) * (nt_ + 2)) / 2 * element_size_);
        memset(ret_, 0, sizeof(cplx) * ((nt_ + 1) * (nt_ + 2)) / 2 * element_size_);
        memset(tv_, 0, sizeof(cplx) * (nt_ + 1) * (ntau_ + 1) * element_size_);
    } else {
        les_ = 0;
        tv_ = 0;
        ret_ = 0;
    }
}
/// @private
template <typename T> unsigned long RAM_herm_pseudo(int nt, int ntau, int dim) {
    int element_size = dim * dim;
    size_t storage = 0;
    storage += sizeof(typename herm_pseudo<T>::cplx) * (ntau + 1) * element_size;
    storage +=
        sizeof(typename herm_pseudo<T>::cplx) * ((nt + 1) * (nt + 2)) / 2 * element_size;
    storage +=
        sizeof(typename herm_pseudo<T>::cplx) * ((nt + 1) * (nt + 2)) / 2 * element_size;
    storage += sizeof(typename herm_pseudo<T>::cplx) * (nt + 1) * (ntau + 1) * element_size;
    return storage;
}
/// @private
template <typename T> herm_pseudo<T>::herm_pseudo(const herm_pseudo &g) {
    nt_ = g.nt_;
    ntau_ = g.ntau_;
    sig_ = g.sig_;
    size1_ = g.size1_;
    size2_ = g.size1_;
    element_size_ = size1_ * size1_;
    if (size1_ > 0) {
        mat_ = new cplx[(ntau_ + 1) * element_size_];
        memcpy(mat_, g.mat_, sizeof(cplx) * (ntau_ + 1) * element_size_);
    } else {
        mat_ = 0;
    }
    if (nt_ >= 0 && size1_ > 0) {
        les_ = new cplx[((nt_ + 1) * (nt_ + 2)) / 2 * element_size_];
        ret_ = new cplx[((nt_ + 1) * (nt_ + 2)) / 2 * element_size_];
        tv_ = new cplx[(nt_ + 1) * (ntau_ + 1) * element_size_];
        memcpy(les_, g.les_, sizeof(cplx) * ((nt_ + 1) * (nt_ + 2)) / 2 * element_size_);
        memcpy(ret_, g.ret_, sizeof(cplx) * ((nt_ + 1) * (nt_ + 2)) / 2 * element_size_);
        memcpy(tv_, g.tv_, sizeof(cplx) * (nt_ + 1) * (ntau_ + 1) * element_size_);
    } else {
        les_ = 0;
        ret_ = 0;
        tv_ = 0;
    }
}
template <typename T> herm_pseudo<T> &herm_pseudo<T>::operator=(const herm_pseudo &g) {
    if (this == &g)
        return *this;
    sig_ = g.sig_;
    if (nt_ != g.nt_ || ntau_ != g.ntau_ || size1_ != g.size1_) {
        delete[] les_;
        delete[] ret_;
        delete[] tv_;
        delete[] mat_;
        nt_ = g.nt_;
        ntau_ = g.ntau_;
        size1_ = g.size1_;
        size2_ = g.size1_;
        element_size_ = size1_ * size1_;
        if (size1_ > 0) {
            mat_ = new cplx[(ntau_ + 1) * element_size_];
        } else {
            mat_ = 0;
        }
        if (size1_ > 0 && nt_ >= 0) {
            les_ = new cplx[((nt_ + 1) * (nt_ + 2)) / 2 * element_size_];
            ret_ = new cplx[((nt_ + 1) * (nt_ + 2)) / 2 * element_size_];
            tv_ = new cplx[(nt_ + 1) * (ntau_ + 1) * element_size_];
        } else {
            les_ = 0;
            ret_ = 0;
            tv_ = 0;
        }
    }
    if (size1_ > 0) {
        memcpy(mat_, g.mat_, sizeof(cplx) * (ntau_ + 1) * element_size_);
        if (nt_ >= 0) {
            memcpy(les_, g.les_, sizeof(cplx) * ((nt_ + 1) * (nt_ + 2)) / 2 * element_size_);
            memcpy(ret_, g.ret_, sizeof(cplx) * ((nt_ + 1) * (nt_ + 2)) / 2 * element_size_);
            memcpy(tv_, g.tv_, sizeof(cplx) * (nt_ + 1) * (ntau_ + 1) * element_size_);
        }
    }
    return *this;
}
/* #######################################################################################
#
#   RESIZE
#
########################################################################################*/
/// @private
template <typename T> void herm_pseudo<T>::resize_discard(int nt, int ntau, int size1) {
    assert(ntau >= 0 && nt >= -1 && size1 >= 0);
    delete[] les_;
    delete[] ret_;
    delete[] tv_;
    delete[] mat_;
    nt_ = nt;
    ntau_ = ntau;
    size1_ = size1;
    size2_ = size1;
    element_size_ = size1 * size1;
    if (size1 > 0) {
        mat_ = new cplx[(ntau_ + 1) * element_size_];
        memset(mat_, 0, sizeof(cplx) * (ntau_ + 1) * element_size_);
    } else {
        mat_ = 0;
    }
    if (nt_ >= 0 && size1_ > 0) {
        les_ = new cplx[((nt_ + 1) * (nt_ + 2)) / 2 * element_size_];
        ret_ = new cplx[((nt_ + 1) * (nt_ + 2)) / 2 * element_size_];
        tv_ = new cplx[(nt_ + 1) * (ntau_ + 1) * element_size_];
        memset(les_, 0, sizeof(cplx) * ((nt_ + 1) * (nt_ + 2)) / 2 * element_size_);
        memset(ret_, 0, sizeof(cplx) * ((nt_ + 1) * (nt_ + 2)) / 2 * element_size_);
        memset(tv_, 0, sizeof(cplx) * (nt_ + 1) * (ntau_ + 1) * element_size_);
    } else {
        les_ = 0;
        tv_ = 0;
        ret_ = 0;
    }
}
/// @private
template <typename T> void herm_pseudo<T>::resize_nt(int nt) {
    int nt1 = (nt_ > nt ? nt : nt_);
    cplx *ret, *les, *tv;
    assert(nt >= -1);
    nt_ = nt;
    if (size1_ == 0)
        return;
    if (nt_ >= 0) {
        les = new cplx[((nt_ + 1) * (nt_ + 2)) / 2 * element_size_];
        ret = new cplx[((nt_ + 1) * (nt_ + 2)) / 2 * element_size_];
        tv = new cplx[(nt_ + 1) * (ntau_ + 1) * element_size_];
        memset(les, 0, sizeof(cplx) * ((nt_ + 1) * (nt_ + 2)) / 2 * element_size_);
        memset(ret, 0, sizeof(cplx) * ((nt_ + 1) * (nt_ + 2)) / 2 * element_size_);
        memset(tv, 0, sizeof(cplx) * (nt_ + 1) * (ntau_ + 1) * element_size_);
        if (nt1 >= 0) {
            memcpy(les, les_, sizeof(cplx) * ((nt1 + 1) * (nt1 + 2)) / 2 * element_size_);
            memcpy(ret, ret_, sizeof(cplx) * ((nt1 + 1) * (nt1 + 2)) / 2 * element_size_);
            memcpy(tv, tv_, sizeof(cplx) * (nt1 + 1) * (ntau_ + 1) * element_size_);
        }
    } else {
        les = 0;
        ret = 0;
        tv = 0;
    }
    delete[] les_;
    delete[] ret_;
    delete[] tv_;
    les_ = les;
    ret_ = ret;
    tv_ = tv;
}
/// @private
template <typename T> void herm_pseudo<T>::resize(int nt, int ntau, int size1) {
    // std::cout  << "herm_pseudo<T>::" << __FUNCTION__ << " " << nt << " " << ntau << " " <<
    // size1 << std::endl;
    if (ntau == ntau_ && size1_ == size1)
        resize_nt(nt);
    else
        resize_discard(nt, ntau, size1);
}
/// @private
template <typename T> void herm_pseudo<T>::clear(void) {
    if (size1_ == 0)
        return;
    memset(mat_, 0, sizeof(cplx) * (ntau_ + 1) * element_size_);
    if (nt_ >= 0) {
        memset(les_, 0, sizeof(cplx) * ((nt_ + 2) * (nt_ + 1)) / 2 * element_size_);
        memset(ret_, 0, sizeof(cplx) * ((nt_ + 2) * (nt_ + 1)) / 2 * element_size_);
        memset(tv_, 0, sizeof(cplx) * (nt_ + 1) * (ntau_ + 1) * element_size_);
    }
}
/* #######################################################################################
#
#   RAW POINTERS TO ELEMENTS
#
########################################################################################*/
/// @private
template <typename T> inline std::complex<T> *herm_pseudo<T>::lesptr(int t, int t1) {
    assert(t >= 0 && t1 >= 0 && t <= t1 && t1 <= nt_);
    return les_ + ((t1 * (t1 + 1)) / 2 + t) * element_size_;
}
/// @private
template <typename T> inline std::complex<T> *herm_pseudo<T>::retptr(int t, int t1) {
    assert(t >= 0 && t1 >= 0 && t <= nt_ && t1 <= t);
    return ret_ + ((t * (t + 1)) / 2 + t1) * element_size_;
}
/// @private
template <typename T> inline std::complex<T> *herm_pseudo<T>::tvptr(int t, int tau) {
    assert(t >= 0 && tau >= 0 && t <= nt_ && tau <= ntau_);
    return tv_ + (t * (ntau_ + 1) + tau) * element_size_;
}
/// @private
template <typename T> inline std::complex<T> *herm_pseudo<T>::matptr(int tau) {
    assert(tau >= 0 && tau <= ntau_);
    return mat_ + tau * element_size_;
}
/* #######################################################################################
#
#   READING ELEMENTS TO ANY MATRIX TYPE
#   OR TO COMPLEX NUMBERS (then only the (0,0) element is addressed for dim>0)
#   note: these are not efficient, in particular the conjugation
#
########################################################################################*/
/// @private
#define herm_pseudo_READ_ELEMENT                                                            \
    {                                                                                       \
        int r, s, dim = size1_;                                                             \
        M.resize(dim, dim);                                                                 \
        for (r = 0; r < dim; r++)                                                           \
            for (s = 0; s < dim; s++)                                                       \
                M(r, s) = x[r * dim + s];                                                   \
    }
/// @private
#define herm_pseudo_READ_ELEMENT_MINUS_CONJ                                                 \
    {                                                                                       \
        cplx w;                                                                             \
        int r, s, dim = size1_;                                                             \
        M.resize(dim, dim);                                                                 \
        for (r = 0; r < dim; r++)                                                           \
            for (s = 0; s < dim; s++) {                                                     \
                w = x[s * dim + r];                                                         \
                M(r, s) = std::complex<T>(-w.real(), w.imag());                             \
            }                                                                               \
    }
/// @private
template <typename T>
template <class Matrix>
void herm_pseudo<T>::get_les(int i, int j, Matrix &M) {
    cplx *x;
    if (i <= j) {
        x = lesptr(i, j);
        herm_pseudo_READ_ELEMENT
    } else {
        x = lesptr(j, i);
        herm_pseudo_READ_ELEMENT_MINUS_CONJ
    }
}
/// @private
template <typename T>
template <class Matrix>
void herm_pseudo<T>::get_ret(int i, int j, Matrix &M) {
    cplx *x;
    if (i >= j) {
        x = retptr(i, j);
        herm_pseudo_READ_ELEMENT
    } else {
        x = retptr(j, i);
        herm_pseudo_READ_ELEMENT_MINUS_CONJ
    }
}
/// @private
template <typename T>
template <class Matrix>
void herm_pseudo<T>::get_tv(int i, int j, Matrix &M) {
    cplx *x = tvptr(i, j);
    herm_pseudo_READ_ELEMENT
}
/// @private
template <typename T>
template <class Matrix>
void herm_pseudo<T>::get_vt(int i, int j, Matrix &M) {
    cplx *x = tvptr(j, ntau_ - i);
    herm_pseudo_READ_ELEMENT_MINUS_CONJ if (sig_ == -1) M = -M;
}
/// @private
template <typename T>
template <class Matrix>
void herm_pseudo<T>::get_mat(int i, Matrix &M) {
    cplx *x = matptr(i);
    herm_pseudo_READ_ELEMENT
}
/// @private
template <typename T>
template <class Matrix>
void herm_pseudo<T>::get_matminus(int i, Matrix &M) {
    cplx *x = matptr(ntau_ - i);
    herm_pseudo_READ_ELEMENT if (sig_ == -1) M = -M;
}
/// @private
template <typename T>
template <class Matrix>
void herm_pseudo<T>::get_gtr(int i, int j, Matrix &M) {
    get_ret(i, j, M);
}
/// @private
template <typename T> inline void herm_pseudo<T>::get_ret(int i, int j, cplx &x) {
    if (i >= j)
        x = *retptr(i, j);
    else {
        x = *retptr(j, i);
        x = -std::conj(x);
    }
}
/// @private
template <typename T> inline void herm_pseudo<T>::get_les(int i, int j, cplx &x) {
    if (i <= j)
        x = *lesptr(i, j);
    else {
        x = *lesptr(j, i);
        x = -std::conj(x);
    }
}
/// @private
template <typename T> inline void herm_pseudo<T>::get_tv(int i, int j, cplx &x) {
    x = *tvptr(i, j);
}
/// @private
template <typename T> inline void herm_pseudo<T>::get_vt(int i, int j, cplx &x) {
    x = *tvptr(j, ntau_ - i);
    if (sig_ == -1)
        x = std::conj(x);
    else
        x = -std::conj(x);
}
/// @private
template <typename T> inline void herm_pseudo<T>::get_mat(int i, cplx &x) { x = *matptr(i); }
/// @private
template <typename T> inline void herm_pseudo<T>::get_matminus(int i, cplx &x) {
    x = *matptr(ntau_ - i);
    if (sig_ == -1)
        x = -x;
}
/// @private
template <typename T> inline void herm_pseudo<T>::get_gtr(int i, int j, cplx &x) {
    get_ret(i, j, x);
}
/* #######################################################################################
#
#   WRITING ELEMENTS FROM ANY MATRIX TYPE
#   OR FROM COMPLEX NUMBERS (then only the (0,0) element is addressed for dim>0)
#
########################################################################################*/
/// @private
#define herm_pseudo_SET_ELEMENT_MATRIX                                                      \
    {                                                                                       \
        int r, s, dim = size1_;                                                             \
        assert(M.rows() == dim && M.cols() == dim);                                         \
        for (r = 0; r < dim; r++)                                                           \
            for (s = 0; s < dim; s++)                                                       \
                x[r * dim + s] = M(r, s);                                                   \
    }
/// @private
template <typename T>
template <class Matrix>
void herm_pseudo<T>::set_ret(int i, int j, Matrix &M) {
    cplx *x = retptr(i, j);
    herm_pseudo_SET_ELEMENT_MATRIX
}
/// @private
template <typename T>
template <class Matrix>
void herm_pseudo<T>::set_les(int i, int j, Matrix &M) {
    cplx *x = lesptr(i, j);
    herm_pseudo_SET_ELEMENT_MATRIX
}
/// @private
template <typename T>
template <class Matrix>
void herm_pseudo<T>::set_tv(int i, int j, Matrix &M) {
    cplx *x = tvptr(i, j);
    herm_pseudo_SET_ELEMENT_MATRIX
}
/// @private
template <typename T>
template <class Matrix>
void herm_pseudo<T>::set_mat(int i, Matrix &M) {
    cplx *x = matptr(i);
    herm_pseudo_SET_ELEMENT_MATRIX
}
/// @private
template <typename T> inline void herm_pseudo<T>::set_les(int i, int j, cplx x) {
    *lesptr(i, j) = x;
}
/// @private
template <typename T> inline void herm_pseudo<T>::set_ret(int i, int j, cplx x) {
    *retptr(i, j) = x;
}
/// @private
template <typename T> inline void herm_pseudo<T>::set_tv(int i, int j, cplx x) {
    *tvptr(i, j) = x;
}
/// @private
template <typename T> inline void herm_pseudo<T>::set_mat(int i, cplx x) { *matptr(i) = x; }
/// @private
template <typename T> std::complex<T> herm_pseudo<T>::density_matrix(int tstp) {
    cplx x1;
    if (tstp == -1) {
        get_mat(ntau_, x1);
        return -x1;
    } else {
        get_les(tstp, tstp, x1);
        return std::complex<T>(0.0, sig_) * x1;
    }
}

// should work with Eigen
/// @private
template <typename T>
template <class Matrix>
void herm_pseudo<T>::density_matrix(int tstp, Matrix &M) {
    if (tstp == -1) {
        get_mat(ntau_, M);
        M *= (-1.0);
    } else {
        get_les(tstp, tstp, M);
        M *= std::complex<T>(0.0, 1.0 * sig_);
    }
}

/* #######################################################################################
#
#   INPUT/OUTPUT FROM/TO FILES
#
########################################################################################*/
/// @private
template <typename T> void herm_pseudo<T>::print_to_file(const char *file, int precision) {
    int i, j, l, sg = element_size_;
    std::ofstream out;
    out.open(file, std::ios::out);
    out.precision(precision);
    out << "# " << nt_ << " " << ntau_ << " " << size1_ << " "
        << " " << sig_ << std::endl;
    for (j = 0; j <= ntau_; j++) {
        out << "mat: " << j;
        for (l = 0; l < sg; l++)
            out << " " << matptr(j)[l].real() << " " << matptr(j)[l].imag();
        out << std::endl;
    }
    out << std::endl;
    if (nt_ >= 0) {
        for (i = 0; i <= nt_; i++) {
            for (j = 0; j <= i; j++) {
                out << "ret: " << i << " " << j;
                for (l = 0; l < sg; l++)
                    out << " " << retptr(i, j)[l].real() << " " << retptr(i, j)[l].imag();
                out << std::endl;
            }
            out << std::endl;
        }
        out << std::endl;
        for (i = 0; i <= nt_; i++) {
            for (j = 0; j <= ntau_; j++) {
                out << "tv: " << i << " " << j;
                for (l = 0; l < sg; l++)
                    out << " " << tvptr(i, j)[l].real() << " " << tvptr(i, j)[l].imag();
                out << std::endl;
            }
            out << std::endl;
        }
        out << std::endl;
        for (j = 0; j <= nt_; j++) {
            for (i = 0; i <= j; i++) {
                out << "les: " << i << " " << j;
                for (l = 0; l < sg; l++)
                    out << " " << lesptr(i, j)[l].real() << " " << lesptr(i, j)[l].imag();
                out << std::endl;
            }
            out << std::endl;
        }
        out << std::endl;
    }
    out.close();
}
/// @private
template <typename T> void herm_pseudo<T>::read_from_file(const char *file) {
    int i, n, m, j, l, size1, sg, sig;
    double real, imag;
    std::string s;
    std::ifstream out;
    out.open(file, std::ios::in);
    if (!(out >> s >> n >> m >> size1 >> sig)) {
        std::cerr << "read G from file " << file << " error in file" << std::endl;
        abort();
    }
    if (n > nt_ || m != ntau_ || size1 != size1_)
        resize(n, m, size1);
    sig_ = sig;
    sg = element_size_;
    for (j = 0; j <= ntau_; j++) {
        out >> s >> s;
        for (l = 0; l < sg; l++) {
            if (!(out >> real >> imag)) {
                std::cerr << "read G from file " << file << " error at mat (" << j << ")"
                          << std::endl;
                abort();
            }
            matptr(j)[l] = std::complex<T>(real, imag);
        }
    }
    if (n >= 0) {
        for (i = 0; i <= n; i++) {
            for (j = 0; j <= i; j++) {
                out >> s >> s >> s;
                for (l = 0; l < sg; l++) {
                    if (!(out >> real >> imag)) {
                        std::cerr << "read G from file " << file << " error at ret (" << i
                                  << "," << j << ")" << std::endl;
                        abort();
                    }
                    retptr(i, j)[l] = std::complex<T>(real, imag);
                }
            }
        }
        for (i = 0; i <= n; i++) {
            for (j = 0; j <= ntau_; j++) {
                out >> s >> s >> s;
                for (l = 0; l < sg; l++) {
                    if (!(out >> real >> imag)) {
                        std::cerr << "read G from file " << file << " error at tv (" << i
                                  << "," << j << ")" << std::endl;
                        abort();
                    }
                    tvptr(i, j)[l] = std::complex<T>(real, imag);
                }
            }
        }
        for (j = 0; j <= n; j++) {
            for (i = 0; i <= j; i++) {
                out >> s >> s >> s;
                for (l = 0; l < sg; l++) {
                    if (!(out >> real >> imag)) {
                        std::cerr << "read G from file " << file << " error at les (" << i
                                  << "," << j << ")" << std::endl;
                        abort();
                    }
                    lesptr(i, j)[l] = std::complex<T>(real, imag);
                }
            }
        }
    }
    out.close();
}
#if CNTR_USE_HDF5 == 1
// identical from herm_matrix<T>
/// @private
template <typename T> void herm_pseudo<T>::write_to_hdf5(hid_t group_id) {
    store_int_attribute_to_hid(group_id, std::string("ntau"), ntau_);
    store_int_attribute_to_hid(group_id, std::string("nt"), nt_);
    store_int_attribute_to_hid(group_id, std::string("sig"), sig_);
    store_int_attribute_to_hid(group_id, std::string("size1"), size1_);
    store_int_attribute_to_hid(group_id, std::string("size2"), size2_);
    store_int_attribute_to_hid(group_id, std::string("element_size"), element_size_);
    hsize_t len_shape = 3, shape[3];
    shape[1] = size1_;
    shape[2] = size2_;
    if (nt_ > -2) {
        shape[0] = ntau_ + 1;
        store_cplx_array_to_hid(group_id, std::string("mat"), matptr(0), shape, len_shape);
    }
    if (nt_ > -1) {
        shape[0] = ((nt_ + 1) * (nt_ + 2)) / 2;
        // CHECK: implement store_cplx_array_to_hid with template typename T
        store_cplx_array_to_hid(group_id, std::string("ret"), retptr(0, 0), shape,
                                len_shape);
        store_cplx_array_to_hid(group_id, std::string("les"), lesptr(0, 0), shape,
                                len_shape);
        shape[0] = (nt_ + 1) * (ntau_ + 1);
        store_cplx_array_to_hid(group_id, std::string("tv"), tvptr(0, 0), shape, len_shape);
    }
}
/// @private
template <typename T>
void herm_pseudo<T>::write_to_hdf5(hid_t group_id, const char *groupname) {
    hid_t sub_group_id = create_group(group_id, groupname);
    this->write_to_hdf5(sub_group_id);
    close_group(sub_group_id);
}
/// @private
template <typename T>
void herm_pseudo<T>::write_to_hdf5(const char *filename, const char *groupname) {
    hid_t file_id = open_hdf5_file(filename);
    this->write_to_hdf5(file_id, groupname);
    close_hdf5_file(file_id);
}
/// @private
template <typename T> void herm_pseudo<T>::read_from_hdf5(hid_t group_id) {
    // -- Read dimensions
    int nt = read_primitive_type<int>(group_id, "nt");
    int ntau = read_primitive_type<int>(group_id, "ntau");
    int sig = read_primitive_type<int>(group_id, "sig");
    int size1 = read_primitive_type<int>(group_id, "size1");
    // RESIZE G
    this->resize(nt, ntau, size1);
    sig_ = sig;
    if (nt > -2) {
        hsize_t mat_size = (ntau + 1) * element_size_;
        read_primitive_type_array(group_id, "mat", mat_size, matptr(0));
    }
    if (nt > -1) {
        hsize_t ret_size = ((nt + 1) * (nt + 2)) / 2 * element_size_;
        hsize_t les_size = ((nt + 1) * (nt + 2)) / 2 * element_size_;
        hsize_t tv_size = ((nt + 1) * (ntau + 1)) * element_size_;
        read_primitive_type_array(group_id, "ret", ret_size, retptr(0, 0));
        read_primitive_type_array(group_id, "les", les_size, lesptr(0, 0));
        read_primitive_type_array(group_id, "tv", tv_size, tvptr(0, 0));
    }
}
/// @private
template <typename T>
void herm_pseudo<T>::read_from_hdf5(hid_t group_id, const char *groupname) {
    hid_t sub_group_id = open_group(group_id, groupname);
    this->read_from_hdf5(sub_group_id);
    close_group(sub_group_id);
}
/// @private
template <typename T>
void herm_pseudo<T>::read_from_hdf5(const char *filename, const char *groupname) {
    hid_t file_id = read_hdf5_file(filename);
    this->read_from_hdf5(file_id, groupname);
    close_hdf5_file(file_id);
}
/// @private
template <typename T> void herm_pseudo<T>::read_from_hdf5(int nt1, hid_t group_id) {
    herm_pseudo<T> gtmp;
    gtmp.read_from_hdf5(group_id);
    assert(nt1 >= -1 && nt1 <= gtmp.nt());
    assert(nt1 >= -1 && nt1 <= nt_);
    assert(gtmp.size1() ==  size1_);
    assert(gtmp.element_size() == element_size_);
    assert(gtmp.ntau() == ntau_);

    for (int n = -1; n <= nt1; n++)
        this->set_timestep(n, gtmp);
}
/// @private
template <typename T>
void herm_pseudo<T>::read_from_hdf5(int nt1, hid_t group_id, const char *groupname) {
    hid_t sub_group_id = open_group(group_id, groupname);
    this->read_from_hdf5(nt1, sub_group_id);
    close_group(sub_group_id);
}
/// @private
template <typename T>
void herm_pseudo<T>::read_from_hdf5(int nt1, const char *filename, const char *groupname) {
    hid_t file_id = read_hdf5_file(filename);
    this->read_from_hdf5(nt1, file_id, groupname);
    close_hdf5_file(file_id);
}
#endif

/* #######################################################################################
#
#   SIMPLE OPERATIONS ON TIMESTEPS
#   NOTE: tstp IS A PHYSICAL TIME, tstp=-1 is the matsubara branch
#
########################################################################################*/
/// @private
template <typename T> void herm_pseudo<T>::set_timestep_zero(int tstp) {
    assert(tstp >= -1 && tstp <= nt_);
    if (tstp == -1) {
        memset(matptr(0), 0, sizeof(cplx) * (ntau_ + 1) * element_size_);
    } else {
        memset(retptr(tstp, 0), 0, sizeof(cplx) * (tstp + 1) * element_size_);
        memset(tvptr(tstp, 0), 0, sizeof(cplx) * (ntau_ + 1) * element_size_);
        memset(lesptr(0, tstp), 0, sizeof(cplx) * (tstp + 1) * element_size_);
    }
}
/// @private
template <typename T> void herm_pseudo<T>::set_timestep(int tstp, herm_pseudo<T> &g1) {
    assert(tstp >= -1 && tstp <= nt_ && tstp <= g1.nt());
    assert(g1.size1() == size1_);
    assert(g1.ntau() == ntau_);
    if (tstp == -1) {
        memcpy(mat_, g1.mat_, sizeof(cplx) * (ntau_ + 1) * element_size_);
    } else {
        memcpy(retptr(tstp, 0), g1.retptr(tstp, 0),
               sizeof(cplx) * (tstp + 1) * element_size_);
        memcpy(tvptr(tstp, 0), g1.tvptr(tstp, 0),
               sizeof(cplx) * (ntau_ + 1) * element_size_);
        memcpy(lesptr(0, tstp), g1.lesptr(0, tstp),
               sizeof(cplx) * (tstp + 1) * element_size_);
    }
}
/// @private
template <typename T>
void herm_pseudo<T>::set_timestep(int tstp, herm_matrix_timestep<T> &timestep) {
    cplx *x = timestep.data_;
    assert(tstp >= -1 && tstp <= nt_);
    assert(timestep.tstp_ == tstp && timestep.ntau_ == ntau_ && timestep.size1_ == size1_);
    if (tstp == -1) {
        memcpy(mat_, x, sizeof(cplx) * (ntau_ + 1) * element_size_);
    } else {
        memcpy(retptr(tstp, 0), x, sizeof(cplx) * (tstp + 1) * element_size_);
        memcpy(tvptr(tstp, 0), x + (tstp + 1) * element_size_,
               sizeof(cplx) * (ntau_ + 1) * element_size_);
        memcpy(lesptr(0, tstp), x + (tstp + 1 + ntau_ + 1) * element_size_,
               sizeof(cplx) * (tstp + 1) * element_size_);
    }
}
/// @private
template <typename T>
void herm_pseudo<T>::get_timestep(int tstp, herm_matrix_timestep<T> &timestep) {
    int len = (2 * (tstp + 1) + ntau_ + 1) * element_size_;
    cplx *x;
    assert(tstp >= -1 && tstp <= nt_);
    if (timestep.total_size_ < len)
        timestep.resize(tstp, ntau_, size1_);
    x = timestep.data_;
    timestep.tstp_ = tstp;
    timestep.ntau_ = ntau_;
    timestep.size1_ = size1_;
    if (tstp == -1) {
        memcpy(x, mat_, sizeof(cplx) * (ntau_ + 1) * element_size_);
    } else {
        memcpy(x, retptr(tstp, 0), sizeof(cplx) * (tstp + 1) * element_size_);
        memcpy(x + (tstp + 1) * element_size_, tvptr(tstp, 0),
               sizeof(cplx) * (ntau_ + 1) * element_size_);
        memcpy(x + (tstp + 1 + ntau_ + 1) * element_size_, lesptr(0, tstp),
               sizeof(cplx) * (tstp + 1) * element_size_);
    }
}
/// @private
#define herm_pseudo_INCR_TSTP                                                               \
    if (alpha == cplx(1.0, 0.0)) {                                                          \
        for (i = 0; i < len; i++)                                                           \
            x0[i] += x[i];                                                                  \
    } else {                                                                                \
        for (i = 0; i < len; i++)                                                           \
            x0[i] += alpha * x[i];                                                          \
    }
/// @private
template <typename T>
void herm_pseudo<T>::incr_timestep(int tstp, herm_matrix_timestep<T> &timestep,
                                   std::complex<T> alpha) {
    int i, len;
    cplx *x, *x0;
    assert(tstp >= -1 && tstp <= nt_);
    assert(timestep.tstp_ == tstp && timestep.ntau_ == ntau_ && timestep.size1_ == size1_);
    if (tstp == -1) {
        len = (ntau_ + 1) * element_size_;
        x0 = matptr(0);
        x = timestep.data_;
        herm_pseudo_INCR_TSTP
    } else {
        len = (tstp + 1) * element_size_;
        x0 = retptr(tstp, 0);
        x = timestep.data_;
        herm_pseudo_INCR_TSTP len = (ntau_ + 1) * element_size_;
        x0 = tvptr(tstp, 0);
        x = timestep.data_ + (tstp + 1) * element_size_;
        herm_pseudo_INCR_TSTP len = (tstp + 1) * element_size_;
        x0 = lesptr(0, tstp);
        x = timestep.data_ + (tstp + 1 + ntau_ + 1) * element_size_;
        herm_pseudo_INCR_TSTP
    }
}
#undef herm_pseudo_INCR_TSTP
/// @private
template <typename T> void herm_pseudo<T>::smul(int tstp, T weight) {
    int m;
    cplx *x0;
    assert(tstp >= -1 && tstp <= nt_);
    if (tstp == -1) {
        x0 = matptr(0);
        for (m = 0; m <= ntau_; m++) {
            element_smul<T, LARGESIZE>(size1_, x0 + m * element_size_, weight);
        }
    } else {
        x0 = retptr(tstp, 0);
        for (m = 0; m <= tstp; m++) {
            element_smul<T, LARGESIZE>(size1_, x0 + m * element_size_, weight);
        }
        x0 = tvptr(tstp, 0);
        for (m = 0; m <= ntau_; m++) {
            element_smul<T, LARGESIZE>(size1_, x0 + m * element_size_, weight);
        }
        x0 = lesptr(0, tstp);
        for (m = 0; m <= tstp; m++) {
            element_smul<T, LARGESIZE>(size1_, x0 + m * element_size_, weight);
        }
    }
}
/// @private
template <typename T> void herm_pseudo<T>::smul(int tstp, std::complex<T> weight) {
    int m;
    cplx *x0;
    assert(tstp >= -1 && tstp <= nt_);
    if (tstp == -1) {
        x0 = matptr(0);
        for (m = 0; m <= ntau_; m++) {
            element_smul<T, LARGESIZE>(size1_, x0 + m * element_size_, weight);
        }
    } else {
        x0 = retptr(tstp, 0);
        for (m = 0; m <= tstp; m++) {
            element_smul<T, LARGESIZE>(size1_, x0 + m * element_size_, weight);
        }
        x0 = tvptr(tstp, 0);
        for (m = 0; m <= ntau_; m++) {
            element_smul<T, LARGESIZE>(size1_, x0 + m * element_size_, weight);
        }
        x0 = lesptr(0, tstp);
        for (m = 0; m <= tstp; m++) {
            element_smul<T, LARGESIZE>(size1_, x0 + m * element_size_, weight);
        }
    }
}

///////////////////////////////////////////////////////////////
// multiply timestep with function ... same as for herm_matrix
// G(t,t') ==> F(t)G(t,t')   ... ft+t*element_size_ points to F(t)
/// @private
template <typename T>
void herm_pseudo<T>::left_multiply(int tstp, std::complex<T> *f0, std::complex<T> *ft,
                                   T weight) {
    int m;
    cplx *xtemp, *ftemp, *x0;
    xtemp = new cplx[element_size_];
    assert(tstp >= -1 && tstp <= nt_);
    if (tstp == -1) {
        x0 = matptr(0);
        for (m = 0; m <= ntau_; m++) {
            element_mult<T, LARGESIZE>(size1_, xtemp, f0, x0 + m * element_size_);
            element_smul<T, LARGESIZE>(size1_, xtemp, weight);
            element_set<T, LARGESIZE>(size1_, x0 + m * element_size_, xtemp);
        }
    } else {
        ftemp = ft + tstp * element_size_;
        x0 = retptr(tstp, 0);
        for (m = 0; m <= tstp; m++) {
            element_mult<T, LARGESIZE>(size1_, xtemp, ftemp, x0 + m * element_size_);
            element_smul<T, LARGESIZE>(size1_, xtemp, weight);
            element_set<T, LARGESIZE>(size1_, x0 + m * element_size_, xtemp);
        }
        x0 = tvptr(tstp, 0);
        for (m = 0; m <= ntau_; m++) {
            element_mult<T, LARGESIZE>(size1_, xtemp, ftemp, x0 + m * element_size_);
            element_smul<T, LARGESIZE>(size1_, xtemp, weight);
            element_set<T, LARGESIZE>(size1_, x0 + m * element_size_, xtemp);
        }
        x0 = lesptr(0, tstp);
        for (m = 0; m <= tstp; m++) {
            element_mult<T, LARGESIZE>(size1_, xtemp, ft + m * element_size_,
                                       x0 + m * element_size_);
            element_smul<T, LARGESIZE>(size1_, xtemp, weight);
            element_set<T, LARGESIZE>(size1_, x0 + m * element_size_, xtemp);
        }
    }
    delete[] xtemp;
}
/// @private
template <typename T>
void herm_pseudo<T>::left_multiply(int tstp, function<T> &ft, T weight) {
    assert(tstp >= -1 && tstp <= nt_ && ft.nt() >= tstp && ft.size1() == size1_ &&
           ft.size2() == size2_);
    if (tstp >= 0)
        left_multiply(tstp, ft.ptr(-1), ft.ptr(0), weight);
    else
        left_multiply(tstp, ft.ptr(-1), 0, weight);
}
// G(t,t') ==> F(t)G(t,t')   ... ft+t*element_size_ points to F(t)
/// @private
template <typename T>
void herm_pseudo<T>::right_multiply(int tstp, std::complex<T> *f0, std::complex<T> *ft,
                                    T weight) {
    int m;
    cplx *xtemp, *ftemp, *x0;
    xtemp = new cplx[element_size_];
    assert(tstp >= -1 && tstp <= nt_);
    if (tstp == -1) {
        x0 = matptr(0);
        for (m = 0; m <= ntau_; m++) {
            element_mult<T, LARGESIZE>(size1_, xtemp, x0 + m * element_size_, f0);
            element_smul<T, LARGESIZE>(size1_, xtemp, weight);
            element_set<T, LARGESIZE>(size1_, x0 + m * element_size_, xtemp);
        }
    } else {
        x0 = retptr(tstp, 0);
        for (m = 0; m <= tstp; m++) {
            element_mult<T, LARGESIZE>(size1_, xtemp, x0 + m * element_size_,
                                       ft + m * element_size_);
            element_smul<T, LARGESIZE>(size1_, xtemp, weight);
            element_set<T, LARGESIZE>(size1_, x0 + m * element_size_, xtemp);
        }
        x0 = tvptr(tstp, 0);
        for (m = 0; m <= ntau_; m++) {
            element_mult<T, LARGESIZE>(size1_, xtemp, x0 + m * element_size_, f0);
            element_smul<T, LARGESIZE>(size1_, xtemp, weight);
            element_set<T, LARGESIZE>(size1_, x0 + m * element_size_, xtemp);
        }
        ftemp = ft + tstp * element_size_;
        x0 = lesptr(0, tstp);
        for (m = 0; m <= tstp; m++) {
            element_mult<T, LARGESIZE>(size1_, xtemp, x0 + m * element_size_, ftemp);
            element_smul<T, LARGESIZE>(size1_, xtemp, weight);
            element_set<T, LARGESIZE>(size1_, x0 + m * element_size_, xtemp);
        }
    }
    delete[] xtemp;
}
/// @private
template <typename T>
void herm_pseudo<T>::right_multiply(int tstp, function<T> &ft, T weight) {
    assert(tstp >= -1 && tstp <= nt_ && ft.nt() >= tstp && ft.size1() == size1_ &&
           ft.size2() == size2_);
    if (tstp >= 0)
        right_multiply(tstp, ft.ptr(-1), ft.ptr(0), weight);
    else
        right_multiply(tstp, ft.ptr(-1), 0, weight);
}

#if 0

template<typename T> void herm_pseudo<T>::get_timestep(int tstp,std::vector<cplx> &data){
	int len=2*(tstp+1)+ntau_,len1,j;
	cplx *x;
	len1=len*element_size_;
	assert(tstp>=-1 && tstp<=nt_);
	if((int)data.size()<len1) data.resize(len1);
	j=(tstp>=0 ? IDXLOWER(tstp) : IDXIMAG(0));
	x=elementptr(j,j);
	for(i=0;i<len1;i++) data[i]=x[i];
}
template<typename T,class Matrix> void herm_pseudo<T>::get_timestep(int tstp,std::vector<Matrix> &data){
	int len=2*(tstp+1)+ntau_,j;
	cplx *x;
	assert(tstp>=-1 && tstp<=nt_);
	if((int)data.size()<len) data.resize(len);
	j=(tstp>=0 ? IDXLOWER(tstp) : IDXIMAG(0));
	x=elementptr(j,j);
	for(i=0;i<len;i++){
	    if(data[i].size1()!=size1_ || data[i].size2()!=size1_) data[i].resize(size1_,size1_);

		data[i]=x[i*element_size_];
	}
}
template<typename T> void herm_pseudo<T>::set_timestep(int tstp,std::vector<cplx> &data){
	int len=2*(tstp+1)+ntau_;
	assert((int)data.size()>=len*element_size_);
	if(tstp>=0){
		j=IDXLOWER(tstp);
		for(i=0;i<len;i++){
		  i1=nt_-tstp+i;
		}
	}else{
	}
}
template<typename T> void herm_pseudo<T>::set_timestep_to_zero(int tstp);
template<typename T> void herm_pseudo<T>::smul_timestep(int tstp,cplx lambda);
template<typename T> void herm_pseudo<T>::increment_timestep(int tstp,std::vector<cplx> &data,cplx alpha=1.0);
template<typename T> void herm_pseudo<T>::set_timestep(int tstp,const herm_pseudo<T> &g1);
template<typename T> void herm_pseudo<T>::set_timestep_pseudo(int tstp,const matrix<T> &g1);
template<class Matrix> void get_timestep(int tstp,std::vector<Matrix> &data);
template<class Matrix> void set_timestep(int tstp,std::vector<Matrix> &data);

/* #######################################################################################
#
#   MANIPULATING TIMESTEPS
#
########################################################################################*/
template <class GG>
void set_timestep_pseudo(int tstp,GG &G,std::complex<double> *timestep){
  int sg=G.element_size(),ntau=G.ntau(),m,l=0;
  std::complex<double> *x1;
  x1=timestep;
  if(tstp>=0){
    for(m=0;m<=tstp;m++){
	  // ret(tstp,m) = gtr(tstp,m)
	  // = -G(m-,tstp-)^*
	  G.element_set(G.retptr(tstp,tstp-m),x1,G);
	  G.element_smul(G.retptr(tstp,tstp-m),-1);
	  G.element_conj(G.retptr(tstp,tstp-m));
	  x1+=sg;l++;
	}
	for(m=0;m<=ntau;m++){
	  // x1 points to Gvt(m,tstp) = -bosefermi *Gtv(tstp,ntau-m)^*
	  G.element_set(G.tvptr(tstp,ntau-m),x1,G);
	  G.element_conj(G.tvptr(tstp,ntau-m));
	  G.element_smul(G.tvptr(tstp,ntau-m),-G.sig());
	  x1+=sg;l++;
	}
	for(m=0;m<=tstp;m++){
	  // x1 points to G(m+,tstp+)=les(m,tstp)
	  G.element_set(G.lesptr(m,tstp),x1,G);
	  x1+=sg;l++;
    }
  }else{
    for(m=0;m<=ntau;m++){
	  // x1 points to Gmat(m,0) = i*Gmat(m)
	  G.element_set(G.matptr(m),x1,G);
	  G.element_smul(G.matptr(m),std::complex<double>(0,-1));
	  x1+=sg;l++;
	}
  }
}
template <class GG>
void get_timestep_pseudo(int tstp,GG &G,std::complex<double> *timestep){
int sg=G.element_size(),ntau=G.ntau(),m,l=0;
  std::complex<double> *x1;
  x1=timestep;
  //std::cout << "get tsp " << tstp << std::endl;
  if(tstp>=0){
    for(m=0;m<=tstp;m++){
	  // ret(tstp,tstp-m) = gtr(tstp,tstp-m)
	  // = -G^gtr(tstp-m,tstp)^*
	  G.element_set(x1,G.retptr(tstp,tstp-m),G);
	  G.element_smul(x1,-1);
	  G.element_conj(x1);
	  x1+=sg;l++;
	}
	for(m=0;m<=ntau;m++){
	  // x1 points to Gvt(m,tstp) = -bosefermi *Gtv(tstp,ntau-m)^*
	  G.element_set(x1,G.tvptr(tstp,ntau-m),G);
	  G.element_conj(x1);
	  G.element_smul(x1,-G.sig());
	  x1+=sg;l++;
	}
	for(m=0;m<=tstp;m++){
	  // x1 points to G(m+,tstp+)=les(m,tstp)
	  G.element_set(x1,G.lesptr(m,tstp),G);
	  x1+=sg;l++;
    }
  }else{
    for(m=0;m<=ntau;m++){
	  // x1 points to Gmat(m,0) = i*Gmat(m)
	  G.element_set(x1,G.matptr(m),G);
	  G.element_smul(x1,std::complex<double>(0,1));
	  x1+=sg;l++;
	}
  }
}
template <typename T,class Matrix>
void set_timestep_matrix_pseudo(int tstp,matrix<T> &G,std::vector<Matrix> &data){
  int i,l,m,sg,dim,ntau=G.ntau(),nmax=(tstp<0 ? ntau+1 : ntau+2*tstp+3);
  std::complex<T> *data1;
  assert((int)data.size()>=nmax);
  dim=G.size1();
  sg=dim*dim;
  data1=new std::complex<T> [nmax*sg];
  for(i=0;i<nmax;i++) for(l=0;l<dim;l++) for(m=0;m<dim;m++) data1[i*sg+l*dim+m]=data[i](l,m);
  set_timestep_pseudo(tstp,G,data1);
  delete [] data1;
}
template <typename T,class Matrix>
void increment_timestep_matrix_pseudo(int tstp,matrix<T> &G,std::complex<T> alpha,std::vector<Matrix> &data){
  int i,l,m,sg,dim,ntau=G.ntau(),nmax=(tstp<0 ? ntau+1 : ntau+2*tstp+3);
  std::complex<T> *data1;
  assert((int)data.size()>=nmax);
  dim=G.size1();
  sg=dim*dim;
  data1=new std::complex<T> [nmax*sg];
  get_timestep_pseudo(tstp,G,data1);
  if(alpha.real()==1 && alpha.imag()==0) for(i=0;i<nmax;i++) for(l=0;l<dim;l++) for(m=0;m<dim;m++) data1[i*sg+l*dim+m] += data[i](l,m);
  else for(i=0;i<nmax;i++) for(l=0;l<dim;l++) for(m=0;m<dim;m++) data1[i*sg+l*dim+m] += alpha*data[i](l,m);
  set_timestep_pseudo(tstp,G,data1);
  delete [] data1;
}
template <typename T,class Matrix>
void get_timestep_matrix_pseudo(int tstp,matrix<T> &G,std::vector<Matrix> &data){
  int i,l,m,sg,dim,ntau=G.ntau(),nmax=(tstp<0 ? ntau+1 : ntau+2*tstp+3);
  std::complex<T> *data1;
  assert((int)data.size()>=nmax);
  dim=G.size1();
  sg=dim*dim;
  data1=new std::complex<T> [nmax*sg];
  get_timestep_pseudo(tstp,G,data1);
  for(i=0;i<nmax;i++) for(l=0;l<dim;l++) for(m=0;m<dim;m++) data[i](l,m)=data1[i*sg+l*dim+m];
  delete [] data1;
}

#endif // 0

} // namespace cntr

#endif  // CNTR_HERM_PSEUDO_IMPL_H
