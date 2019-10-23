#ifndef CNTR_FUNCTION_IMPL_H
#define CNTR_FUNCTION_IMPL_H

#include "cntr_function_decl.hpp"
#include "linalg.hpp"
//#include "cntr_exception.hpp"
#include "cntr_elements.hpp"

namespace cntr {

template <typename T>
function<T>::function() {
    data_ = 0;
    nt_ = -2;
    size1_ = 0;
    size2_ = 0;
    element_size_ = 0;
}
template <typename T>
function<T>::~function() {
    if (data_ != 0)
        delete[] data_;
}
/** \brief <b> Initializes the `function` class for a square-matrix-valued function of time.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Initializes the `function` class for a square-matrix-valued function of time.
* > `nt = 0` leads to only one element (with time = -1), which is the value of the function
* > on the Matsubara axis (initial equilibirum)
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param nt
* > Number of time steps
* @param size1
* > Leading dimension of matrix
*/

template <typename T>
function<T>::function(int nt, int size1) {
    int len = (nt + 2) * size1 * size1;
    assert(size1 >= 0 && nt >= -1);
    if (len == 0)
        data_ = 0;
    else {
        data_ = new cplx[len];
    }
    size1_ = size1;
    size2_ = size1;
    element_size_ = size1 * size1_;
    nt_ = nt;
    total_size_ = (nt_ + 2) * size1_ * size2_;
}

/** \brief <b> Initializes the `function` class for a non-square-matrix-valued function of time.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Initializes the `function` class for a non-square-matrix valued function of time.
* > `nt = 0` leads to only one element (with time = -1), which is the value of function
* > on the Matsubara axis (initial equilibirum)
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param nt
* > Number of time steps
* @param size1
* > Leading dimension of matrix
* @param size2
* > Second dimension of matrix
*/
template <typename T> function<T>::function(int nt,int size1,int size2){
   int len=(nt+2)*size1*size2;
   assert(size1>=0 && nt>=-1);
   if(len==0) data_=0;
   else{
      data_ = new cplx [len];
	  memset(data_, 0, sizeof(cplx)*len);
   }
   size1_=size1;
   size2_=size2;
   element_size_=size1*size2_;
   nt_=nt;
   total_size_ = len;
}
/** \brief <b> Initializes the `function` class from an existing function object</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Initializes the `function` class from an existing `function` object. It copies all
* > data from the existing `function` object.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param g
* > the `function` object one wants to copy
*/

template <typename T>
function<T>::function(const function &g) {
    int len;
    nt_ = g.nt_;
    size1_ = g.size1_;
    size2_ = g.size2_;
    element_size_ = size1_ * size2_;
    len = (nt_ + 2) * size1_ * size2_;
    if (len > 0) {
	total_size_ = len;
        data_ = new cplx[len];
        memcpy(data_, g.data_, sizeof(cplx) * len);
    } else {
        data_ = 0;
    }
}
#if __cplusplus >= 201103L
/** \brief <b> Initializes the `function` class from an existing `function` object(right-value reference)</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Initializes the `function` class from an existing `function` object. It copies all
* > data from the existing `function` object. The method is specified for right value reference when
* > C++11 is used.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param g
* > the `function` object one wants to copy.
*/

template <typename T>
function<T>::function(function &&g) noexcept
    : data_(g.data_),
      nt_(g.nt_),
      size1_(g.size1_),
      size2_(g.size2_),
      element_size_(g.element_size_) {
    g.data_ = nullptr;
    g.nt_ = -2;
    g.size1_ = 0;
    g.size2_ = 0;
    g.element_size_ = 0;
}
/** \brief <b> Overloaded operator '=', which copies data from an existing `function` object(right-value reference)</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Overloaded operator '='. It copies data from an existing `function` object(right-value reference).
* > It copies all data from the existing `function` object. The method is specified for right value
* > reference when C++11 is used.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param g
* > the `function` object one wants to copy.
*/

template <typename T>
function<T> &function<T>::operator=(function &&g) noexcept {
    if (&g == this)
        return *this;

    data_ = g.data_;
    nt_ = g.nt_;
    size1_ = g.size1_;
    size2_ = g.size2_;
    element_size_ = g.element_size_;

    g.data_ = nullptr;
    g.nt_ = -2;
    g.size1_ = 0;
    g.size2_ = 0;
    g.element_size_ = 0;

    return *this;
}
#endif
/** \brief <b> Overloaded operator '='. It copies data from an existing `function` object</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Overloaded operator '='. It copies data from an existing `function` object.
* > It copies all data from the existing `function` object. The method is specified for right value
* > reference when C++11 is used.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param g
* > the `function` object one wants to copy.
*/
template <typename T>
function<T> &function<T>::operator=(const function &g) {
    int len;
    if (this == &g)
        return *this;
    if (data_ != 0)
        delete[] data_;
    nt_ = g.nt_;
    size1_ = g.size1_;
    size2_ = g.size2_;
    element_size_ = size1_ * size2_;
    len = (nt_ + 2) * size1_ * size2_;
    if (len > 0) {
    	total_size_ = len;
        data_ = new cplx[len];
        memcpy(data_, g.data_, sizeof(cplx) * len);
    } else {
        data_ = 0;
    }
    return *this;
}
/** \brief <b> Resize the time length and/or matrix dimension of a square-matrix function</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Resize the time length and/or matrix dimension of a square-matrix function
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param nt
* > the new number of time-steps.
* @param size1
* > the new leading dimension of the matrix
*/

template <typename T>
void function<T>::resize(int nt, int size1) {
    int len = (nt + 2) * size1 * size1;
    assert(nt >= -1 && size1 >= 0);
    if (data_ != 0)
        delete[] data_;
    if (len == 0)
        data_ = 0;
    else {
        data_ = new cplx[len];
    }
    size1_ = size1;
    size2_ = size1;
    element_size_ = size1_ * size1_;
    nt_ = nt;
    total_size_ = len;
}
/** \brief <b> Set all data to zero for the `function` class</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Set all data to zero for the `function` class.
*
* <!-- ARGUMENTS
*      ========= -->
*
*/

template <typename T>
void function<T>::set_zero(void) {
    int len = element_size_ * (nt_ + 2);
    if (len > 0)
        memset(data_, 0, sizeof(cplx) * len);
}
/** \brief <b> Set all data to a constant(scalar or matrix) for the `function` class</b>

* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Set all data to a constant(scalar or matrix) for the `function` class
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param f0
* > the pointer of an array storing the constant data in the form `a[i * size1_ + j] = matrix(i,j)`.
* > For a scalar, the size of array should be 1.
*/
template <typename T>
void function<T>::set_constant(std::complex<T> *f0) {
    int t;
    for (t = -1; t <= nt_; t++)
        element_set<T, LARGESIZE>(size1_, ptr(t), f0);
}
/** \brief <b> Set all data to a constant matrix for the `function` class</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Set all data to a constant matrix for the `function` class
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param M
* > The constant matrix as EigenMatrix
*/
template <typename T> template<class EigenMatrix>
void function<T>::set_constant(EigenMatrix &M){
  assert(M.rows() == size1_);
	int t;
	for(t=-1;t<=nt_;t++){
		set_value(t,M);
	}
}

//template <class EigenMatrix>
/** \brief <b> Set scalar value at a specific time point</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Set scalar value at a specific time point. The value is assumed to be a complex number.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > the time step to be set
* @param x
* > the new complex value of function at time step tstp.
*/
template <typename T> template<class EigenMatrix>
void function<T>::set_value(int tstp, cplx x) {
    int i1, i2;
    cplx *ft;
    assert(1 == size1_);
    assert(1 == size2_);
    *ptr(tstp) =x;
}
/** \brief <b> Set matrix value at a specific time point</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Set matrix value at a specific time point
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > the time step to be set
* @param M
* > the new matrix value of function at time step tstp as an EigenMatrix.
*/

template <typename T>
template <class EigenMatrix>
void function<T>::set_value(int tstp, EigenMatrix &M) {
    int i1, i2;
    cplx *ft;
    assert(M.rows() == size1_);
    assert(M.cols() == size2_);
  
    ft = ptr(tstp);
    for (i1 = 0; i1 < size1_; i1++)
        for (i2 = 0; i2 < size2_; i2++)
            ft[i1 * size2_ + i2] = M(i1, i2);
}
/** \brief <b> Get matrix value of this `function` object at a specific time point</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Get matrix value of this `function` object at a specific time point.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > the time step
* @param M
* > the Eigen matrix in which the value will be stored
*/

template <typename T>
template <class EigenMatrix>
void function<T>::get_value(int tstp, EigenMatrix &M) const {
	int i1, i2;
	const cplx *ft;
	M.resize(size1_, size2_);
	ft = ptr(tstp);
	for (i1 = 0; i1 < size1_; i1++)
		for (i2 = 0; i2 < size2_; i2++)
			M(i1, i2) = ft[i1 * size2_ + i2];
}
/** \brief <b> Multiply the function with a scalar</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Multiply the function with a scalar. Function values at all times would be multiplied
* > with the same scalar value.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param weight
* > the scalar to be multiplied with all function values
*/
template <typename T>
void function<T>::smul(T weight)
{
	for(int m = 0; m < total_size_; m++)
		data_[m] *= weight;
}
/** \brief <b> Multiply function values with a matrix-valued function from the left</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Multiply the function with another matrix-valued function from the left time-pointwisely.
* > The shape of the other matrix-valued function must be the same as the original function.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param f
* > the array storing the matrix-valued function.
* @param weight
* > the factor that function values are multiplied with in addition to the matrix.
*/
template <typename T>
void function<T>::left_multiply(cplx* f, T weight)
{
	cplx *mytemp, *xtemp, *ftemp;
	xtemp = new cplx[element_size_];
	for (int tstp = -1; tstp <= nt_; tstp++) {
		ftemp = f + (tstp + 1) * element_size_;
		mytemp = ptr(tstp);
		element_mult<T, LARGESIZE>(size1_, xtemp,
				ftemp, mytemp);
		element_smul<T, LARGESIZE>(size1_, xtemp, weight);
		element_set<T, LARGESIZE>(size1_, mytemp, xtemp);
	}
	delete[] xtemp;
}
/** \brief <b> Multiply the function with another matrix-valued function from the right</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Multiply the function with another matrix-valued function from the right time-pointwisely.
* > The shape of the other matrix-valued function is assumed to be the same as the original function.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param f
* > the array storing the matrix-valued function.
* @param weight
* > the factor that function values are multiplied with in addition to the matrix.
*/

template <typename T>
void function<T>::right_multiply(cplx* f, T weight)
{
	cplx *mytemp, *xtemp, *ftemp;
	xtemp = new cplx[element_size_];
	for (int tstp = -1; tstp <= nt_; tstp++) {
		ftemp = f + (tstp + 1) * element_size_;
		mytemp = ptr(tstp);
		element_mult<T, LARGESIZE>(size1_, xtemp,
				mytemp, ftemp);
		element_smul<T, LARGESIZE>(size1_, xtemp, weight);
		element_set<T, LARGESIZE>(size1_, mytemp, xtemp);
	}
	delete[] xtemp;
}
/** \brief <b> Multiply function values with a matrix-valued function from the left</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Multiply the function with another matrix-valued function from the left time-pointwisely.
* > The shape of the other matrix-valued function must be the same as the original function.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param f
* > the function object to be multiplied with this object.
* @param weight
* > the extra scalar factor that function values are multiplied with.
*/
template <typename T>
void function<T>::left_multiply(function<T> &f, T weight) {
  assert(this -> nt_ == f.nt_);
  assert(this -> size1_ == f.size1_);
  assert(this -> size2_ == f.size2_);
  this -> left_multiply(f.ptr(-1), weight);
}
/** \brief <b> Multiply function values with a matrix-valued function from the right</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Multiply function values with a matrix-valued function from the right. Function values
* > would be multiplied with the other function at the same time-step.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param f
* > the function object to be multiplied with this object.
* @param weight
* > the extra scalar factor that function values are multiplied with.
*/
template <typename T>
void function<T>::right_multiply(function<T> &f, T weight) {
  assert(this -> nt_ == f.nt_);
  assert(this -> size1_ == f.size1_);
  assert(this -> size2_ == f.size2_);
  this->right_multiply(f.ptr(-1), weight);
}
/** \brief <b> Set a matrix element from an Eigen matrix</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Set a matrix element of the `function` object at all time steps from an Eigen matrix.
* > The matrix element of the function at each time step would be set to the same matrix
* > element of the Eigen matrix.
*
* <!-- ARGUMENTS
*      ========= -->
* @param tstp
* > time step of the matrix element to be set
* @param i1
* > the first index of the matrix element to be set
* @param i2
* > the second index of the matrix element to be set
* @param M
* > the matrix used to set the value
* @param j1
* > the first index of the matrix element used to set
* @param j2
* > the second index of the matrix element used to set
*/

template <typename T> template<class EigenMatrix>
void function<T>::set_matrixelement(int tstp,int i1,int i2,EigenMatrix &M,int j1,int j2){
	assert(0<=i1 && i1<size1_ && 0<=i2 && i2<size2_);
	assert(0<=j1 && j1<M.rows() && 0<=j2 && j2<M.cols());
	cplx *ft;
	ft=ptr(tstp);
	ft[i1*size2_+i2]=M(j1,j2);
}
/** \brief <b> Get a matrix element of the `function` object</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Get a matrix element of the `function` object at all time steps from an Eigen matrix.
* > The matrix element of the function at each time step will be stored in another `function`
* > object.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i1
* > the first index of the matrix element to be set
* @param i2
* > the second index of the matrix element to be set
* @param g
* > the `function` object that the matrix element will be stored. If it is matrix-valued, only one
* > matrix element will be changed.
*/
template <typename T>
void function<T>::get_matrixelement(int i1, int i2, function<T> &g)
{
  assert(this -> nt_ == g.nt_);
  assert(i1 <= this -> size1_);
  assert(i2 <= this -> size2_);
		cplx *target_position, *my_data;
	for(int tstp = -1; tstp <= nt_; tstp++)
	{
		my_data = this -> ptr(tstp) + i1 * size2_ + i2;
		target_position = g.ptr(tstp);
		element_set <T, 1> (1, target_position, my_data);
		//*target_position = *my_data;
	}
}
/** \brief <b> Increase the `function` object by another function.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Increase the `function` object by another function which is firstly multiplied with a weight..
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param g
* > the `function` object whose data will be added to the original `function` object.
* @param weight
* > the weight to be multiplied to the g `function`.
*/
template <typename T>
void function<T>::incr(function<T> &g, T weight)
{
  assert(this -> nt_ == g.nt_);
		for(int m = 0; m < total_size_; m++)
		{
			this -> data_[m] += g.data_[m] * weight;
		}
}
/** \brief <b> Set a matrix element from a `function` object</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Set a matrix element from a `function` object at all time steps. The `function` should have
* > the same size as the original `function`.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i1
* > the first index of the matrix element to be set
* @param i2
* > the second index of the matrix element to be set
* @param g
* > the function object used to set the elements.
** @param j1
* > the first index of the matrix element used to set
* @param j2
* > the second index of the matrix element used to set
*/
template <typename T>
void function<T>::set_matrixelement(int i1,int i2,function &g,int j1,int j2){
    assert(0<=i1 && i1<size1_ && 0<=i2 && i2<size2_);
    assert(0<=j1 && j1<g.size1() && 0<=j2 && j2<g.size2());

    for(int tstp=-1;tstp<=nt_;tstp++){
        cdmatrix M;
        g.get_value(tstp,M);
        set_matrixelement(tstp,i1,i2,M,j1,j2);
    }

}
/** \brief <b> Store the `function` data to a txt file</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Store the `function` data to a file. It will be stored in txt form.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param file
* > the file name
* @param precision
* > the precision of stored `function` data
*/
template <typename T>
void function<T>::print_to_file(const char *file, int precision) const {
    int i, l, sg = element_size_;
    std::ofstream out;
    out.open(file, std::ios::out);
    out.precision(precision);
    out << "# " << nt_ << " " << size1_ << std::endl;
    for (i = -1; i <= nt_; i++) {
        for (l = 0; l < sg; l++)
            out << ptr(i)[l].real() << " " << ptr(i)[l].imag() << " ";
        out << std::endl;
    }
    out.close();
}
/** \brief <b> Read the `function` data from a txt file</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Read the `function` from a txt file. The data file is probably written by print_to_file method,
* > or is assumed to have the same format.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param file
* > the file name
*/

template <typename T>
void function<T>::read_from_file(const char *file) {
    int i, n, l, size1, sg;
    double real, imag;
    std::string s;
    std::ifstream out;
    out.open(file, std::ios::in);
    if (!(out >> s >> n >> size1)) {
        std::cerr << "read function from file " << file << " error in file"
                  << std::endl;
        abort();
    }
    if (n > nt_ || size1 != size1_)
        resize(n, size1);
    sg = element_size_;
    for (i = -1; i <= n; i++) {
        for (l = 0; l < sg; l++) {
            if (!(out >> real >> imag)) {
                std::cerr << "read function from file " << file
                          << " error at  t= " << i << std::endl;
                abort();
            }
            ptr(i)[l] = std::complex<T>(real, imag);
        }
    }
    out.close();
}
/** \brief <b> Read the `function` data from a txt file</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Read the `function` data from a txt file. The file is probably written by print_to_file method,
* > or is assumed to have the same format.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param nt1
* > the number of time steps to be read
* @param file
* > the file name
*/
template <typename T>
void function<T>::read_from_file(int nt1, const char *file) {
    int i, n, l, size1, sg;
    double real, imag;
    std::string s;
    std::ifstream out;
    out.open(file, std::ios::in);
    if (!(out >> s >> n >> size1)) {
        std::cerr << "read function from file " << file << " error in file"
                  << std::endl;
        abort();
    }

    assert(nt1 >= -1 && nt1 <= nt_);
    assert(nt1 <= n);
    assert(size1 == size1_);
    
    sg = size1 * size1;
    for (i = -1; i <= n; i++) {
        for (l = 0; l < sg; l++) {
            if (!(out >> real >> imag)) {
                std::cerr << "read function from file " << file
                          << " error at  t= " << i << std::endl;
                abort();
            }
            if (i <= nt1)
                ptr(i)[l] = std::complex<T>(real, imag);
        }
    }
    out.close();
}

#if CNTR_USE_HDF5 == 1
// identical from herm_matrix<T>
/** \brief <b> Write the `function` data to a hdf5 file</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Write the `function` data to a hdf5 file. Four parameters: nt, size1, size2 and element_size are stored.
* > data array is stored as complex numbers.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param group_id
* > the group_id of the hdf5 data
*/
template <typename T>
void function<T>::write_to_hdf5(hid_t group_id) const {
    store_int_attribute_to_hid(group_id, std::string("nt"), nt_);
    store_int_attribute_to_hid(group_id, std::string("size1"), size1_);
    store_int_attribute_to_hid(group_id, std::string("size2"), size2_);
    store_int_attribute_to_hid(group_id, std::string("element_size"),
                               element_size_);
    hsize_t len_shape = 3, shape[3];
    shape[1] = size1_;
    shape[2] = size2_;
    if (nt_ > -2) {
        shape[0] = nt_ + 2;
        store_cplx_array_to_hid(group_id, std::string("data"), this->data_,
                                shape, len_shape);
    }
}
// identical from herm_matrix<T>
/** \brief <b> Write the `function` data to a hdf5 file</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Write the `function` data to a hdf5 file. Four parameters: nt, size1, size2 and element_size are stored.
* > data array is stored as complex numbers. The data is stored under a specified group name.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param group_id
* > the group_id of the hdf5 data
* @param groupname
* > the specified (sub)group name of stored function
*/

template <typename T>
void function<T>::write_to_hdf5(hid_t group_id, const char *groupname) const {
    hid_t sub_group_id = create_group(group_id, groupname);
    this->write_to_hdf5(sub_group_id);
    close_group(sub_group_id);
}
// identical from herm_matrix<T>
/** \brief <b> Write the `function` data to a hdf5 file</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Write the `function` data to a hdf5 file. Four parameters: nt, size1, size2 and element_size are stored.
* > data array is stored as complex numbers. The data is stored to a file with given name and under a specified group name.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param filename
* > the name of the file in which the function would be stored
* @param groupname
* > the specified (sub)group name of stored function
*/
template <typename T>
void function<T>::write_to_hdf5(const char *filename,
                                const char *groupname) const {
    hid_t file_id = open_hdf5_file(filename);
    this->write_to_hdf5(file_id, groupname);
    close_hdf5_file(file_id);
}
/** \brief <b> Read the `function` data from a hdf5 file</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Read the `function` data from a hdf5 file. The file is probably written by write_to_hdf5 method,
* > or is assumed to have the same format.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param group_id
* > the group_id of the hdf5 data
*/

template <typename T>
void function<T>::read_from_hdf5(hid_t group_id) {
    // -- Read dimensions
    int nt = read_primitive_type<int>(group_id, "nt");
    int size1 = read_primitive_type<int>(group_id, "size1");
    this->resize(nt, size1);
    if (nt > -2) {
        hsize_t data_size = (nt + 2) * element_size_;
        read_primitive_type_array(group_id, "data", data_size, this->data_);
    }
}
/** \brief <b> Read the `function` data from a hdf5 file</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Read the `function` data from a hdf5 file. The file is probably written by write_to_hdf5 method,
* > or is assumed to have the same format.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param group_id
* > the group_id of the hdf5 data
* @param groupname
* > the specified (sub)group name of stored function
*/

template <typename T>
void function<T>::read_from_hdf5(hid_t group_id, const char *groupname) {
    hid_t sub_group_id = open_group(group_id, groupname);
    this->read_from_hdf5(sub_group_id);
    close_group(sub_group_id);
}
/** \brief <b> Read the `function` data from a hdf5 file</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Read the `function` data from a hdf5 file. The file is probably written by write_to_hdf5 method,
* > or is assumed to have the same format.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param filename
* > the name of the file in which the function is stored
* @param groupname
* > the specified (sub)group name of stored function
*/

template <typename T>
void function<T>::read_from_hdf5(const char *filename,
                                 const char *groupname) {
    hid_t file_id = read_hdf5_file(filename);
    this->read_from_hdf5(file_id, groupname);
    close_hdf5_file(file_id);
}
/** \brief <b> Read the `function` data from a hdf5 file</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Read the `function` data from a hdf5 file. The file is probably written by write_to_hdf5 method,
* > or is assumed to have the same format.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param nt1
* > the number of time steps to be read
* @param group_id
* > the group id of the hdf5 data
*/
template <typename T>
void function<T>::read_from_hdf5(int nt1, hid_t group_id) {
    function<T> ftmp;
    ftmp.read_from_hdf5(group_id);

    assert(nt1 >= -1 && nt1 <= ftmp.nt());
    assert(nt1 >= -1 && nt1 <= nt_);
    assert(ftmp.size1() == size1_);
    assert(ftmp.element_size() == element_size_);
  
    memcpy(data_, ftmp.data_, sizeof(cplx) * (nt1 + 2) * element_size_);
}
/** \brief <b> Read the `function` data from a hdf5 file</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Read the `function` data from a hdf5 file. The file is probably written by write_to_hdf5 method,
* > or is assumed to have the same format.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param nt1
* > the number of time steps to be read
* @param group_id
* > the group id of the hdf5 data
* @param groupname
* > the specified group in which the function data is stored
*/
template <typename T>
void function<T>::read_from_hdf5(int nt1, hid_t group_id,
                                 const char *groupname) {
    hid_t sub_group_id = open_group(group_id, groupname);
    this->read_from_hdf5(nt1, sub_group_id);
    close_group(sub_group_id);
}
/** \brief <b> Read the `function` data from a hdf5 file</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Read the `function` data from a hdf5 file. The file is probably written by write_to_hdf5 method,
* > or is assumed to have the same format.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param nt1
* > the number of time steps to be read
* @param filename
* > the name of the file from which the function data is read
* @param groupname
* > the specified group in which the function data is stored
*/
template <typename T>
void function<T>::read_from_hdf5(int nt1, const char *filename,
                                 const char *groupname) {
    hid_t file_id = read_hdf5_file(filename);
    this->read_from_hdf5(nt1, file_id, groupname);
    close_hdf5_file(file_id);
}
#endif

/* #######################################################################################
#
#   MPI UTILS
#
########################################################################################*/
#if CNTR_USE_MPI==1
/** \brief <b> Broadcast a time step of the `function` using MPI</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Broadcast a time step of the `function` using MPI. The value at one time step
* > is broadcasted from the `root` process to all other processes.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > the time step. The value at time point `tstp` is broadcasted
* @param root
* > the `root` index. The process indexed by `root` will send its data. Others will receive data.
*/
template <typename T>
void function<T>::Bcast_timestep(int tstp,int root){
  int numtasks,taskid;
  cdmatrix ftemp;
  ftemp.resize(size1_,size2_);
  numtasks=MPI::COMM_WORLD.Get_size();
  taskid=MPI::COMM_WORLD.Get_rank();
    if(taskid==root) this->get_value(tstp,ftemp);
    if(sizeof(T)==sizeof(double)) MPI::COMM_WORLD.Bcast(ftemp.data(),ftemp.size(),MPI::DOUBLE_COMPLEX,root);
    else MPI::COMM_WORLD.Bcast(ftemp.data(),ftemp.size(),MPI::COMPLEX,root);
    if(taskid!=root) this->set_value(tstp,ftemp);
}
#endif


} // namespace cntr

#endif  // CNTR_FUNCTION_IMPL_H
