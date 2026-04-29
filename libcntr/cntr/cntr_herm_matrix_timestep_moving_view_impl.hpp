#ifndef CNTR_HERM_TIMESTEP_MOVING_VIEW_IMPL_H
#define CNTR_HERM_TIMESTEP_MOVING_VIEW_IMPL_H

#include "cntr_herm_matrix_timestep_moving_view_decl.hpp"
#include "cntr_elements.hpp"
//#include "cntr_exception.hpp"

namespace cntr {

  
#define CPLX std::complex<T>
  template <typename T>
  herm_matrix_timestep_moving_view<T>::herm_matrix_timestep_moving_view() {
    ret_=0;
    les_=0;
    tc_=-1;
    size1_=0;
    size2_=0;
    element_size_=0;
    sig_= -1;
  }
  
  template <typename T>
  herm_matrix_timestep_moving_view<T>::~herm_matrix_timestep_moving_view() {
    // donothing
  }
  
/**
 * @brief Copy constructor for herm_matrix_timestep_moving_view.
 *
 * This constructor creates a new instance of `herm_matrix_timestep_moving_view`
 * by copying the internal data pointers and attributes from another instance.
 * Note that this is a shallow copy, meaning that the new object shares the same 
 * data pointers (`ret_` and `les_`) as the original, rather than duplicating the data.
 *
 * @tparam T Must be either `double` or `float`.
 * @param g The `herm_matrix_timestep_moving_view` object to copy from.
 */
  template <typename T>
  herm_matrix_timestep_moving_view<T>::herm_matrix_timestep_moving_view(const herm_matrix_timestep_moving_view &g) {
    tc_=g.tc_;
    size1_=g.size1_;
    size2_=g.size2_;
    element_size_=g.element_size_;
    sig_= g.sig_;
    ret_=g.ret_;
    les_=g.les_;
  }

/**
 * @brief Assignment operator for herm_matrix_timestep_moving_view.
 *
 * This operator assigns the values of an existing `herm_matrix_timestep_moving_view`
 * object to the current instance. It performs a **shallow copy**, meaning that 
 * the internal data pointers (`ret_` and `les_`) are shared rather than duplicated.
 *
 * @tparam T Must be either `double` or `float`.
 * @param g The `herm_matrix_timestep_moving_view` object to copy from.
 * @return A reference to the assigned `herm_matrix_timestep_moving_view` object.
 */  
  template <typename T>
  herm_matrix_timestep_moving_view<T> &herm_matrix_timestep_moving_view<T>::operator=(const herm_matrix_timestep_moving_view &g){
    if(this == &g) return *this;
    tc_=g.tc_;
    size1_=g.size1_;
    size2_=g.size2_;
    element_size_=size1_*size2_;
    sig_= g.sig_;
    ret_=g.ret_;
    les_=g.les_;
    return *this;
  }

  
/**
 * @brief Constructs a `herm_matrix_timestep_moving_view` from a `herm_matrix_timestep_moving` object.
 *
 * This constructor creates a **view** of an existing `herm_matrix_timestep_moving` object 
 * without copying its data. The internal pointers (`ret_` and `les_`) are directly assigned 
 * to the corresponding data in `g`, allowing for efficient access while avoiding duplication.
 *
 * @tparam T Must be either `double` or `float`.
 * @param g The `herm_matrix_timestep_moving` object from which to create the view.
 *
 * @note Since this is a shallow copy, modifications to `g` will affect this view and vice versa.
 * @warning The lifetime of `g` must exceed that of this view to prevent dangling pointers.
 */
  template<typename T>
  herm_matrix_timestep_moving_view<T>::herm_matrix_timestep_moving_view(herm_matrix_timestep_moving<T> &g){
    tc_=g.tc();
    size1_=g.size1();
    size2_=g.size2();
    element_size_=size1_*size2_;
    sig_= g.sig();
    ret_=g.retptr(0);
    les_=g.lesptr(0);
  }

  
/**
 * @brief Constructs a `herm_matrix_timestep_moving_view` from a `herm_matrix_moving` object 
 *        at a specific relative time slice.
 *
 * This constructor creates a **view** of the Green’s function stored at a given time slice `deltG` 
 * of a `herm_matrix_moving`object . The internal pointers (`ret_` and `les_`)
 * reference the corresponding data in `g`.
 *
 * @tparam T Must be either `double` or `float`.
 * @param deltG The relative time step at which to create the view. Must be in the range `[0, g.tc()]`.
 * @param g The `herm_matrix_moving` object storing the moving window of the Green's function.
 *
 * @note This is a **shallow copy**; no new memory is allocated. Modifications to `g` affect this view.
 * @warning Ensure that `g` remains valid while this view is in use to prevent dangling pointers.
 */
  template <typename T>
  herm_matrix_timestep_moving_view<T>::herm_matrix_timestep_moving_view(herm_matrix_moving<T> &g,int deltG){
    assert(deltG>=0 && deltG<=g.tc());
    tc_=g.tc();
    size1_=g.size1();
    size2_=g.size2();
    element_size_=size1_*size2_;
    sig_= g.sig();
    ret_=g.retptr(deltG,0);
    les_=g.lesptr(deltG,0);
  }

    template <typename T>
  herm_matrix_timestep_moving_view<T>::herm_matrix_timestep_moving_view(herm_matrix_moving<T> &g){
    tc_=g.tc();
    size1_=g.size1();
    size2_=g.size2();
    element_size_=size1_*size2_;
    sig_= g.sig();
    ret_=g.retptr(0,0);
    les_=g.lesptr(0,0);
  }

/**
 * @brief Sets a specific matrix element from another `herm_matrix_timestep_moving_view` object.
 *
 * This function copies the matrix element `(j1, j2)` in the source object `g`
 * to the position `(i1, i2)` in the current object  for all times of the timeslice
 * `this^[X](t_lead,t_lead-j)_{i1,i2}` set to `g^[X](t_lead,t_lead-j)_{j1,j2}` for `j=0...tc` and `X=ret,les`
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param i1 Row index of the destination matrix element.
 * @param i2 Column index of the destination matrix element.
 * @param g Source `herm_matrix_timestep_moving_view` object.
 * @param j1 Row index of the source matrix element.
 * @param j2 Column index of the source matrix element.
 *
 * @note The time cutoffs (`tc_`) of both objects must match.
 */
  template <typename T>
  void herm_matrix_timestep_moving_view<T>::set_matrixelement(int i1, int i2, herm_matrix_timestep_moving_view<T> &g,int j1, int j2){
    int i, sij = i1 * size2_ + i2, tij = j1 * g.size2() + j2;
    assert(tc_==g.tc_);
    assert(i1>=0 && i1 <=size1_ -1);
    assert(i2>=0 && i2 <=size2_ -1);
    assert(j1>=0 && j1 <=g.size1_-1);
    assert(j2>=0 && j2 <=g.size2_-1);
    for(i =0; i<=tc_;i++)
      retptr(i)[sij]=g.retptr(i)[tij];
    for(i =0; i<=tc_;i++)
      lesptr(i)[sij]=g.lesptr(i)[tij];
  }

/**
 * @brief Sets a specific matrix element from another `herm_matrix_timestep_moving` object.
 *
 * This function copies the matrix element `(j1, j2)` in the source object `g`
 * to the position `(i1, i2)` in the current object  for all times of the timeslice
 * `this^[X](t_lead,t_lead-j)_{i1,i2}` set to `g^[X](t_lead,t_lead-j)_{j1,j2}` for `j=0...tc` and `X=ret,les`
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param i1 Row index of the destination matrix element.
 * @param i2 Column index of the destination matrix element.
 * @param g Source `herm_matrix_timestep_moving` object.
 * @param j1 Row index of the source matrix element.
 * @param j2 Column index of the source matrix element.
 *
 * @note The time cutoffs (`tc_`) of both objects must match.
 */
  template <typename T>
  void herm_matrix_timestep_moving_view<T>::set_matrixelement(int i1, int i2, herm_matrix_timestep_moving<T> &g,int j1, int j2){
    herm_matrix_timestep_moving_view<T> tmp(g);
    this->set_matrixelement(i1,i2,tmp,j1,j2);
  }

/**
 * @brief Sets a specific matrix element from a given  timeslice deltG of another `herm_matrix_moving` object 
 *
 * This function copies the matrix element `(j1, j2)` in the source object `g`
 * to the position `(i1, i2)` in the current object  for all times of the timeslice
 * `this^[X](t_lead,t_lead-j)_{i1,i2}` set to `g^[X](t_lead-deltG,t_lead-deltG-j)_{j1,j2}` for `j=0...tc` and `X=ret,les`
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param i1 Row index of the destination matrix element.
 * @param i2 Column index of the destination matrix element.
 * @param g Source `herm_matrix_moving` object.
 * @param deltG The relative time step at which to create the view. Must be in the range `[0, g.tc()]`.
 * @param j1 Row index of the source matrix element.
 * @param j2 Column index of the source matrix element.
 *
 * @note The time cutoffs (`tc_`) of both objects must match.
 */
 
  template <typename T>
  void herm_matrix_timestep_moving_view<T>::set_matrixelement(int i1, int i2, herm_matrix_moving<T> &g,int deltG,int j1, int j2){
    herm_matrix_timestep_moving_view<T> tmp(g,deltG);
    this->set_matrixelement(i1,i2,tmp,j1,j2);
  }


/**
 * @brief Copies the data from another `herm_matrix_timestep_moving_view` object.
 *
 * This function copies the retarded and lesser Green’s function components 
 * from the given `herm_matrix_timestep_moving_view` object into the current instance.
 * `this^[X](t_lead,t_lead-j)` set to `g^[X](t_lead,t_lead-j)` for `j=0...tc` and `X=ret,les`
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param g1 Source `herm_matrix_timestep_moving_view` object.
 *
 * @note The time cut-off (`tc_`) and matrix sizes (`size1_`) of both objects must be the same.
 */
  template <typename T>
void herm_matrix_timestep_moving_view<T>::set_timestep(herm_matrix_timestep_moving_view<T> &g1) {
  assert(tc_== g1.tc() && "tc_== g1.tc()");
  assert(g1.size1() == size1_ && "g1.size1() == size1_");
  memcpy(retptr(0), g1.retptr(0),sizeof(cplx) * (tc_ + 1) * element_size_);
  memcpy(lesptr(0), g1.lesptr(0),sizeof(cplx) * (tc_ + 1) * element_size_);
}

/**
 * @brief Copies the data from another `herm_matrix_timestep_moving` object.
 *
 * This function copies the retarded and lesser Green’s function components 
 * from the given `herm_matrix_timestep_moving` object into the current instance.
 * `this^[X](t_lead,t_lead-j)` set to `g^[X](t_lead,t_lead-j)` for `j=0...tc` and `X=ret,les`
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param g1 Source `herm_matrix_timestep_moving_view` object.
 *
 * @note The time cut-off (`tc_`) and matrix sizes (`size1_`) of both objects must be the same.
 */
template <typename T>
void herm_matrix_timestep_moving_view<T>::set_timestep(herm_matrix_timestep_moving<T> &g) {
  assert(tc_== g.tc() && "tc_== g.tc()");
  assert(g.size1() == size1_ && "g.size1() == size1_");
  memcpy(retptr(0), g.retptr(0),sizeof(cplx) * (tc_ + 1) * element_size_);
  memcpy(lesptr(0), g.lesptr(0), sizeof(cplx) * (tc_ + 1) * element_size_);
}

/**
 * @brief Copies the data from the timeslice deltG of another `herm_matrix_moving_view` object
 *
 * This function copies the retarded and lesser Green’s function components 
 * from the given `herm_matrix_moving` object into the current instance.
 * `this^[X](t_lead,t_lead-j)` set to `g^[X](t_lead-deltG,t_lead-deltG-j)` for `j=0...tc` and `X=ret,les`
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param g1 Source `herm_matrix_timestep_moving_view` object.
 * @param deltG The relative time step at which to create the view. Must be in the range `[0, g.tc()]`.
 *
 * @note The time cut-off (`tc_`) and matrix sizes (`size1_`) of both objects must be the same.
 */
template <typename T>
void herm_matrix_timestep_moving_view<T>::set_timestep(herm_matrix_moving<T> &g,int deltG) {
  assert(tc_== g.tc() && "tc_== g.tc()");
  assert(deltG >= 0 && deltG <= g.tc() && "i >= 0 && i <= g.tc()");
  assert(g.size1() == size1_ && "g.size1() == size1_");
  memcpy(retptr(0), g.retptr(deltG, 0), sizeof(cplx) * (tc_ + 1) * element_size_);
  memcpy(lesptr(0), g.lesptr(deltG, 0),sizeof(cplx) * (tc_ + 1) * element_size_);	
}
 

/**
 * @brief Clears the timestep data by setting all elements to zero.
 *
 * This function sets all elements of the retarded and lesser Green’s function 
 * components to zero for the entire range of the timestep:
 * `this^[X](t_lead,t_lead-j)` set to 0 for `j=0...tc` and `X=ret,les`
 *
 * @tparam T Floating-point type (`double` or `float`).
 *
 */
  template <typename T>
void herm_matrix_timestep_moving_view<T>::clear_timestep(void) {
        memset(retptr(0), 0, sizeof(cplx) * (tc_ + 1) * element_size_);
        memset(lesptr(0), 0, sizeof(cplx) * (tc_ + 1) * element_size_);
}


/**
 * @brief Clears the timestep data by setting all elements to zero.
 *
 * This function sets all elements of the retarded and lesser Green’s function 
 * components to zero for the entire range of the timestep
 * `this^[X](t_lead,t_lead-j)` set to 0 for `j=0...tc` and `X=ret,les`
 *
 * @tparam T Floating-point type (`double` or `float`).
 *
 */
  template <typename T>
void herm_matrix_timestep_moving_view<T>::set_timestep_zero(void) {
        memset(retptr(0), 0, sizeof(cplx) * (tc_ + 1) * element_size_);
        memset(lesptr(0), 0, sizeof(cplx) * (tc_ + 1) * element_size_);
}

  
/* #######################################################################################
#
#   WRITING ELEMENTS FROM ANY MATRIX TYPE
#   OR FROM COMPLEX NUMBERS (then only the (0,0) element is addressed for dim>0)
#
########################################################################################*/
#define herm_matrix_SET_ELEMENT_MATRIX                                      \
      {                                                                     \
         int r, s;                                                          \
         for (r = 0; r < size1_; r++)                                       \
            for (s = 0; s < size2_; s++)                                    \
               x[r * size2_ + s] = M(r, s);                                 \
         }


/**
 * @brief Sets the retarded Green's function component for a given relative time step.
 *
 * This function assigns values from the provided matrix `M` to the retarded component 
 * at a specific time difference `deltj`. The matrix elements are copied into the 
 * corresponding positions of the internal storage: retarded component at `(t_lead,t_lead-deltj)` set to `M`
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @tparam Matrix Type of the input matrix, which must support element access via `operator()`.
 * @param deltj Relative time index (must be within the time cutoff `tc_`).
 * @param M Input matrix containing the new values.
 *
 * `deltj` must be within the allowed range (`0 <= deltj <= tc_`).
 * The dimensions of `M` must match `size1_ × size2_`.
 */  
template<typename T> template <class Matrix> void herm_matrix_timestep_moving_view<T>::set_ret(int deltj,Matrix &M){
         assert(deltj <= tc_);
         cplx *x=retptr(deltj);
         herm_matrix_SET_ELEMENT_MATRIX
      }

/**
 * @brief Sets the retarded Green's function component for a given relative time step.
 *
 * This function assigns the scalar value `x` to the retarded component at a specific 
 * time difference `deltj`. Specifically, it sets the retarded component at `(t_lead, t_lead - deltj)`
 * to `x`. This function is intended for cases where the matrix size is 1×1.
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param deltj Relative time index (must be within the time cutoff `tc_`).
 * @param x Scalar value to assign.
 *
 * `deltj` must be within the allowed range (`0 <= deltj <= tc_`).
 */  
template<typename T> void herm_matrix_timestep_moving_view<T>::set_ret(int deltj, cplx x){
      assert(deltj <= tc_);
      retptr(deltj)[0]=x;
   }
  

/**
 * @brief Sets the lesser Green's function component for a given relative time step.
 *
 * This function assigns values from the provided matrix `M` to the lesser component 
 * at a specific time difference `deltj`. The matrix elements are copied into the 
 * corresponding positions of the internal storage: lesser component at `(t_lead,t_lead-deltj)` set to `M`
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @tparam Matrix Type of the input matrix, which must support element access via `operator()`.
 * @param deltj Relative time index (must be within the time cutoff `tc_`).
 * @param M Input matrix containing the new values.
 *
 * `deltj` must be within the allowed range (`0 <= deltj <= tc_`).
 * The dimensions of `M` must match `size1_ × size2_`.
 */  
  template<typename T> template <class Matrix> void herm_matrix_timestep_moving_view<T>::set_les(int deltj, Matrix &M){
      assert(deltj <= tc_);
      cplx *x=lesptr(deltj);
      herm_matrix_SET_ELEMENT_MATRIX
   }
/**
 * @brief Sets the lesser Green's function component for a given relative time step.
 *
 * This function assigns the scalar value `x` to the lesser component at a specific 
 * time difference `deltj`. Specifically, it sets the lesser component at `(t_lead, t_lead - deltj) `
 * to `x`. This function is intended for cases where the matrix size is 1×1.
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param deltj Relative time index (must be within the time cutoff `tc_`).
 * @param x Scalar value to assign.
 *
 * `deltj` must be within the allowed range (`0 <= deltj <= tc_`).
 */  
  template<typename T> void herm_matrix_timestep_moving_view<T>::set_les(int deltj, cplx x){
      assert(deltj <= tc_);
      lesptr(deltj)[0]=x;
   }
  
////////////////////////////////////////////////////////////////////////////////////////
// the following routines are not very "efficent" but sometimes simple to
// implement
#define herm_matrix_READ_ELEMENT                                             \
{                                                                        \
   int r, s, dim = size1_;                                              \
   M.resize(dim, dim);                                                  \
   for (r = 0; r < dim; r++)                                            \
      for (s = 0; s < dim; s++)                                        \
         M(r, s) = x[r * dim + s];                                    \
   }


/**
 * @brief Retrieves the lesser Green's function component for a given relative time step.
 *
 * This function extracts the lesser component at a specific time difference `deltj` 
 * and stores it in the provided matrix `M`. The function resizes `M` to match 
 * the matrix dimensions (`size1_ × size2_`) before copying the values:
 * M set to lesser component at times `(t_lead,t_lead-deltj)`
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @tparam Matrix Type of the output matrix, which must support `resize()` and element access via `operator()`.
 * @param deltj Relative time index (must be within the time cutoff `tc_`).
 * @param M Output matrix where the extracted values will be stored.
 *
 * `deltj` must be within the allowed range (`0 <= deltj <= tc_`).
 */
  
  template <typename T> template <class Matrix>
      void herm_matrix_timestep_moving_view<T>::get_les(int deltj, Matrix &M){
         assert(deltj <= tc_);
         cplx *x;
         x = lesptr(deltj);
         herm_matrix_READ_ELEMENT
      }
/**
 * @brief Retrieves the lesser Green's function component for a given relative time step.
 *
 * This function extracts the lesser component at a specific time difference `deltj` 
 * and stores it in the provided variable `x`. The value of `x` corresponds to the 
 * lesser component at the time `(t_lead,t_lead-deltj)`.
 * This function is intended for cases where the matrix size is 1×1.
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param deltj Relative time index (must be within the time cutoff `tc_`).
 * @param x Output variable where the extracted value will be stored.
 *
 * `deltj` must be within the allowed range (`0 <= deltj <= tc_`).
 */
  template<typename T> void herm_matrix_timestep_moving_view<T>::get_les(int deltj, cplx &x){
      assert(deltj <= tc_);
      x=lesptr(deltj)[0];
   }
  
/**
 * @brief Retrieves the retarded Green's function component for a given relative time step.
 *
 * This function extracts the retarded component at a specific time difference `deltj` 
 * and stores it in the provided matrix `M`. The function resizes `M` to match 
 * the matrix dimensions (`size1_ × size2_`) before copying the values:
 * M set to retarded component at times `(t_lead,t_lead-deltj)`
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @tparam Matrix Type of the output matrix, which must support `resize()` and element access via `operator()`.
 * @param deltj Relative time index (must be within the time cutoff `tc_`).
 * @param M Output matrix where the extracted values will be stored.
 *
 * `deltj` must be within the allowed range (`0 <= deltj <= tc_`).
 */
  
template <typename T> template <class Matrix>
      void herm_matrix_timestep_moving_view<T>::get_ret(int deltj, Matrix &M) {
         assert(deltj <= tc_);
         cplx *x;
         x = retptr(deltj);
         herm_matrix_READ_ELEMENT
      }
/**
 * @brief Retrieves the retarded Green's function component for a given relative time step.
 *
 * This function extracts the retarded component at a specific time difference `deltj` 
 * and stores it in the provided variable `x`. The value of `x` corresponds to the 
 * retarded component at the time `(t_lead,t_lead-deltj)`.
 * This function is intended for cases where the matrix size is 1×1.
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param deltj Relative time index (must be within the time cutoff `tc_`).
 * @param x Output variable where the extracted value will be stored.
 *
 * `deltj` must be within the allowed range (`0 <= deltj <= tc_`).
 */  
  template<typename T> void herm_matrix_timestep_moving_view<T>::get_ret(int deltj, cplx &x){
      assert(deltj <= tc_);
      x=retptr(deltj)[0];
   }
  

  /**
 * @brief Retrieves the greater Green's function component for a given relative time step.
 *
 * This function extracts the greater component at a specific time difference `deltj` 
 * and stores it in the provided matrix `M`. The function resizes `M` to match 
 * the matrix dimensions (`size1_ × size2_`) before copying the values:
 * M set to greater component at times `(t_lead,t_lead-deltj)`
 * We use the lelation  `greater = retarded + lesser`.
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @tparam Matrix Type of the output matrix, which must support `resize()` and element access via `operator()`.
 * @param deltj Relative time index (must be within the time cutoff `tc_`).
 * @param M Output matrix where the extracted values will be stored.
 *
 * `deltj` must be within the allowed range (`0 <= deltj <= tc_`).
 */
  
template <typename T> template <class Matrix>
void herm_matrix_timestep_moving_view<T>::get_gtr(int deltj, Matrix &M) {
  Matrix M1;
  get_ret(deltj,M);
  get_les(deltj,M1);
  M+=M1;
}

  /**
 * @brief Retrieves the greater Green's function component for a given relative time step.
 *
 * This function extracts the greater component at a specific time difference `deltj` 
 * and stores it in the provided variable `x`. The value of `x` corresponds to the 
 * greater component at the time `(t_lead,t_lead-deltj)`.
 * We use the lelation  `greater = retarded + lesser`.
 * This function is intended for cases where the matrix size is 1×1.
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param deltj Relative time index (must be within the time cutoff `tc_`).
 * @param x Output variable where the extracted value will be stored.
 *
 * `deltj` must be within the allowed range (`0 <= deltj <= tc_`).
 */  
    template<typename T> void herm_matrix_timestep_moving_view<T>::get_gtr(int deltj, cplx &x){
      assert(deltj <= tc_);
      x=lesptr(deltj)[0]+retptr(deltj)[0];
   }


/**
 * @brief Computes the density matrix from the lesser Green's function component.
 *
 * This function retrieves the lesser Green's function component at the leading time step 
 * `t_lead` (i.e., with `deltj = 0`). The density matrix is computed as: 
 * `density_matrix = imag_i * (Bose/Fermi) * lesser(t_lead,t_lead)`.
 * This function is intended for cases where the matrix size is 1×1.
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @return The computed density matrix as a complex number.
 *
 */
  template <typename T>
  std::complex<T> herm_matrix_timestep_moving_view<T>::density_matrix(void){
    CPLX x1;
    get_les(0,x1);
    return CPLX(0.0,1.0*sig_)*x1;
  }
/**
 * @brief Computes the density matrix from the lesser Green's function component.
 *
 * This function retrieves the lesser Green's function component at the leading time step 
 * `t_lead` (i.e., with `deltj = 0`). The density matrix is computed as: 
 * `density_matrix = imag_i * (Bose/Fermi) * lesser(t_lead,t_lead)`. 
 * The result is stored in the provided matrix `M`.
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @tparam Matrix Type of the output matrix, which must support element-wise multiplication.
 * @param M Output matrix where the computed density matrix will be stored.
 *
 * @post The matrix `M` is updated with the computed density matrix.
 */
  template <typename T> template <class Matrix>
  void herm_matrix_timestep_moving_view<T>::density_matrix(Matrix &M){
    get_les(0,M);
    M *=CPLX(0.0,1.0*sig_);   
  }

/**
 * @brief Scales the Green's function components by a given scalar factor.
 *
 * This function multiplies all elements of the retarded and lesser Green's function 
 * components by the provided scalar `weight`.
 * `this^[X](t_lead,t_lead-j)` set to `eight*this^[X](t_lead,t_lead-j)` for `j=0...tc` and `X=ret,les`
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param weight The scalar factor by which all elements are multiplied.
 *
 */
template <typename T>
void herm_matrix_timestep_moving_view<T>::smul(T weight) {
    int m;
    CPLX *x0;
    x0 = ret_;
    for (m = 0; m <= tc_; m++) {
      element_smul<T, LARGESIZE>(size1_, x0 + m * element_size_,
				 weight);
    }
    x0 = les_;
    for (m = 0; m <= tc_; m++) {
      element_smul<T, LARGESIZE>(size1_, x0 + m * element_size_,
				 weight);
    }
    
}

/**
 * @brief Scales the Green's function components by a given scalar factor.
 *
 * This function multiplies all elements of the retarded and lesser Green's function 
 * components by the provided scalar `weight`.
 * `this^[X](t_lead,t_lead-j)` set to `weight*this^[X](t_lead,t_lead-j)` for `j=0...tc` and `X=ret,les`
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param weight The scalar factor by which all elements are multiplied.
 *
 * @post Both the retarded and lesser Green's function components are scaled by `weight`.
 */

template <typename T>
void herm_matrix_timestep_moving_view<T>::smul(std::complex<T> weight) {
    int m;
    CPLX *x0;
    x0 = ret_;
    for (m = 0; m <= tc_; m++) {
      element_smul<T, LARGESIZE>(size1_, x0 + m * element_size_,
				 weight);
    }
    x0 = les_;
    for (m = 0; m <= tc_; m++) {
      element_smul<T, LARGESIZE>(size1_, x0 + m * element_size_,
				 weight);
    }
    
}

  
/**
 * @brief Left-multiplies  the Green's function timeslice by a time-dependent contour function.
 *
 * This function performs the operation \f$ G(t,t') \rightarrow F(t) G(t,t') \f$
 * at the leading time step, where `F(t)` is a time-dependent contour function.
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param ft The time-dependent contour function `F(t)`, which must have the same time cutoff (`tc_`).
 *
 * The time cutoff `tc_` of `ft` must match `tc_` of the Green’s function.
 */
template <typename T>
void herm_matrix_timestep_moving_view<T>::left_multiply(function_moving<T> &ft) {
   assert( ft.tc() == tc_);
    int m;
    cplx *xtemp, *ftemp, *x0;
    xtemp = new cplx[element_size_];
    ftemp = ft.ptr(0);
    x0 = retptr(0);
    for (m = 0; m <= tc_; m++) {
      element_mult<T, LARGESIZE>(size1_, xtemp, ftemp, x0 + m * element_size_);
      element_set<T, LARGESIZE>(size1_, x0 + m * element_size_, xtemp);
    }
    x0 = lesptr(0);
    for (m = 0; m <= tc_; m++) {
      element_mult<T, LARGESIZE>(size1_, xtemp, ftemp,x0 + m * element_size_);
      element_set<T, LARGESIZE>(size1_, x0 + m * element_size_, xtemp);
    }    
    delete[] xtemp;
    
}
/**
 * @brief Left-multiplies the Green's function timeslice by the hermitian conjugate of a  time-dependent contour function.
 *
 * This function performs the operation \f$ G(t,t') \rightarrow F(t)^\dagger G(t,t') \f$
 * at the leading time step, where `F(t)` is a time-dependent contour function.
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param ft The time-dependent contour function `F(t)`, which must have the same time cutoff (`tc_`).
 *
 * The time cutoff `tc_` of `ft` must match `tc_` of the Green’s function.
 */
template <typename T>
void herm_matrix_timestep_moving_view<T>::left_multiply_hermconj(function_moving<T> &ft) {
   assert( ft.tc() == tc_);
    int m;
    cplx *xtemp, *x0, *fcc;
    xtemp = new cplx[element_size_];
    fcc = new cplx[element_size_];    
    element_conj<T, LARGESIZE>(size1_, fcc, ft.ptr(0));    
    x0 = retptr(0);
    for (m = 0; m <= tc_; m++) {
      element_mult<T, LARGESIZE>(size1_, xtemp, fcc, x0 + m * element_size_);
      element_set<T, LARGESIZE>(size1_, x0 + m * element_size_, xtemp);
    }
    x0 = lesptr(0);
    for (m = 0; m <= tc_; m++) {
      element_mult<T, LARGESIZE>(size1_, xtemp, fcc,x0 + m * element_size_);
      element_set<T, LARGESIZE>(size1_, x0 + m * element_size_, xtemp);
    }    
    delete[] xtemp;
    delete[] fcc;

}

/**
 * @brief Right-multiplies  the Green's function timeslice by a time-dependent contour function.
 *
 * This function performs the operation \f$ G(t,t') \rightarrow G(t,t') F(t')  \f$
 * at the leading time step, where `F(t')` is a time-dependent contour function
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param ft The time-dependent contour function `F(t')`, which must have a time cutoff (`tc_`) 
 *           greater than or equal to that of the Green's function.
 *
 * The time cutoff `tc_` of `ft` must be at least as large as `tc_` of the Green’s function.
 */
template <typename T>
void herm_matrix_timestep_moving_view<T>::right_multiply(function_moving<T> &ft) {
   assert( ft.tc() >= tc_);
    int m;
    cplx *xtemp, *x0,*ftemp;
    ftemp=ft.ptr(0);
    xtemp = new cplx[element_size_];
    x0 = retptr(0);
    for (m = 0; m <= tc_; m++) {
      element_mult<T, LARGESIZE>(size1_, xtemp, x0 + m * element_size_,ftemp + m * element_size_);
      element_set<T, LARGESIZE>(size1_, x0 + m * element_size_, xtemp);
    }
    x0 = lesptr(0);
    for (m = 0; m <= tc_; m++) {
      element_mult<T, LARGESIZE>(size1_, xtemp, x0 + m * element_size_,ftemp + m * element_size_);
      element_set<T, LARGESIZE>(size1_, x0 + m * element_size_, xtemp);
    }
    delete[] xtemp;
}


/**
 * @brief Right-multiplies the Green's function timeslice by the conjugate of time-dependent contour function.
 *
 * This function performs the operation \f$ G(t,t') \rightarrow G(t,t') F(t')^\dagger  \f$
 * at the leading time step, where `F(t')` is a time-dependent contour function. 
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param ft The time-dependent contour function `F(t')`
 *
 * The time cutoff `tc_` of `ft` must be at least as large as `tc_` of the Green’s function.
 */
template <typename T>
void herm_matrix_timestep_moving_view<T>::right_multiply_hermconj(function_moving<T> &ft) {
   assert( ft.tc() >= tc_);
    int m;
    cplx *xtemp, *x0,*fcc;
    fcc = new cplx[element_size_];    
    xtemp = new cplx[element_size_];

    x0 = retptr(0);
    for (m = 0; m <= tc_; m++) {
      element_conj<T, LARGESIZE>(size1_, fcc, ft.ptr(m));    
      element_mult<T, LARGESIZE>(size1_, xtemp, x0 + m * element_size_,fcc);
      element_set<T, LARGESIZE>(size1_, x0 + m * element_size_, xtemp);
    }
    x0 = lesptr(0);
    for (m = 0; m <= tc_; m++) {
      element_conj<T, LARGESIZE>(size1_, fcc, ft.ptr(m));          
      element_mult<T, LARGESIZE>(size1_, xtemp, x0 + m * element_size_,fcc);
      element_set<T, LARGESIZE>(size1_, x0 + m * element_size_, xtemp);
    }
    delete[] xtemp;
    delete[] fcc;
}



#define HERM_MATRIX_INCR_TSTP                                                \
    if (alpha == CPLX(1.0, 0.0)) {                                           \
        for (i = 0; i < len; i++)                                            \
            x0[i] += x[i];                                                   \
    } else {                                                                 \
        for (i = 0; i < len; i++)                                            \
            x0[i] += alpha * x[i];                                           \
    }

/**
 * @brief Increments the Green's function by another Green's function at the given timestep, scaled by a factor.
 *
 * This function performs the operation:
 * `this^[X](t_lead,t_lead-j) += \alpha * g(t_lead,t_lead-j)` for `j=0...tc` and `X=ret,les`
 * where g is another Green’s function of the same dimensions, and `\alpha` is a complex scaling factor.
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param g The Green’s function `g(t,t')` to be added.
 * @param alpha the complex scaling factor
 *
 * `g` must have the same `tc_`, `size1_`, and `size2_` as this object.
 */
template <typename T>
void herm_matrix_timestep_moving_view<T>::incr_timestep(herm_matrix_timestep_moving_view<T> &g,cplx alpha) {
  assert(tc_==g.tc() && size1_ == g.size1() && size2_ == g.size2());

    int i, len;
    CPLX *x, *x0;
    len = (tc_ + 1) * element_size_;
    x0 = ret_;
    x = g.ret_;
    HERM_MATRIX_INCR_TSTP
    x0 = les_;
    x = g.les_;
    HERM_MATRIX_INCR_TSTP
}
#undef HERM_MATRIX_INCR_TSTP

/**
 * @brief Increments the Green's function by another Green's function at the given timestep, scaled by a factor.
 *
 * This function performs the operation:
 * `this^[X](t_lead,t_lead-j) += \alpha * g(t_lead,t_lead-j)` for `j=0...tc` and `X=ret,les`
 * where g is another Green’s function of the same dimensions, and `\alpha` is a complex scaling factor.
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param g The Green’s function `g(t,t')` to be added.
 * @param alpha the complex scaling factor
 *
 * `g` must have the same `tc_`, `size1_`, and `size2_` as this object.

 */
template <typename T>
void herm_matrix_timestep_moving_view<T>::incr_timestep(herm_matrix_timestep_moving<T> &g, cplx alpha) {
    herm_matrix_timestep_moving_view<T> tmp(g);
    incr_timestep(tmp, alpha);
}

/**
 * @brief Increments the Green's function by timeslide deltG of another Green's function at the given timestep, scaled by a factor.
 *
 * This function performs the operation:
 * `this^[X](t_lead,t_lead-j) += \alpha * g(t_lead-deltG,t_lead-deltG-j)` for `j=0...tc` and `X=ret,les`
 * where g is another Green’s function of the same dimensions, and `\alpha` is a complex scaling factor.
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param g The Green’s function `g(t,t')` to be added.
 * @param alpha the complex scaling factor
 *
 * `g` must have the same `tc_`, `size1_`, and `size2_` as this object.
 */
template <typename T>
void herm_matrix_timestep_moving_view<T>::incr_timestep(herm_matrix_moving<T> &g,int deltG, cplx alpha) {
  herm_matrix_timestep_moving_view<T> tmp(g,deltG);
  incr_timestep(tmp, alpha);
}

  
/**
 * @brief Copies the data from a physical `herm_matrix_timestep_view` object.
 *
 * This function copies the retarded and lesser Green’s function components 
 * form the Greenfunction G an its physical timestep tstp
 * `this^[X](t_lead,t_lead-j)` set to `g^[X](tstp,tstp-j)` for all `j=0....min(tc,tstp)` and `X=ret,les`
 * if `tstp<tc`, the remaining data `this(t_lead,t_lead-j), j=tstp+1,...,tc` are set to zero
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param g Source `herm_matrix_timestep_view` object, pointing to timestep tstp of a Greenfunction
 * @param gcc hermitioan conjugate of g. Use g=gcc if g is hermitian
 *
 * The matrix sizes (`size1_`) must match.
 */
  template <typename T>
void herm_matrix_timestep_moving_view<T>::set_timestep(herm_matrix_timestep_view<T> &g,herm_matrix_timestep_view<T> &gcc){
    assert(size1_==g.size1());
    int t1;
    int tstp=g.tstp();
    clear_timestep();
    int smax=(tstp > tc_ ? tc_ : tstp); // note that at this point t>=0
    for(t1=0;t1<=smax;t1++){
    //for(t1=0;t1<=tc_;t1++){
      element_set<T,LARGESIZE>(size1_,retptr(t1),g.retptr(tstp-t1));
      element_minusconj<T,LARGESIZE>(size1_,lesptr(t1),gcc.lesptr(tstp-t1));
    }  
}

/**
 * @brief Copies the data from a physical `herm_matrix_timestep` object.
 *
 * This function copies the retarded and lesser Green’s function components 
 * form the Greenfunction G an its physical timestep tstp.
 * `this^[X](t_lead,t_lead-j)` set to `g^[X](tstp,tstp-j)` for all `j=0....min(tc,tstp)` and `X=ret,les`.
 * If `tstp<tc`, the remaining data `this(t_lead,t_lead-j), j=tstp+1,...,tc` are set to zero
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param g Source `herm_matrix_timestep_view` object, pointing to timestep tstp of a Greenfunction
 * @param gcc hermitian conjugate of g. Use g=gcc if g is hermitian
 *
 * The matrix sizes (`size1_`) must match. */
  template <typename T>
void herm_matrix_timestep_moving_view<T>::set_timestep(herm_matrix_timestep<T> &g,herm_matrix_timestep<T> &gcc){
    int tstp=g.tstp();
    herm_matrix_timestep_view<T> tmp(tstp,g);
    herm_matrix_timestep_view<T> tmp1(tstp,gcc);
    this->set_timestep(tmp,tmp1);
  }
  
/**
  * @brief Copies the data from timeslice tstp of a a physical `herm_matrix` object.
 *
 * This function copies the retarded and lesser Green’s function components 
 * form the Greenfunction G an its physical timestep tstp.
 * `this^[X](t_lead,t_lead-j)` set to `g^[X](tstp,tstp-j)` for all `j=0....min(tc,tstp)` and `X=ret,les`.
 * If `tstp<tc`, the remaining data t`his(t_lead,t_lead-j), j=tstp+1,...,tc` are set to zero
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param g Source `herm_matrix_timestep_view` object, pointing to timestep tstp of a Greenfunction
 * @param gcc hermitian conjugate of g. Use g=gcc if g is hermitian
 * @param tstp the physical timeslice of g. Must mach 0<=tstp<=g.nt()
 *
 * The matrix sizes (`size1_`) must match.
 */
  template <typename T>
void herm_matrix_timestep_moving_view<T>::set_timestep(herm_matrix<T> &g,herm_matrix<T> &gcc,int tstp){
    assert(0<=tstp && tstp<=g.nt());
    herm_matrix_timestep_view<T> tmp(tstp,g);
    herm_matrix_timestep_view<T> tmp1(tstp,gcc);
    this->set_timestep(tmp,tmp1);
}
    
#if CNTR_USE_MPI ==1
  /// @private
  template <typename T>
  void my_mpi_reduce_moving(std::complex<T> *data, int len, int root){
    std::cerr << __PRETTY_FUNCTION__ << ", LEN=" << len
	      << " ... NOT DEFINED FOR THIS TYPE " << std::endl;
    exit(0);
  }
  /// @private
  template<>
  inline void my_mpi_reduce_moving<double>(std::complex<double> *data, int len, int root) {
    int tid,ntasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &tid);
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
    assert(root>=0 && root <= ntasks -1);
    assert(len>=0);
    if (tid == root) {
      MPI_Reduce(MPI_IN_PLACE, (double *)data, 2 * len,
		 MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
    } else {
      MPI_Reduce((double *)data, (double *)data, 2 * len,
		 MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
    }
  }
  /// @private
  template <typename T>
  void herm_matrix_timestep_moving_view<T>::MPI_Reduce(int root) {
    my_mpi_reduce_moving<T>(les_, (tc_ + 1) * element_size_, root);
    my_mpi_reduce_moving<T>(ret_, (tc_ + 1) * element_size_, root);
  }  
#endif


  #if CNTR_USE_HDF5 == 1
  
  template <typename T>
  void herm_matrix_timestep_moving_view<T>::write_to_hdf5(hid_t group_id) {
    store_int_attribute_to_hid(group_id, std::string("herm_matrix_timestep_moving"), 1);
    store_int_attribute_to_hid(group_id, std::string("tc"), tc_);
    store_int_attribute_to_hid(group_id, std::string("size1"), size1_);
    store_int_attribute_to_hid(group_id, std::string("size2"), size2_);
    store_int_attribute_to_hid(group_id, std::string("element_size"),element_size_);
    store_int_attribute_to_hid(group_id, std::string("sig"),sig_);
    hsize_t len_shape = 3, shape[3],start[3];
    shape[1] = size1_;
    shape[2] = size2_;
    if (tc_==-1) {
      // continue;
    }else{
      shape[0] = (tc_+1)*element_size_;
      //store_cplx_array just to init datastructure
      store_cplx_array_to_hid(group_id, std::string("ret"), retptr(0), shape, len_shape);
      store_cplx_array_to_hid(group_id, std::string("les"), lesptr(0), shape, len_shape);
    }
  }
  template <typename T>
  void herm_matrix_timestep_moving_view<T>::write_to_hdf5(hid_t group_id, const char *groupname) {
    hid_t sub_group_id = create_group(group_id, groupname);
    this->write_to_hdf5(sub_group_id);
    close_group(sub_group_id);
  }
  template <typename T>
  void herm_matrix_timestep_moving_view<T>::write_to_hdf5(const char *filename,
					    const char *groupname, h5_mode mode) {
    hid_t file_id = open_hdf5_file(filename,mode);
    this->write_to_hdf5(file_id, groupname);
    close_hdf5_file(file_id);
  }
  
#endif


  

#undef CPLX

  
}// namespace cntr
#endif
