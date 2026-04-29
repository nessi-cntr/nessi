#ifndef GF_IMPL
#define GF_IMPL
#define EIGEN_MULT

#include "ness.hpp"


namespace ness
{

/** \brief <b> Set a submatrix of the Green's function equal to a smaller Green's function</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > Set a submatrix of this `GF` Green's function \f$ A \f$, specified by multiindices `subid1`, `subid2`, equal to another Green's function \f$ B\f$
* > which has the submatrix size, for all times /frequencies \f$ t \f$ and for lesser and retarded component.
* > \f$ A^{R,<}_{{\tt subid1},{\tt subid2}}(t) = B^{R,<}(t) \f$
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param  subid1
* [std::vector<int>] Vector of integers specifying the rows of the submatrix
* @param  subid2
* [std::vector<int>] Vector of integers specifying the columns of the submatrix
* @param  gf
* [GF] Smaller Green's function which is used to set the submatrix values of the larger one
*/
	int GF::set_submatrix(const std::vector<int> & subid1, const std::vector<int> & subid2, const GF& gf)
	{
		CNTR_ASSERT_EQ(CNTR_ASSERT_LEVEL_0, gf.size1_, subid1.size(),__PRETTY_FUNCTION__);
		CNTR_ASSERT_EQ(CNTR_ASSERT_LEVEL_0, gf.size2_, subid2.size(),__PRETTY_FUNCTION__);
		CNTR_ASSERT_EQ(CNTR_ASSERT_LEVEL_0, this -> ngrid_, gf.ngrid_,__PRETTY_FUNCTION__);

		for(int w = 0; w < this -> ngrid_; w++)
		{
			for(int i1 = 0; i1 < gf.size1_; i1++)
			{
				for(int i2 = 0; i2 < gf.size2_; i2++)
				{
					g_ret[w * el_size_ + subid1[i1] * size2_ + subid2[i2]] = gf.g_ret[w * gf.el_size_ + i1 * gf.size2_ + i2];
					g_les[w * el_size_ + subid1[i1] * size2_ + subid2[i2]] = gf.g_les[w * gf.el_size_ + i1 * gf.size2_ + i2];
				}
			}
		}
		return 0;
	}

/** \brief <b> Set the Green's function equal to a submatrix of a larger Green's function</b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Set this `GF` Green's function \f$ A \f$ equal to a submatrix of a larger `GF` Green's function \f$ B\f$.
 * > The submatrix, specified by multiindices `subid1`, `subid2`, has the size of the smaller Green's function and is set for all times /frequencies \f$ t \f$ as well as for retarded and lesser component.
 * > \f$ A^{R,<}(t) = B^{R,<}_{{\tt subid1},{\tt subid2}}(t) \f$
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param  subid1
 * [std::vector<int>] Vector of integers specifying the rows of the submatrix
 * @param  subid2
 * [std::vector<int>] Vector of integers specifying the columns of the submatrix
 * @param  gf
 * [GF] Larger Green's function to extract the submatrix from and put it into the smaller one
 */
	int GF::get_submatrix(const std::vector<int> & subid1, const std::vector<int> & subid2, const GF& gf)
	{
		CNTR_ASSERT_EQ(CNTR_ASSERT_LEVEL_0, size1_, subid1.size(),__PRETTY_FUNCTION__);
		CNTR_ASSERT_EQ(CNTR_ASSERT_LEVEL_0, size2_, subid2.size(),__PRETTY_FUNCTION__);
		CNTR_ASSERT_EQ(CNTR_ASSERT_LEVEL_0, this -> ngrid_, gf.ngrid_,__PRETTY_FUNCTION__);

		for(int w = 0; w < this -> ngrid_; w++)
		{
			for(int i1 = 0; i1 < size1_; i1++)
			{
				for(int i2 = 0; i2 < size2_; i2++)
				{
					g_ret[w * el_size_ + i1 * size2_ + i2] = gf.g_ret[w * gf.el_size_ + subid1[i1] * gf.size2_ + subid2[i2]];
					g_les[w * el_size_ + i1 * size2_ + i2] = gf.g_les[w * gf.el_size_ + subid1[i1] * gf.size2_ + subid2[i2]];
				}
			}
		}
		return 0;
	}

/** \brief <b> Set an element of the Green's function equal to a scalar valued Green's function</b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Set a single element `(i1,i2)` of this `GF` Green's function \f$ A \f$ equal to a scalar valued (`size=1`) `GF` Green's function \f$ B \f$ for all times /frequencies \f$ t\f$.
 * > \f$ A^{R,<}_{{\tt i1},{\tt i2}}(t) = B^{R,<}(t) \f$
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param i1
 *  [int] Row index of element to be set
 * @param i2
 *  [int] Column index of element to be set
 * @param  gf
 *  [GF] Scalar valued Green's function
 */
	int GF::set_element(int i1, int i2, const GF& gf)
	{
		CNTR_ASSERT_EQ(NESS_ASSERT_0, gf.size1_, 1, __PRETTY_FUNCTION__);
		CNTR_ASSERT_EQ(NESS_ASSERT_0, gf.size2_, 1, __PRETTY_FUNCTION__);
		CNTR_ASSERT_EQ(NESS_ASSERT_0, gf.ngrid_, this -> ngrid_, __PRETTY_FUNCTION__);
		for(int i = 0; i < ngrid_; i++)
		{
			g_ret[i * el_size_ + i1 * size2_ + i2] = gf.g_ret[i];
			g_les[i * el_size_ + i1 * size2_ + i2] = gf.g_les[i];
		}
		return 0;
	}

/** \brief <b> Set an element of the Green's function equal to an element of another Green's function</b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Set an element `(i1,i2)` of this `GF` Green's function \f$ A \f$ equal to an element `(s1,s2)` of another `GF` Green's function \f$ B \f$ for all times /frequencies \f$ t\f$ and retarded and lesser component.
 * > \f$ A^{R,<}_{{\tt i1},{\tt i2}}(t) = B^{R,<}_{{\tt s1},{\tt s2}}(t) \f$
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param i1
 *  [int] Row index of element to be set
 * @param i2
 *  [int] Column index of element to be set
 * @param  gf
 *  [GF] Green's function to take element from
 * @param s1
 *  [int] Row index of element to extracted
 * @param s2
 *  [int] Column index of element to be extracted
 */
	int GF::set_element(int i1, int i2, const GF& gf, int s1, int s2)
	{
		CNTR_ASSERT_EQ(NESS_ASSERT_0, gf.ngrid_, this -> ngrid_, __PRETTY_FUNCTION__);
		for(int i = 0; i < ngrid_; i++)
		{
			g_ret[i * el_size_ + i1 * size2_ + i2] = gf.g_ret[i * gf.el_size_ + s1 * gf.size2_ + s2];
			g_les[i * el_size_ + i1 * size2_ + i2] = gf.g_les[i * gf.el_size_ + s1 * gf.size2_ + s2];
		}
		return 0;
	}

/** \brief <b> Set a scalar valued Green's function equal to an element of the Green's function </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Set a scalar valued (`size=1`) `GF` Green's function \f$ B \f$ equal to a single element `(i1,i2)` of this `GF` Green's function \f$ A \f$ for all times /frequencies \f$ t\f$ and retarded and lesser component.
 * > \f$ B^{R,<}(t) = A^{R,<}_{{\tt i1},{\tt i2}}(t) \f$
 *
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param i1
 * [int] Row index of element to be extracted
 * @param i2
 * [int] Column index of element to be extracted
 * @param  gf
 * [GF] Scalar valued Green's function to be set
 */
	int GF::get_element(int i1, int i2, GF& gf) const
	{
		CNTR_ASSERT_EQ(NESS_ASSERT_0, gf.size1_, 1, __PRETTY_FUNCTION__);
		CNTR_ASSERT_EQ(NESS_ASSERT_0, gf.size2_, 1, __PRETTY_FUNCTION__);
		CNTR_ASSERT_EQ(NESS_ASSERT_0, gf.ngrid_, this -> ngrid_, __PRETTY_FUNCTION__);
		for(int i = 0; i < ngrid_; i++)
		{
			gf.g_ret[i] = g_ret[i * el_size_ + i1 * size2_ + i2];
			gf.g_les[i] = g_les[i * el_size_ + i1 * size2_ + i2];
		}
		return 0;
	}


/** \brief <b> Set an element of a Green's function object equal to an element of the Green's function</b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Set a single element `(s1,s2)` of another `GF` Green's function \f$ B \f$  equal to a single element
 * > `(i1,i2)` of this `GF` Green's function \f$ A \f$  for all times /frequencies \f$ t \f$ and for retarded and lesser component.
 * > \f$ B^{R,<}_{{\tt s1},{\tt s2}}(t) = A^{R,<}_{{\tt i1},{\tt i2}}(t) \f$
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param i1
 *  [int] Row index of element to be set
 * @param i2
 *  [int] Column index of element to be set
 * @param  gf
 *  [GF] Green's function to take element from
 * @param s1
 *  [int] Row index of element to extracted
 * @param s2
 *  [int] Column index of element to be extracted
 */
    int GF::get_element(int i1, int i2, const GF& gf, int s1, int s2)
    {
        CNTR_ASSERT_EQ(NESS_ASSERT_0, gf.ngrid_, this -> ngrid_, __PRETTY_FUNCTION__);
        for(int i = 0; i < ngrid_; i++)
        {
            gf.g_ret[i * gf.el_size_ + s1 * gf.size2_ + s2] = g_ret[i * el_size_ + i1 * size2_ + i2];
            gf.g_les[i * gf.el_size_ + s1 * gf.size2_ + s2] = g_les[i * el_size_ + i1 * size2_ + i2];
        }
        return 0;
    }

/** \brief <b> Clear all data of a Green's function </b> */
	int GF::clear()
	{
		std::memset(g_ret, 0, sizeof(cplx) * data_size_);
		std::memset(g_les, 0, sizeof(cplx) * data_size_);
		return 0;
	}

/** \brief <b> Multiplies the Green's function with a matrix from the left </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Multiplies the `GF` Green's function \f$ A \f$  with a matrix \f$ g \f$  from the left for all times /frequencies \f$ t \f$ and for retarded and lesser component,
 * > matrix structure has to match.
 * > \f$ A^{r,<}(t) \to A^{r,<}(t) \cdot g \f$
 *
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param g
 * [cdmatrix] Matrix to multiply the Green's function with
 */
	int GF::left_multiply(const cdmatrix &g)
	{
		CNTR_ASSERT_EQ(NESS_ASSERT_0, g.rows(), size1_, __PRETTY_FUNCTION__);
		CNTR_ASSERT_EQ(NESS_ASSERT_0, g.cols(), size2_, __PRETTY_FUNCTION__);
        
		for (int w = 0; w < ngrid_; w++)
		{
			Retarded[w] = (Retarded[w] * g).eval();
			Lesser[w] = (Lesser[w] * g).eval();
		}
		return 0;
	}

/** \brief <b> Multiplies the Green's function with a matrix from the right </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Multiplies the `GF` Green's function \f$ A \f$  with a matrix \f$ g \f$  from the right for all times /frequencies \f$ t \f$ and for retarded and lesser component,
 * > matrix structure has to match.
 * > \f$ A^{r,<}(t) \to g \cdot A^{r,<}(t) \f$
 *
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param g
 * [cdmatrix] Matrix to multiply the Green's function with
 */
	int GF::right_multiply(const cdmatrix &g)
	{
		for (int w = 0; w < ngrid_; w++)
		{
			Retarded[w] = (g * Retarded[w]).eval();
			Lesser[w] = (g * Lesser[w]).eval();
		}
		return 0;
	}

/** \brief <b> Shift the Green's function regarding its argument </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Shift the `GF` Green's function \f$ A \f$ by a specified amount of grid points `shift`, total energy shift is \f$ {\tt shift} \cdot d\omega \f$ (both retarded and lesser component).
 * >\f$ A^{R,<}(t) \to A^{R,<}(t+{\tt shift}) \f$
 * \Note Should only be applied to frequency Green's functions, considers order of frequencies.
 *
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param shift
 * [cdmatrix] Amount to shift the frequency index of the GF
 */
	void GF::set_shiftGF(int shift)
	{
		cdmatrix zerom = Eigen::MatrixXcd::Zero(size1_, size2_);
		int ashift = std::abs(shift);
		if (shift > 0)
		{
			for(int w = ngrid_ / 2 + 1; w < ngrid_ - shift; w++)
			{
				Retarded[w] = Retarded[w + shift];
				Lesser[w] = Lesser[w + shift];
			}
			for(int w = ngrid_ - shift; w < ngrid_; w++)
			{
				Retarded[w] = Retarded[w - ngrid_ + shift];
				Lesser[w] = Lesser[w - ngrid_ + shift];
			}

			for(int w = 0; w < ngrid_ / 2 + 1 - shift; w++)
			{
				Retarded[w] = Retarded[w + shift];
				Lesser[w] = Lesser[w + shift];
			}
			for(int w = ngrid_ / 2 - shift + 1; w < ngrid_ / 2 + 1; w++)
			{
				Retarded[w] = zerom;
				Lesser[w] = zerom;
			}



		}
		else if (shift < 0)
		{
			//note shift is neg.

			for(int w = ngrid_ / 2; w >= ashift; w--)
			{
				Retarded[w] = Retarded[w - ashift];
				Lesser[w] = Lesser[w - ashift];
			}
			for(int w = ashift - 1; w >= 0; w--)
			{
				Retarded[w] = Retarded[ngrid_ + w - ashift];
				Lesser[w] = Lesser[ngrid_ + w - ashift];
			}

			for(int w = ngrid_ - 1; w >= ngrid_ / 2 + 1 + ashift; w--)
			{
				Retarded[w] = Retarded[w - ashift];
				Lesser[w] = Lesser[w - ashift];
			}
			for(int w = ngrid_ / 2 + ashift ; w >= ngrid_ / 2 + 1; w--)
			{
				Retarded[w] = zerom;
				Lesser[w] = zerom;
			}


		}
	}

/** \brief <b> Destructor </b> */
	GF::~GF(){
		delete[] g_ret;
		delete[] g_les;
		delete[] grid_;
	}

/** \brief <b> Default constructor </b>*/
	GF::GF(){
		dgrid_ = 0.0;
		ngrid_ = 0;
		size1_ = 0;
		size2_ = 0;
		el_size_ = 0;
		data_size_ = 0;
		g_ret = NULL;
		g_les = NULL;
		grid_ = NULL;
		gf_type_ = 0;
	}

/** \brief <b> Adds another Green's function multiplied by a factor to the Green's function </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Adds another `GF` Green's function \f$ B \f$ multiplied by a complex weight factor `weight` to this `GF` Green's function \f$ A \f$ (both retarded and lesser component).
 * >  For all times /frequencies \f$ t \f$, matrix structure has to match.
 * > \f$ A^{R,<} \to A^{R,<} + {\tt weight} \cdot B^{R,<} \f$
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param gf
 * [GF] Green's function for incrementation
 * @param weight
 * [cplx] complex weight factor
 */
	int GF::incr(GF &gf, cplx weight)
	{
		assert(this -> gf_type_ == gf.gf_type_);
		assert(this -> el_size_ == gf.el_size_);
		assert(this -> data_size_ == gf.data_size_);
		for (int idx = 0; idx < data_size_; idx++)
		{
			g_ret[idx] += gf.g_ret[idx] * weight;
			g_les[idx] += gf.g_les[idx] * weight;
		}
		return 0;
	}

/** \brief <b> Multiplies the Green's function  by a factor </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Multiplies the `GF` Green's function \f$ A \f$ by a complex factor  `weight` for all times / frequencies and for the retarded and lesser component.
 * > \f$ A^{R,<} \to {\tt weight} \cdot A^{R,<} \f$
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param weight
 * [cplx] complex  factor
*/
	int GF::smul(cplx weight)
	{
		for (int idx = 0; idx < data_size_; idx++)
		{
			g_ret[idx] *= weight;
			g_les[idx] *= weight;
		}
		return 0;
	}

/** \brief <b> Copy constructor </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > Constructs an empty `GF` with all parameters given by another `GF`,
* > so that the new object has the same shape and type as the parameter `GF`.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param g
* [GF] Green's function to copy parameters from
*/
	GF::GF(const GF &g)
	{
		size1_ = g.size1_;
		size2_ = g.size2_;
		el_size_ = g.el_size_;
		dgrid_ = g.dgrid_;
		ngrid_ = g.ngrid_;
		gf_type_ = g.gf_type_;
		data_size_ = g.data_size_;
		sign_ = g.sign_;

        if (data_size_ > 0){
            g_ret = new cplx[data_size_];
            g_les = new cplx[data_size_];
            grid_ = new double[ngrid_];
            std::memcpy(g_ret, g.g_ret, sizeof(cplx) * data_size_);
            std::memcpy(g_les, g.g_les, sizeof(cplx) * data_size_);
            std::memcpy(grid_, g.grid_, sizeof(double) * ngrid_);
        }else{
            g_ret = 0;
            g_les = 0;
            grid_ = 0;
        }
		Retarded.clear();
		Lesser.clear();

		Retarded.reserve(ngrid_);
		Lesser.reserve(ngrid_);
		for(int w = 0; w < ngrid_; w++)
		{
			Retarded.emplace_back(Eigen::Map<Eigen::Matrix<cplx, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > (p_ret(w,0,0), size1_, size2_));
			Lesser.emplace_back(Eigen::Map<Eigen::Matrix<cplx, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > (p_les(w,0,0), size1_, size2_));

		}
	}

/** \brief <b> Move  constructor </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > Constructs a `GF` with all parameters and contents moved from another `GF`.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param g
* [GF] Green's function
*/
#if __cplusplus > 199711L
	GF::GF(GF &&g)
	{
		size1_ = g.size1_;
		size2_ = g.size2_;
		el_size_ = g.el_size_;
		dgrid_ = g.dgrid_;
		ngrid_ = g.ngrid_;
		gf_type_ = g.gf_type_;
		data_size_ = g.data_size_;
		sign_ = g.sign_;

		g_ret = g.g_ret;
		g.g_ret = nullptr;
		g_les = g.g_les;
		g.g_les = nullptr;
		grid_ = g.grid_;
		g.grid_ = nullptr;

		Retarded = std::move(g.Retarded);
		Lesser = std::move(g.Lesser);
	}
#endif

/** \brief <b> Operator = for Green's functions </b> */
	GF& GF::operator=(const GF &g)
	{
		size1_ = g.size1_;
		size2_ = g.size2_;
		el_size_ = g.el_size_;
		dgrid_ = g.dgrid_;
		gf_type_ = g.gf_type_;
		sign_ = g.sign_;
		if (this -> ngrid_ != g.ngrid_)
		{
			if (ngrid_ != 0)
				delete[] grid_;
			ngrid_ = g.ngrid_;
			grid_ = new double[ngrid_];

		}
		Retarded.clear();
		Lesser.clear();
		Retarded.reserve(ngrid_);
		Lesser.reserve(ngrid_);

		if (this -> data_size_!= g.data_size_) 
		{
			if (data_size_ != 0)
			{
				delete[] g_ret;
				delete[] g_les;
			}
			data_size_ = g.data_size_;
			g_ret = new cplx[data_size_];
			g_les = new cplx[data_size_];
		}
		std::memcpy(g_ret, g.g_ret, sizeof(cplx) * data_size_);
		std::memcpy(g_les, g.g_les, sizeof(cplx) * data_size_);
		std::memcpy(grid_, g.grid_, sizeof(double) * ngrid_);
		for(int w = 0; w < ngrid_; w++)
		{
			Retarded.emplace_back(Eigen::Map<Eigen::Matrix<cplx, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > (p_ret(w,0,0), size1_, size2_));
			Lesser.emplace_back(Eigen::Map<Eigen::Matrix<cplx, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > (p_les(w,0,0), size1_, size2_));

		}


		return *this;
	}

/** \brief <b> Standard constructor for a `GF` object </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Constructs a `GF` with all grid (spacing `dgrid`, length `ngrid`), shape (orbital size `size1 x size2`) and type (fermion/boson `sign`, time/freqeuncy `gf_type`) parameters specified manually.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param dgrid
 * [double] Grid spacing
 * @param ngrid
 * [int] Number of grid points
 * @param size1
 * [int] First orbital dimension, `GF` matrix rows
 * @param size2
 * [int] Second orbital dimension, `GF` matrix columns
 * @param sign
 * [int] Sign of particle type, +1 for boson,-1 for fermion.
 * @param gf_type
 * [int] `GF` type, 0: time, 1: frequency
 */
	GF::GF(double dgrid, int ngrid, int size1, int size2, int sign, int gf_type):
		size1_(size1), size2_(size2), sign_(sign),
		dgrid_(dgrid), ngrid_(ngrid), gf_type_(gf_type)
	{
		Retarded.reserve(ngrid_);
		Lesser.reserve(ngrid_);
		int halfn = ngrid / 2 + 1;
		el_size_ = size1_ * size2_;
		grid_ = new double[ngrid_];
		reset_grid(dgrid);
		data_size_ = ngrid_ * el_size_;
		g_ret = new cplx[data_size_];
		g_les = new cplx[data_size_];
		std::memset(g_ret, 0, sizeof(cplx) * data_size_);
		std::memset(g_les, 0, sizeof(cplx) * data_size_);
		for(int w = 0; w < ngrid_; w++)
		{
			Retarded.emplace_back(Eigen::Map<Eigen::Matrix<cplx, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > (p_ret(w,0,0), size1_, size2_));
			Lesser.emplace_back(Eigen::Map<Eigen::Matrix<cplx, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > (p_les(w,0,0), size1_, size2_));
		}


	}

/** \brief <b> Standard constructor for a `GF` object (square matrix) </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > > Constructs a `GF` with all grid (spacing `dgrid`, length `ngrid`), shape (orbital size `size x size`) and type (fermion/boson `sign`, time/freqeuncy `gf_type`) parameters specified manually.
 * > Here assuming a square orbital matrix form.
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param dgrid
 * [double] Grid spacing
 * @param ngrid
 * [int] Number of grid points
 * @param size1
 * [int] Orbital dimension, `GF` is a square matrix
 * @param sign
 * [int] Sign of particle type, +1 for boson,-1 for fermion
 * @param gf_type
 * [int] `GF` type, 0: time, 1: frequency
 */
	GF::GF(double dgrid, int ngrid, int size1, int sign, int gf_type):
		size1_(size1), size2_(size1), sign_(sign),
		dgrid_(dgrid), ngrid_(ngrid), gf_type_(gf_type)
	{
		Retarded.reserve(ngrid_);
		Lesser.reserve(ngrid_);
		int halfn = ngrid / 2 + 1;
		el_size_ = size1_ * size2_;
		//same way of constructing grid as above
		grid_ = new double[ngrid_];
		reset_grid(dgrid);
		data_size_ = ngrid_ * el_size_;
		g_ret = new cplx[data_size_];
		g_les = new cplx[data_size_];
		std::memset(g_ret, 0, sizeof(cplx) * data_size_);
		std::memset(g_les, 0, sizeof(cplx) * data_size_);
		for(int w = 0; w < ngrid_; w++)
		{
			Retarded.emplace_back(Eigen::Map<Eigen::Matrix<cplx, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > (p_ret(w,0,0), size1_, size2_));
			Lesser.emplace_back(Eigen::Map<Eigen::Matrix<cplx, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > (p_les(w,0,0), size1_, size2_));

		}

	}

/** \brief <b> Get the greater component of the Green's function </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Utility to extract the greater component of the `GF` Green's function for a given time / frequency \f$ w \f$.
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param w
 * [int] time / frequency point
 */
	cdmatrix GF::Greater(int w) const
	{
		cdmatrix gtr(size1_, size2_);
		assert(size1_ == size2_);
		assert(w < ngrid_);
		if(gf_type_ == freq_gf)
		{
			for(int idx1 = 0; idx1 < size1_; idx1++)
			{
				for(int idx2 = 0; idx2 < size2_; idx2++)
				{
					gtr(idx1, idx2) = g_ret[w * el_size_ + idx1 * size2_ + idx2].real() 
						- g_ret[w * el_size_ + idx2 * size2_ + idx1].real() +
						NESS_II * (g_ret[w * el_size_ + idx1 * size2_ + idx2].imag() 
								+ g_ret[w * el_size_ + idx2 * size2_ + idx1].imag()) 
						+ g_les[w * el_size_ + idx1 * size2_ + idx2];
				}
			}
		}
		else if(gf_type_ == time_gf)
		{
			if(w > 0)
			{
				for(int idx1 = 0; idx1 < size1_; idx1++)
				{
					for(int idx2 = 0; idx2 < size2_; idx2++)
					{
						gtr(idx1, idx2) = g_ret[w * el_size_ + idx1 * size2_ + idx2] + g_les[w * el_size_ + idx1 * size2_ + idx2];
					}
				}

			}
			else if(w == 0)
			{
				for(int idx1 = 0; idx1 < size1_; idx1++)
				{
					for(int idx2 = 0; idx2 < size2_; idx2++)
					{
						gtr(idx1, idx2) = 2.0 * g_ret[w * el_size_ + idx1 * size2_ + idx2] + g_les[w * el_size_ + idx1 * size2_ + idx2];
					}
				}
			}
		}

		return gtr;
	}

/** \brief <b> Reset the grid of the Green's function </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > Utility to reset the time / frequency grid spacing `dgrid` of the `GF` Green's function to a given value.
* <!-- ARGUMENTS
*      ========= -->
*
* @param dgrid
* [double] New time / frequency grid spacing
*/
	int GF::reset_grid(double dgrid)
	{
		dgrid_ = dgrid;
		int halfn = ngrid_ / 2 + 1;
		if(gf_type_ == freq_gf)
		{
			for (int w_ind = 0; w_ind < halfn; w_ind++)
			{
				grid_[w_ind] = w_ind * dgrid_;
			}
			for (int w_ind = halfn; w_ind < ngrid_; w_ind++)
			{
				grid_[ngrid_ - 1 - (w_ind - halfn)] =  -(w_ind - halfn + 1) * dgrid_;
			}

		}
		else if(gf_type_ == time_gf)
		{
			for (int w_ind = 0; w_ind < ngrid_; w_ind++)
			{
				grid_[w_ind] = w_ind * dgrid_;
			}
		}

		return 0;
	}

};//end of namespace
#endif
