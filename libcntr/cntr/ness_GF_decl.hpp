#ifndef ness_GF
#define ness_GF

#include "ness_global_settings.hpp"


namespace ness{
enum {time_gf, freq_gf};

/*------ Green function class-------------*/

/** \brief <b> Class `GF` represents Green's functions in either time or frequency domain. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > This class contains the data structures for representing complex matrix valued Green's functions
 * > on a time or frequency grid, stored in two data arrays for retarded and lesser component.
 *
 */
class GF {
	public:
	~GF();
	GF();
	GF(double , int , int , int , int, int gf_type);
	GF(double , int , int , int, int gf_type = freq_gf);
	GF(const GF& g);
#if __cplusplus > 199711L
	GF(GF &&gf);
#endif
    /** \brief <b> Operator = for Green's functions </b> */
	GF& operator=(const GF &g);
    /** \brief <b> Operator + for Green's functions </b> */
	GF operator+(const GF &g);
	void set_shiftGF(int shift);
	int incr(GF &gf, cplx weight=1.0);
	int smul(cplx weight=1.0);
	int clear();

	int left_multiply(const cdmatrix & mat);

	int right_multiply(const cdmatrix & mat);
	int set_element(int, int, const GF&);
	int set_element(int, int, const GF&, int, int);
	int get_element(int, int, GF&) const;
    int get_element(int, int, const GF&, int, int);
	int set_submatrix(const std::vector<int> &, const std::vector<int> &, const GF&);
	int get_submatrix(const std::vector<int> &, const std::vector<int> &, const GF&);

	friend int dyson(GF &invg0, GF &self, GF &g);
	friend int dyson_from_inv(GF &g, GF &self);
    /** \brief <b> Return frequency, not stored in adcending order </b> */
	int gf_idx(int w) { return (gf_type_ == freq_gf ? w - ngrid_ / 2 : w); };
	int reset_grid(double dgrid);

	cdmatrix Greater(int) const;

    /** \brief <b> Grid spacing </b> */
	double dgrid_;
    /** \brief <b> Number of grid points </b> */
	long ngrid_;
    /** \brief <b> Orbital dimension 1 </b> */
	int size1_, size2_; /**< \brief <b> Orbital dimension 2 </b> */
    /// @private
	int el_size_, data_size_;
    /** \brief <b> Green's function type, 0: time, 1: frequency </b> */
	int gf_type_;
    /** \brief <b> Sign +1 for bosons, -1 for fermions. </b> */
	int sign_;
    /// @private
    /** \brief <b> Return pointer to specific data in pointer for retarded component </b> */
	cplx *p_ret(int w_idx, int i, int j)const { return g_ret + w_idx * el_size_ + i * size2_ + j; }
    /// @private
    /** \brief <b> Return pointer to specific data in pointer for lesser component </b> */
	cplx *p_les(int w_idx, int i, int j)const { return g_les + w_idx * el_size_ + i * size2_ + j; }

    /// @private
    /** \brief <b> Raw data pointer for grid </b> */
	double* grid_;
    /// @private
    /** \brief <b> Raw data pointers for Green's function </b> */
	cplx *g_ret, *g_les;
    /** \brief <b> Save raw data pointer in a nice eigen3 structure </b> */
	std::vector<Eigen::Map<Eigen::Matrix<cplx, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>> Retarded, Lesser; /**< \brief <b> Save raw data pointer in a nice eigen3 structure </b> */

};
/** \brief <b> Define a vector of Green's functions </b> */
typedef std::vector<GF> GFs;
/** \brief <b> Define a vector of references of Green's functions </b> */
typedef std::vector<GF & > GF_refs;

};// end of the namespace
#endif

