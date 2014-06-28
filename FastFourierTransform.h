#ifndef FastFourierTransform_h
#define FastFourierTransform_h

/*
 * That code provides fast fourier transformation functionality for series those length N is equal 4*2^n*3^m (n, m >= 0)
 * and a set of functions commonly using in the fourier analysis.
 * 
 * General logic & common functions are placed in the GeneralFourierTransformer class
 * Optimized fast fourier transformation is moved into separate FastFourierTransformer class
 *
 * CalculateFourierCoefficients and BuildCrossFourierCoefficients return fourier coefficients in the internal format:
 * FCos[0], FCos[1], ..., FCos[N/2], FSin[1], ..., FSin[N/2-1]
 * That internal format can be used directly or
 * expand into two sets of FCos[0], ..., FCos[N-1] and FSin[0], ..., FSin[N-1] using BuildFullCosSinCoeffs
 *
 * If you are looking for tutorial or deep knowladge of fourier transformation technique, I recommend
 * "Spectral Analysis and its applications" by G.Jenkins and D.Watts
 */

namespace NumericUtils
{
class GeneralFourierTransformer
{
protected:
	int dim;

	double *p_memory_pool, *p_free_memory;
	double *get_mem (int elements) {double *temp=p_free_memory; p_free_memory+=elements; return temp;};
	void ret_mem (int elements) {p_free_memory -= elements;};
public:
	GeneralFourierTransformer (int length, int extra_pool_rows=0);
	virtual ~GeneralFourierTransformer ();

	virtual bool IsValidObject () const {return dim > 0;};

	int GetSize() const {return dim;};
	int GetAmplSize () const {return (dim > 0) ? (dim/2+1) : 0;};

	virtual bool IsAllowedLength (int length) const {return (length > 0) && !(length%2);};
	virtual bool CalculateFourierCoefficients (const double *p_data, double *p_FX_coeffs) = 0;

	bool BuildFullCosSinCoeffs (const double *p_FX_coeffs, double *p_cos_coeffs, double *p_sin_coeffs) const;

	bool ReverseFourierTransformation (const double *p_FX_coeffs, double *p_data);
	bool ReverseFourierTransformation (const double *p_cos_coeffs, const double *p_sin_coeffs, double *p_data);

	bool BuildAmplitudeSpectrum (const double *p_FX_coeffs, double *p_ampl, double *p_phase) const;
	bool BuildCrossFourierCoefficients (const double *p_FX_coeffs, const double *p_FY_coeffs, double *p_XY_coeffs) const;
};

class FastFourierTransformer : public GeneralFourierTransformer
{
private:
	const double *p_sintbl, *end_sin_tbl;

	int  n1,m;
	double *p_xc_1, *p_xs_1, *p_xc_2, *p_xs_2, *p_xc_3, *p_xs_3;
	double *p_yc, *p_ys, *p_zc, *p_zs, *p_uc, *p_us;
	double val_yc, val_ys, val_zc, val_zs, val_uc, val_us;
	const double *pc_01, *ps_01, *pc_02, *ps_02, *pc_21, *ps_21, *pc_12, *ps_12, *pc_11, *ps_11, *pc_22, *ps_22;
	double vc_01, vs_01, vc_02, vs_02, vc_21, vs_21, vc_12, vs_12, vc_11, vs_11, vc_22, vs_22;

	void recursive_fast_fourier (const double *data, double *f_x, int n, int step);
public:
	FastFourierTransformer (int length);
	virtual ~FastFourierTransformer () {};

	virtual bool IsAllowedLength (int length) const;
	virtual bool CalculateFourierCoefficients (const double *p_data, double *p_FX_coeffs);
};

};

#endif	// FastFourierTransform_h
