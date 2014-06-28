#include <stddef.h>		// NULL definition
#include <math.h>

#include "FastFourierTransform.h"

using namespace NumericUtils;

NumericUtils::GeneralFourierTransformer::GeneralFourierTransformer (int length, int extra_pool_rows) :
						dim (IsAllowedLength(length) ? length : 0),
						p_memory_pool (NULL),
						p_free_memory (NULL)
{
	if (dim > 0)
	{
		if (extra_pool_rows < 0)
			extra_pool_rows = 0;

		p_memory_pool = new double [(extra_pool_rows+6)*dim];
		p_free_memory = p_memory_pool;
	}
}

NumericUtils::GeneralFourierTransformer::~GeneralFourierTransformer ()
{
	if (p_memory_pool)
	{
		delete [] p_memory_pool;
	}
}

bool NumericUtils::GeneralFourierTransformer::BuildFullCosSinCoeffs (const double *p_FX_coeffs, double *p_cos_coeffs, double *p_sin_coeffs) const
{
	if (!IsValidObject() || !p_FX_coeffs)
		return false;

	if (p_cos_coeffs)
	{
		const double *p_FX_curr=p_FX_coeffs;
		double *p_left = p_cos_coeffs;
		double *p_right = p_cos_coeffs + dim - 1;

		*p_left++ = *p_FX_curr++;

		for (int cnt=1; cnt < dim/2; ++cnt)
		{
			const double curr_val = *p_FX_curr++;

			*p_left++ = curr_val;
			*p_right-- = curr_val;
		}

		*p_left = *p_FX_curr;
	}

	if (p_sin_coeffs)
	{
		const double *p_FX_curr = p_FX_coeffs + dim/2 + 1;
		double *p_left = p_sin_coeffs;
		double *p_right = p_sin_coeffs + dim - 1;

		*p_left++ = .0;

		for (int cnt=1; cnt < dim/2; ++cnt)
		{
			const double curr_val = *p_FX_curr++;

			*p_left++ = curr_val;
			*p_right-- = -curr_val;
		}

		*p_left = .0;
	}

	return true;
}

bool NumericUtils::GeneralFourierTransformer::ReverseFourierTransformation (const double *p_cos_coeffs, const double *p_sin_coeffs, double *p_data)
{
	if (!IsValidObject() || !p_cos_coeffs || !p_sin_coeffs || !p_data)
		return false;

	double *p_U = get_mem (4*dim);
	double *p_V = p_U + dim;
	double *p_uc = p_V +dim;
	double *p_vs = p_uc + dim;

	const bool ret_val = CalculateFourierCoefficients (p_cos_coeffs, p_U) &&
						 CalculateFourierCoefficients (p_sin_coeffs, p_V) &&
						 BuildFullCosSinCoeffs (p_U, p_uc, NULL) &&
						 BuildFullCosSinCoeffs (p_V, NULL, p_vs);

	if (ret_val)
	{
		for (int cnt=0; cnt<dim; ++cnt)
			*p_data++ = dim * (*p_uc++ + *p_vs++);
	}

	ret_mem (4*dim);

	return ret_val;
}

bool NumericUtils::GeneralFourierTransformer::ReverseFourierTransformation (const double *p_FX_coeffs, double *p_data)
{
	if (!IsValidObject() || !p_FX_coeffs || !p_data)
		return false;

	double *p_XC = get_mem (2*dim);
	double *p_XS = p_XC + dim;

	const bool ret_val = BuildFullCosSinCoeffs (p_FX_coeffs, p_XC, p_XS) &&
						 ReverseFourierTransformation (p_XC, p_XS, p_data);

	ret_mem (2*dim);

	return ret_val;
}

bool NumericUtils::GeneralFourierTransformer::BuildAmplitudeSpectrum (const double *p_FX_coeffs, double *p_ampl, double *p_phase) const
{
	if (!IsValidObject() || !p_FX_coeffs)
		return false;

	const double *p_XC = p_FX_coeffs;
	const double *p_XS = p_FX_coeffs + dim/2 + 1;

	double curr_xc = *p_XC++;
	double curr_xs = .0;

	if (curr_xc < .0)
	{
		if (p_ampl)
			*p_ampl++ = -curr_xc;
		if (p_phase)
			*p_phase++ = M_PI;
	}
	else
	{
		if (p_ampl)
			*p_ampl++ = curr_xc;
		if (p_phase)
			*p_phase++ = .0;
	}

	for (int cnt=1; cnt<dim/2; ++cnt)
	{
		curr_xc = *p_XC++;
		curr_xs = *p_XS++;

		if (p_ampl)
			*p_ampl++ = M_SQRT2 * sqrt (curr_xc*curr_xc + curr_xs*curr_xs);
		if (p_phase)
			*p_phase++ = atan2 (curr_xc, curr_xs);
	}

	curr_xc = *p_XC++;

	if (curr_xc < .0)
	{
		if (p_ampl)
			*p_ampl++ = -curr_xc;
		if (p_phase)
			*p_phase++ = M_PI;
	}
	else
	{
		if (p_ampl)
			*p_ampl++ = curr_xc;
		if (p_phase)
			*p_phase++ = .0;
	}

	return true;
}

bool NumericUtils::GeneralFourierTransformer::BuildCrossFourierCoefficients (const double *p_FX_coeffs, const double *p_FY_coeffs, double *p_XY_coeffs) const
{
	if (!IsValidObject() || !p_FX_coeffs || !p_FY_coeffs || !p_XY_coeffs)
		return false;

	const double *p_XC = p_FX_coeffs;
	const double *p_XS = p_FX_coeffs + dim/2 + 1;
	const double *p_YC = p_FY_coeffs;
	const double *p_YS = p_FY_coeffs + dim/2 + 1;
	double *p_XYC = p_XY_coeffs;
	double *p_XYS = p_XY_coeffs + dim/2 + 1;

	double curr_xc = *p_XC++;
	double curr_xs = .0;
	double curr_yc = *p_YC++;
	double curr_ys = .0;

	*p_XYC++ = curr_xc * curr_yc;

	for (int cnt=1; cnt < dim/2 + 1; ++cnt)
	{
		curr_xc = *p_XC++;
		curr_xs = *p_XS++;
		curr_yc = *p_YC++;
		curr_ys = *p_YS++;

		*p_XYC++ = curr_xc * curr_yc + curr_xs * curr_ys;
		*p_XYS++ = curr_xc * curr_ys - curr_xs * curr_yc;
	}

	*p_XYC = *p_XC * *p_YC;

	return true;
}

NumericUtils::FastFourierTransformer::FastFourierTransformer (int length) : GeneralFourierTransformer (length, 3),
																			p_sintbl (NULL),
																			end_sin_tbl (NULL)
{
	if (!IsAllowedLength (length))
		dim = 0;

	if (dim > 0)
	{
		double *p_init_sintbl = get_mem (dim);
		p_init_sintbl[0] = .0;
		p_init_sintbl[dim/4] = 1.;
		p_init_sintbl[dim/2] = .0;
		p_init_sintbl[3*dim/4] = -1.;

		for (int cnt=1; cnt<dim/4; ++cnt)
		{
			const double temp = sin (2.*M_PI*cnt/dim);

			p_init_sintbl[cnt] = temp;
			p_init_sintbl[dim/2-cnt] = temp;
			p_init_sintbl[dim/2+cnt] = -temp;
			p_init_sintbl[dim-cnt] = -temp;
		}

		p_sintbl = p_init_sintbl;
		end_sin_tbl = p_sintbl + dim - 1;
	}
}

bool NumericUtils::FastFourierTransformer::IsAllowedLength (int length) const
{
	if ((length < 4) || (length % 4))
		return false;

	length /= 4;

	while (length > 4)
	{
		if (!(length % 3))
			length /= 3;
		else if (!(length % 2))
			length /= 2;
		else
			return false;
	}

	return true;
}

bool NumericUtils::FastFourierTransformer::CalculateFourierCoefficients (const double *p_data, double *p_FX_coeffs)
{
	if (!IsValidObject() || !p_data || !p_FX_coeffs)
		return false;

	recursive_fast_fourier (p_data, p_FX_coeffs, dim, 1);

	return true;
}
static const double val_sv=.288675134594812866;  // val_sv = sqrt(3.)/6

void NumericUtils::FastFourierTransformer::recursive_fast_fourier (const double *data, double *f_x, int n, int step)
{
  if (n==2)
  {
    f_x[0] =(data[0] + data[step]) * 0.5;  // FXc[0]=(x[0]+x[1])/2
    f_x[1] =(data[0] - data[step]) * 0.5;  // FXc[1]=(x[0]-x[1])/2
    return;
  }

  double *f_y=get_mem (n);

  if (!(n%3))
  {
    recursive_fast_fourier (data, f_y, n/3, 3*step);
    recursive_fast_fourier (data+step, f_y+n/3, n/3, 3*step);
    recursive_fast_fourier (data+2*step, f_y+2*n/3, n/3, 3*step);

/*
 *        p_xc_1 points to FXc[m]		p_xs_1 points to FXs[m]
 *        p_xc_2 points to FXc[n/3-m]	p_xs_2 points to FXs[n/3-m]
 *        p_xc_3 points to FXc[n/3+m]	p_xs_3 points to FXs[n/3+m]
 */
    n1 = n/6; m = 2*n1;
    p_xc_1 = f_x; p_xc_2 = p_xc_1 + m; p_xc_3 = p_xc_2 + 1;
    p_xs_1 = p_xc_3 + n1; p_xs_2 = p_xs_1 + m - 1; p_xs_3 = p_xs_2 +1;

/*
 *        FXc[0]	= (FYc[0] + FZc[0] + FUc[0])/3.
 *        FXc[n/3]	= FYc[0]/3. - FZc[0]/6. - FUc[0]/6.
 *        FXs[n/3]	= (FZc[0] - FUc[0]) * sqrt(3.)/6.
 */
    p_yc = f_y; p_zc = p_yc + m; p_uc = p_zc + m;
    val_yc = *p_yc++; p_ys = p_yc + n1;
    val_zc = *p_zc++; p_zs = p_zc + n1;
    val_uc = *p_uc++; p_us = p_uc + n1;
    
    *p_xc_1++ = (val_yc + val_zc + val_uc) /3.;
    *p_xc_2-- = val_yc/3. - val_zc/6. - val_uc/6.;
    *p_xs_2-- = (val_zc - val_uc) *val_sv;

/*
 *        ps_ij points to sin(2*PI/n*(i*n/3 + j*m))
 *        pc_ij points to cos(2*PI/n*(i*n/3 + j*m))
 */
    n1 = dim/3; m = dim/4;
    ps_01 = p_sintbl + step; pc_01 = ps_01 + m;
    ps_11 = ps_01 + n1; pc_11 = ps_11 + m;
    ps_21 = ps_11 + n1; pc_21 = ps_21 + m;
    if (pc_21 > end_sin_tbl) pc_21 -= dim;
    ps_02 = ps_01 + step; pc_02 = ps_02 + m;
    ps_12 = ps_02 + n1; pc_12 = ps_12 + m;
    if (pc_12 > end_sin_tbl) pc_12 -= dim;
    ps_22 = ps_12 + n1; pc_22 = ps_22 + m;
    if (ps_22 > end_sin_tbl) ps_22 -= dim;
    if (pc_22 > end_sin_tbl) pc_22 -= dim;

    n1 = n/6;
    for (m=1; m<n1; ++m)
    {
      val_yc = *p_yc++; val_ys = *p_ys++;
      val_zc = *p_zc++; val_zs = *p_zs++;
      val_uc = *p_uc++; val_us = *p_us++;
      vc_01 = *pc_01; pc_01 += step; if (pc_01 > end_sin_tbl) pc_01 -= dim;
      vs_01 = *ps_01; ps_01 += step; if (ps_01 > end_sin_tbl) ps_01 -= dim;
      vc_11 = *pc_11; pc_11 += step; if (pc_11 > end_sin_tbl) pc_11 -= dim;
      vs_11 = *ps_11; ps_11 += step; if (ps_11 > end_sin_tbl) ps_11 -= dim;
      vc_21 = *pc_21; pc_21 += step; if (pc_21 > end_sin_tbl) pc_21 -= dim;
      vs_21 = *ps_21; ps_21 += step; if (ps_21 > end_sin_tbl) ps_21 -= dim;
      vc_02 = *pc_02; pc_02 += 2*step; if (pc_02 > end_sin_tbl) pc_02 -= dim;
      vs_02 = *ps_02; ps_02 += 2*step; if (ps_02 > end_sin_tbl) ps_02 -= dim;
      vc_12 = *pc_12; pc_12 += 2*step; if (pc_12 > end_sin_tbl) pc_12 -= dim;
      vs_12 = *ps_12; ps_12 += 2*step; if (ps_12 > end_sin_tbl) ps_12 -= dim;
      vc_22 = *pc_22; pc_22 += 2*step; if (pc_22 > end_sin_tbl) pc_22 -= dim;
      vs_22 = *ps_22; ps_22 += 2*step; if (ps_22 > end_sin_tbl) ps_22 -= dim;
/*
 *        3*FXc[m]		= FYc[m] + FZc[m]*vc_01 - FZs[m]*vs_01 + FUc[m]*vc_02 - FUs[m]*vs_02
 *        3*FXs[m]		= FYs[m] + FZs[m]*vc_01 + FZc[m]*vs_01 + FUs[m]*vc_02 + FUc[m]*vs_02
 *        3*FXc[n/3-m]	= FYc[m] + FZc[m]*vc_21 - FZs[m]*vs_21 + FUc[m]*vc_12 - FUs[m]*vs_12
 *        3*FXs[n/3-m]	=-FYs[m] - FZs[m]*vc_21 - FZc[m]*vs_21 - FUs[m]*vc_12 - FUc[m]*vs_12
 *        3*FXc[n/3+m]	= FYc[m] + FZc[m]*vc_11 - FZs[m]*vs_11 + FUc[m]*vc_22 - FUs[m]*vs_22
 *        3*FXs[n/3+m]	= FYs[m] + FZs[m]*vc_11 + FZc[m]*vs_11 + FUs[m]*vc_22 + FUc[m]*vs_22
 */      
      *p_xc_1++ = (val_yc + val_zc*vc_01 - val_zs*vs_01 + val_uc*vc_02 - val_us*vs_02) / 3.;
      *p_xs_1++ = (val_ys + val_zs*vc_01 + val_zc*vs_01 + val_us*vc_02 + val_uc*vs_02) / 3.;
      *p_xc_2-- = (val_yc + val_zc*vc_21 - val_zs*vs_21 + val_uc*vc_12 - val_us*vs_12) / 3.;
      *p_xs_2-- =(-val_ys - val_zs*vc_21 - val_zc*vs_21 - val_us*vc_12 - val_uc*vs_12) / 3.;
      *p_xc_3++ = (val_yc + val_zc*vc_11 - val_zs*vs_11 + val_uc*vc_22 - val_us*vs_22) / 3.;
      *p_xs_3++ = (val_ys + val_zs*vc_11 + val_zc*vs_11 + val_us*vc_22 + val_uc*vs_22) / 3.;
    }
/*
 *        FXc[n/6] = FYc[n/6]/3. + FZc[n/6]/6. - FUc[n/6]/6.
 *        FXs[n/6] = (FZc[n/6] + FUc[n/6]) * sqrt(3.)/6.
 *        FXc[n/2] = (FYc[n/6] - FZc[n/6] + FUc[n/6])/3.
 */
    val_yc = *p_yc; val_zc = *p_zc; val_uc = *p_uc;
    *p_xc_1 = val_yc/3. + val_zc/6. - val_uc/6.;
    *p_xc_3 = (val_yc - val_zc + val_uc) / 3.;
    *p_xs_1 = (val_zc + val_uc)*val_sv;
  }
  else      // n=2*k
  {
    recursive_fast_fourier (data, f_y, n/2, 2*step);
    recursive_fast_fourier (data+step, f_y+n/2, n/2, 2*step);

/*
 *        p_xc_1 points to FXc[m]		p_xs_1 points to FXs[m]
 *        p_xc_2 points to FXc[n/2-m]	p_xs_2 points to FXs[n/2-m]
 *
 *        FXc[0] = (FYc[0] + FZc[0])/2.
 *        FXc[n/2] = (FYc[0] - FZc[0])/2.
 */
    n1 = n/4; m = 2*n1;
    p_xc_1 = f_x; p_xc_2 = p_xc_1 + m;
    p_yc = f_y; p_zc = p_yc + m;
    val_yc = *p_yc++; p_ys = p_yc + n1;
    val_zc = *p_zc++; p_zs = p_zc + n1;

    *p_xc_1++ = (val_yc + val_zc)*.5; p_xs_1 = p_xc_1 + m;
    *p_xc_2-- = (val_yc - val_zc)*.5; p_xs_2 = p_xc_2 + m;

    ps_01 = p_sintbl + step; pc_01 = ps_01 + dim/4;
    
    for (m=1; m<n1; ++m)
    {
      val_yc = *p_yc++; val_ys = *p_ys++;
      val_zc = *p_zc++; val_zs = *p_zs++;
      vc_01 = *pc_01; pc_01 += step; if (pc_01 > end_sin_tbl) pc_01 -= dim;
      vs_01 = *ps_01; ps_01 += step; if (ps_01 > end_sin_tbl) ps_01 -= dim;

/*
 *        2*FXc[m]     = FYc[m] + FZc[m]*vc_01 - FZs[m]*vs_01
 *        2*FXs[m]     = FYs[m] + FZs[m]*vc_01 + FZc[m]*vs_01    
 *        2*FXc[n/2-m] = FYc[m] - FZc[m]*vc_01 + FZs[m]*vs_01
 *        2*FXs[n/2-m] =-FYs[m] + FZs[m]*vc_01 + FZc[m]*vs_01
 */
      *p_xc_1++ = (val_yc + val_zc*vc_01 - val_zs*vs_01)*.5;
      *p_xs_1++ = (val_ys + val_zs*vc_01 + val_zc*vs_01)*.5;
      *p_xc_2-- = (val_yc - val_zc*vc_01 + val_zs*vs_01)*.5;
      *p_xs_2-- = (-val_ys + val_zs*vc_01 + val_zc*vs_01)*.5;
    }
    *p_xc_1 = *p_yc*.5;    // FXc[n/4] = FYc[n/4]/2.
    *p_xs_1 = *p_zc*.5;    // FXs[n/4] = FZc[n/4]/2.
  }

  ret_mem (n);

  return;
}
