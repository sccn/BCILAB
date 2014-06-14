/*! \file init.c
    \brief Parse input from MATLAB and do some sanity checking

*/

#include "rwt_init.h"
#include <math.h>

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
/*!
 * Checks for correct # of input variables based on type of transform.
 *
 * @param nrhs number of items on right hand side of matlab call
 * @param transform_type 
 *
 */
int rwt_check_parameter_count(int nrhs, transform_t transform_type) {
  if (transform_type == INVERSE_REDUNDANT_DWT) {
    if (nrhs > 4) {
      rwt_errormsg("There are at most 4 input parameters allowed!");
      return 1;
    }
    if (nrhs < 3) {
      rwt_errormsg("There are at least 3 input parameters required!");
      return 1;
    }
  }
  else {
    if (nrhs > 3) {
      rwt_errormsg("There are at most 3 input parameters allowed!");
      return 1;
    }
    if (nrhs < 2) {
      rwt_errormsg("There are at least 2 input parameters required!");
      return 1;
    }
  }
  return 0;
}


/*!
 * For the inverse redundant transform check that the dimensions of the low and high inputs match
 *
 * @param prhs
 * @param params
 *
 */
int rwt_check_yl_matches_yh(const mxArray *prhs[], size_t nrows, size_t ncols, int levels) {
  size_t mh = mxGetM(prhs[1]);
  size_t nh = mxGetN(prhs[1]);
  if (min(nrows, ncols) > 1) {
    if ((nrows != mh) | (3 * ncols * levels != nh)) {
      return 0;
    }
  }
  else {
    if ((nrows != mh) | (ncols * levels != nh)) {
      return 0;
    }
  }
  return 1;
}
#endif


/*!
 * Find L, the number of levels
 *
 * @param m the number of rows in the input
 * @param n the number of columns in the input
 *
 * L is the exponent of the largest power of 2 that is a factor of all input dimensions
 * 
 */
int rwt_find_levels(size_t m, size_t n) {
  size_t i, j, L;
  i = n ; j = 0;
  while (even(i)) {
    i = (i >> 1);
    j++;
  }
  L = m; i = 0;
  while (even(L)) {
    L = (L >> 1);
    i++;
  }
  if (min(m, n) == 1)
    L = max(i, j);
  else
    L = min(i, j);
  if (L == 0) {
    rwt_errormsg("Maximum number of levels is zero; no decomposition can be performed!");
    return -1;
  }
  else return L;
}


/*!
 * Check that length is divisble by 2^L
 *
 * @param length the number of rows or number of columns
 * @param L the number of levels
 *
 */
int rwt_check_dimensions(size_t length, int L) {
  double test = (double) length / pow(2.0, (double) L);
  if ((test - floor(test)) > 0.0) {
    return -1;
  }
  return 0;
}


/*!
 * Sanity check the levels parameter
 *
 * @param levels the number of levels specified or calculated for the input
 * @param rows the number of rows of input
 * @param cols the number of columns of input
 *
 */
int rwt_check_levels(int levels, size_t rows, size_t cols) {
  if (levels < 1) {
    rwt_errormsg("The number of levels, L, must be a positive integer");
    return -1;
  }

  /*! Check that both the rows and columns are divisible by 2^L */
  if ((rows > 1 && rwt_check_dimensions(rows, levels)) || (cols > 1 && rwt_check_dimensions(cols, levels))) {
    rwt_errormsg("All dimensions must be divisible by 2^L");
    return -1;
  }

  return 0;
}


#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
/*!
 * Parse input from MATLAB and do some sanity checking
 *
 * @param nlhs number of items on left hand side of matlab call
 * @param plhs pointer to left hand side data structure
 * @param nrhs number of items on right hand side of matlab call
 * @param prhs pointer to right hand side data structure
 * @param transform_type which transform are we setting up to do
 *
 */
rwt_init_params rwt_matlab_init(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[], transform_t transform_type) {
  rwt_init_params params;
  int argNumL = (transform_type == INVERSE_REDUNDANT_DWT) ? 3 : 2;
  /*! Check for correct # of input parameters */
  if (rwt_check_parameter_count(nrhs, transform_type) != 0) return params;
  /*! Check that we don't have more than two dimensions in the input since that is currently unsupported. */
  if (mxGetNumberOfDimensions(prhs[0]) > 2) {
    rwt_errormsg("Matrix must have fewer than 3 dimensions!");
    return params;
  }
  /*! Get the number of rows and columns in the input matrix. */
  params.nrows = mxGetM(prhs[0]);
  params.ncols = mxGetN(prhs[0]);

  if (params.nrows == 0 && params.ncols == 0) {
    rwt_errormsg("The input matrix cannot be empty");
    return params;
  }

  /*! Read the number of levels, L, from the input values if it was given, otherwise calculate L. Sanity check L */
  if ((argNumL + 1) == nrhs)
    params.levels = (int) *mxGetPr(prhs[argNumL]);
  else
    params.levels = rwt_find_levels(params.nrows, params.ncols);

  if (rwt_check_levels(params.levels, params.nrows, params.ncols)) {
    return params;
  }

  /*! Read the scaling coefficients, h, from the input and find their length, ncoeff. 
   *  In the case of the redundant transform, the scalings are found one further position to the right, 
   *  and also we check for matching dimensions in the low and high inputs
   */
  if (transform_type == INVERSE_REDUNDANT_DWT) {
    params.scalings = mxGetPr(prhs[2]);
    params.ncoeff = max(mxGetM(prhs[2]), mxGetN(prhs[2]));
    if (!rwt_check_yl_matches_yh(prhs, params.nrows, params.ncols, params.levels)) {
      rwt_errormsg("Dimensions of first two input matrices not consistent!");
      return params;
    }
  }
  else {
    params.scalings = mxGetPr(prhs[1]);
    params.ncoeff = max(mxGetM(prhs[1]), mxGetN(prhs[1]));
  }
  /*! Create the first item in the output array as a double matrix with the same dimensions as the input. */
  plhs[0] = mxCreateDoubleMatrix(params.nrows, params.ncols, mxREAL);
  return params;
}
#endif
