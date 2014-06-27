/*! \file midwt.c
    \brief MATLAB gateway for the inverse discrete wavelet transform

    This file is used to produce a MATLAB MEX binary for the inverse discrete wavelet transform

%y = midwt(x,h,L);
% 
% function computes the inverse discrete wavelet transform y for a 1D or 2D
% input signal x.
%
%    Input:
%	x    : finite length 1D or 2D input signal (implicitely periodized)
%       h    : scaling filter
%       L    : number of levels. in case of a 1D signal length(x) must be
%              divisible by 2^L; in case of a 2D signal the row and the
%              column dimension must be divisible by 2^L.
%
% see also: mdwt, mrdwt, mirdwt
*/

#include "mex.h"
#include "rwt_init.h"
#include "rwt_transforms.h"

/*!
 * Matlab MEX definition for the inverse discrete wavelet transform.
 *
 * @param nlhs number of items on left hand side of matlab call
 * @param plhs pointer to left hand side data structure
 * @param nrhs number of items on right hand side of matlab call
 * @param prhs pointer to right hand side data structure
 *
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  double *x, *y;
  rwt_init_params params = rwt_matlab_init(nlhs, plhs, nrhs, prhs, INVERSE_DWT);
  y = mxGetPr(prhs[0]);
  x = mxGetPr(plhs[0]);
  plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
  *mxGetPr(plhs[1]) = params.levels;
  idwt(x, params.nrows, params.ncols, params.scalings, params.ncoeff, params.levels, y);
}

