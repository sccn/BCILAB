/*! \file mirdwt.c
    \brief MATLAB gateway for the inverse redundant discrete wavelet transform

    This file is used to produce a MATLAB MEX binary for the inverse redundant discrete wavelet transform

%function x = mirdwt(y_low,y_high,h,L);
% 
% function computes the inverse redundant discrete wavelet transform y for a
% 1D or  2D input signal. redundant means here that the subsampling after
% each stage of the forward transform has been omitted. y_low contains the
% lowpass and y_low the highpass components as computed, e.g., by mrdwt. In
% case of a 2D signal the ordering in y_high is [ncoeff hl hh ncoeff hl ... ] (first
% letter refers to row, second to column filtering).  
%
%    Input:
%       y_low   : lowpass component
%       y_high   : highpass components
%       h    : scaling filter
%       L    : number of levels. in case of a 1D signal length(y_low) must be
%              divisible by 2^L; in case of a 2D signal the row and the
%              column dimension must be divisible by 2^L.
%   
%    Output:
%	x    : finite length 1D or 2D signal
%
% see also: mdwt, midwt, mrdwt
*/

#include "mex.h"
#include "rwt_init.h"
#include "rwt_transforms.h"

/*!
 * Matlab MEX definition for the redundant discrete wavelet transform.
 *
 * @param nlhs number of items on left hand side of matlab call
 * @param plhs pointer to left hand side data structure
 * @param nrhs number of items on right hand side of matlab call
 * @param prhs pointer to right hand side data structure
 *
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  double *x, *yl, *yh;
  rwt_init_params params = rwt_matlab_init(nlhs, plhs, nrhs, prhs, INVERSE_REDUNDANT_DWT);
  yl = mxGetPr(prhs[0]);
  yh = mxGetPr(prhs[1]);
  x = mxGetPr(plhs[0]);
  plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
  *mxGetPr(plhs[1]) = params.levels;
  irdwt(x, params.nrows, params.ncols, params.scalings, params.ncoeff, params.levels, yl, yh);
}

