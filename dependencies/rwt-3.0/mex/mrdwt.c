/*! \file mrdwt.c
    \brief MATLAB gateway for the redundant discrete wavelet transform

    This file is used to produce a MATLAB MEX binary for the redundant discrete wavelet transform

%[yl,yh] = mrdwt(x,h,L);
% 
% function computes the redundant discrete wavelet transform y for a 1D or
% 2D input signal . redundant means here that the subsampling after each
% stage is omitted. yl contains the lowpass and yl the highpass
% components. In case of a 2D signal the ordering in yh is [ncoeff hl hh ncoeff hl
% ... ] (first letter refers to row, second to column filtering). 
%
%    Input:
%	x    : finite length 1D or 2D signal (implicitely periodized)
%       h    : scaling filter
%       L    : number of levels. in case of a 1D signal length(x) must be
%              divisible by 2^L; in case of a 2D signal the row and the
%              column dimension must be divisible by 2^L.
%   
%    Output:
%       yl   : lowpass component
%       yh   : highpass components
%
% see also: mdwt, midwt, mirdwt
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
  rwt_init_params params = rwt_matlab_init(nlhs, plhs, nrhs, prhs, REDUNDANT_DWT);
  if (min(params.nrows, params.ncols) == 1)
    plhs[1] = mxCreateDoubleMatrix(params.nrows, params.levels*params.ncols, mxREAL);
  else
    plhs[1] = mxCreateDoubleMatrix(params.nrows, 3*params.levels*params.ncols, mxREAL);
  x = mxGetPr(prhs[0]);
  yl = mxGetPr(plhs[0]);
  yh = mxGetPr(plhs[1]);
  plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
  *mxGetPr(plhs[2]) = params.levels;
  rdwt(x, params.nrows, params.ncols, params.scalings, params.ncoeff, params.levels, yl, yh);
}

