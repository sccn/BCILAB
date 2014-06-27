/*! \file mdwt.c
    \brief MATLAB gateway for the discrete wavelet transform

    This file is used to produce a MATLAB MEX binary for the discrete wavelet transform

%y = mdwt(x,h,L);
% 
% function computes the discrete wavelet transform y for a 1D or 2D input
% signal x.
%
%    Input:
%	x    : finite length 1D or 2D signal (implicitely periodized)
%       h    : scaling filter
%       L    : number of levels. in case of a 1D signal length(x) must be
%              divisible by 2^L; in case of a 2D signal the row and the
%              column dimension must be divisible by 2^L.
%
% see also: midwt, mrdwt, mirdwt
*/

#include "mex.h"
#include "rwt_init.h"
#include "rwt_transforms.h"

/*!
 * Matlab MEX definition for the discrete wavelet transform.
 *
 * @param nlhs number of items on left hand side of matlab call
 * @param plhs pointer to left hand side data structure
 * @param nrhs number of items on right hand side of matlab call
 * @param prhs pointer to right hand side data structure
 *
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  rwt_init_params params = rwt_matlab_init(nlhs, plhs, nrhs, prhs, NORMAL_DWT);     /*! Check input and determine the parameters for dwt() */
  plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);                                                               /*! Create the output matrix */
  *mxGetPr(plhs[1]) = params.levels;                                                  /*! The second returned item is the number of levels */
  dwt(mxGetPr(prhs[0]), params.nrows, params.ncols, params.scalings, params.ncoeff, params.levels, mxGetPr(plhs[0]));  /*! Perform the DWT */
}

