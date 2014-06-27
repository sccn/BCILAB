/*! \file rwt_init.h
    \brief Header for matlab init functions in init.c
*/
#ifndef RWT_INIT_H_
#define RWT_INIT_H_

#include "rwt_platform.h"

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
  #include "mex.h"
  #ifndef OCTAVE_MEX_FILE
    #include "matrix.h"
  #endif
  typedef struct {
    size_t nrows;     /*!< The number of rows in the input matrix. Output matrix will match.  */
    size_t ncols;     /*!< The number of columns in the input matrix. Output matrix will match. */
    int levels;       /*!< L, the number of levels for the transform. */
    int ncoeff;       /*!< Length of h / the number of scaling coefficients */
    double *scalings; /*!< Wavelet scaling coefficients */
  } rwt_init_params;
  typedef enum {NORMAL_DWT, REDUNDANT_DWT, INVERSE_DWT, INVERSE_REDUNDANT_DWT} transform_t;
#endif

#ifdef __cplusplus
extern "C" {
#endif

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
  rwt_init_params rwt_matlab_init(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[], transform_t dwtType);
#else
  int rwt_find_levels(size_t m, size_t n);
  int rwt_check_levels(int levels, size_t rows, size_t cols);
#endif

#ifdef __cplusplus
}
#endif

#endif /* RWT_INIT_H_ */
