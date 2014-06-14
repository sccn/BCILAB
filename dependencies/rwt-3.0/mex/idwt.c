/*! \file idwt.c
    \brief Implementation of the inverse discrete wavelet transform

*/

#include "rwt_platform.h"

/*!
 * Perform convolution for idwt
 *
 * @param x_out
 * @param lx
 * @param coeff_low
 * @param coeff_high
 * @param ncoeff_minus_one
 * @param ncoeff_halved_minus_one
 * @param x_in_low
 * @param x_in_high
 * 
 */
void idwt_convolution(double *x_out, size_t lx, double *coeff_low, double *coeff_high, int ncoeff_minus_one, int ncoeff_halved_minus_one, double *x_in_low, double *x_in_high) {
  int k;
  size_t i, j, ind, tj;
  double x0, x1;

  for (k=ncoeff_halved_minus_one-1; k > -1; k--) {
    x_in_low[k]  = x_in_low[lx+k];
    x_in_high[k] = x_in_high[lx+k];
  }

  ind = 0;
  for (i=0; i<(lx); i++) {
    x0 = 0;
    x1 = 0;
    tj = 0;
    for (j=0; j<=ncoeff_halved_minus_one; j++) {
      x0 = x0 + (x_in_low[i+j] * coeff_low[ncoeff_minus_one-1-tj]) + (x_in_high[i+j] * coeff_high[ncoeff_minus_one-1-tj]);
      x1 = x1 + (x_in_low[i+j] * coeff_low[ncoeff_minus_one-tj])   + (x_in_high[i+j] * coeff_high[ncoeff_minus_one-tj]);
      tj += 2;
    }
    x_out[ind++] = x0;
    x_out[ind++] = x1;
  }
}


/*!
 * Allocate memory for idwt
 *
 * @param m the number of rows of the input matrix
 * @param n the number of columns of the input matrix
 * @param ncoeff the number of scaling coefficients
 * @param x_dummy
 * @param y_dummy_low
 * @param y_dummy_high
 * @param coeff_low
 * @param coeff_high
 *
 */
void idwt_allocate(size_t m, size_t n, int ncoeff, double **x_dummy, double **y_dummy_low, double **y_dummy_high, double **coeff_low, double **coeff_high) {
  *x_dummy      = (double *) rwt_calloc(max(m,n),            sizeof(double));
  *y_dummy_low  = (double *) rwt_calloc(max(m,n)+ncoeff/2-1, sizeof(double));
  *y_dummy_high = (double *) rwt_calloc(max(m,n)+ncoeff/2-1, sizeof(double));
  *coeff_low    = (double *) rwt_calloc(ncoeff,              sizeof(double));
  *coeff_high   = (double *) rwt_calloc(ncoeff,              sizeof(double));
}


/*!
 * Free memory we allocated for idwt
 *
 * @param x_dummy
 * @param y_dummy_low
 * @param y_dummy_high
 * @param coeff_low
 * @param coeff_high
 *
 */
void idwt_free(double **x_dummy, double **y_dummy_low, double **y_dummy_high, double **coeff_low, double **coeff_high) {
  rwt_free(*x_dummy);
  rwt_free(*y_dummy_low);
  rwt_free(*y_dummy_high);
  rwt_free(*coeff_low);
  rwt_free(*coeff_high);
}


/*!
 * Put the scaling coeffients into a form ready for use in the convolution function
 *
 * @param ncoeff length of h / the number of scaling coefficients
 * @param h  the wavelet scaling coefficients
 * @param coeff_low same as h
 * @param coeff_high reversed h, even values are sign flipped
 *
 */
void idwt_coefficients(int ncoeff, double *h, double **coeff_low, double **coeff_high) {
  int i;
  for (i=0; i<ncoeff; i++) {
    (*coeff_low)[i] = h[i];
    (*coeff_high)[i] = h[ncoeff-i-1];
  }
  for (i=1; i<=ncoeff; i+=2)
    (*coeff_high)[i] = -((*coeff_high)[i]);
}


/*!
 * Perform the inverse discrete wavelet transform
 *
 * @param x      the output signal with the inverse wavelet transform applied
 * @param nrows  number of rows in the input
 * @param ncols  number of columns in the input
 * @param h      wavelet scaling coefficients
 * @param ncoeff the number of scaling coefficients
 * @param levels the number of levels
 * @param y      the input signal
 *
 */
void idwt(double *x, size_t nrows, size_t ncols, double *h, int ncoeff, int levels, double *y) {
  double  *coeff_low, *coeff_high, *y_dummy_low, *y_dummy_high, *x_dummy;
  long i;
  int current_level, ncoeff_minus_one, ncoeff_halved_minus_one, sample_f;
  size_t current_rows, current_cols, row_cursor, column_cursor, idx_rows, idx_cols;

  idwt_allocate(nrows, ncols, ncoeff, &x_dummy, &y_dummy_low, &y_dummy_high, &coeff_low, &coeff_high);
  idwt_coefficients(ncoeff, h, &coeff_low, &coeff_high);

  if (ncols==1) {
    ncols = nrows;
    nrows = 1;
  }
  
  ncoeff_minus_one = ncoeff - 1;
  ncoeff_halved_minus_one = ncoeff/2 - 1;
  /* 2^levels */
  sample_f = 1;
  for (i=1; i<levels; i++)
    sample_f = sample_f*2;
  
  if (nrows>1)
    current_rows = nrows/sample_f;
  else 
    current_rows = 1;
  current_cols = ncols/sample_f;

  for (i=0; i<(nrows*ncols); i++)
    x[i] = y[i];
  
  /* main loop */
  for (current_level=levels; current_level >= 1; current_level--) {
    row_cursor = current_rows/2;
    column_cursor = current_cols/2;
    
    /* go by columns in case of a 2D signal*/
    if (nrows>1) {
      for (idx_cols=0; idx_cols<current_cols; idx_cols++) {         /* loop over columns */
	/* store in dummy variables */
	idx_rows = row_cursor;
	for (i=0; i<row_cursor; i++){    
	  y_dummy_low[i+ncoeff_halved_minus_one]  = mat(x, i,          idx_cols, nrows, ncols);  
	  y_dummy_high[i+ncoeff_halved_minus_one] = mat(x, idx_rows++, idx_cols, nrows, ncols);  
	}
	/* perform filtering lowpass and highpass*/
	idwt_convolution(x_dummy, row_cursor, coeff_low, coeff_high, ncoeff_minus_one, ncoeff_halved_minus_one, y_dummy_low, y_dummy_high); 
	/* restore dummy variables in matrix */
	for (i=0; i<current_rows; i++)
	  mat(x, i, idx_cols, nrows, ncols) = x_dummy[i];  
      }
    }
    /* go by rows */
    for (idx_rows=0; idx_rows<current_rows; idx_rows++) {           /* loop over rows */
      /* store in dummy variable */
      idx_cols = column_cursor;
      for  (i=0; i<column_cursor; i++){    
	y_dummy_low[i+ncoeff_halved_minus_one]  = mat(x, idx_rows, i,          nrows, ncols);  
	y_dummy_high[i+ncoeff_halved_minus_one] = mat(x, idx_rows, idx_cols++, nrows, ncols);  
      } 
      /* perform filtering lowpass and highpass*/
      idwt_convolution(x_dummy, column_cursor, coeff_low, coeff_high, ncoeff_minus_one, ncoeff_halved_minus_one, y_dummy_low, y_dummy_high); 
      /* restore dummy variables in matrices */
      for (i=0; i<current_cols; i++)
        mat(x, idx_rows, i, nrows, ncols) = x_dummy[i];  
    }  
    if (nrows==1)
      current_rows = 1;
    else
      current_rows = current_rows*2;
    current_cols = current_cols*2;
  }
  idwt_free(&x_dummy, &y_dummy_low, &y_dummy_high, &coeff_low, &coeff_high);
}

