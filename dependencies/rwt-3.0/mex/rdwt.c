/*! \file rdwt.c
    \brief Implementation of the redundant discrete wavelet transform

*/

#include "rwt_platform.h"

/*!
 * Perform convolution for rdwt
 *
 * @param x_in input signal values
 * @param lx the length of x
 * @param coeff_low the low pass coefficients
 * @param coeff_high the high pass coefficients
 * @param ncoeff the number of scaling coefficients
 * @param x_out_low low pass results
 * @param x_out_high high pass results
 * 
 */
void rdwt_convolution(double *x_in, size_t lx, double *coeff_low, double *coeff_high, int ncoeff, double *x_out_low, double *x_out_high) {
  size_t i, j;
  double x0, x1;

  for (i=lx; i < lx+ncoeff-1; i++)
    x_in[i] = x_in[i-lx];
  for (i=0; i<lx; i++) {
    x0 = 0;
    x1 = 0;
    for (j=0; j<ncoeff; j++) {
      x0 = x0 + x_in[j+i] * coeff_low[ncoeff-1-j];
      x1 = x1 + x_in[j+i] * coeff_high[ncoeff-1-j];
    }
    x_out_low[i] = x0;
    x_out_high[i] = x1;
  }
}


/*!
 * Allocate memory for rdwt
 *
 * @param m the number of rows of the input matrix
 * @param n the number of columns of the input matrix
 * @param ncoeff the number of scaling coefficients
 * @param x_dummy_low
 * @param x_dummy_high
 * @param y_dummy_low_low
 * @param y_dummy_low_high
 * @param y_dummy_high_low
 * @param y_dummy_high_high
 * @param coeff_low
 * @param coeff_high
 *
 */
void rdwt_allocate(size_t m, size_t n, int ncoeff, double **x_dummy_low, double **x_dummy_high, double **y_dummy_low_low, 
  double **y_dummy_low_high, double **y_dummy_high_low, double **y_dummy_high_high, double **coeff_low, double **coeff_high) {
  *x_dummy_low       = (double *) rwt_calloc(max(m,n)+ncoeff-1, sizeof(double));
  *x_dummy_high      = (double *) rwt_calloc(max(m,n)+ncoeff-1, sizeof(double));
  *y_dummy_low_low   = (double *) rwt_calloc(max(m,n),          sizeof(double));
  *y_dummy_low_high  = (double *) rwt_calloc(max(m,n),          sizeof(double));
  *y_dummy_high_low  = (double *) rwt_calloc(max(m,n),          sizeof(double));
  *y_dummy_high_high = (double *) rwt_calloc(max(m,n),          sizeof(double));
  *coeff_low         = (double *) rwt_calloc(ncoeff,            sizeof(double));
  *coeff_high        = (double *) rwt_calloc(ncoeff,            sizeof(double));
}


/*!
 * Free memory that we allocated for dwt
 *
 * @param x_dummy_low
 * @param x_dummy_high
 * @param y_dummy_low_low
 * @param y_dummy_low_high
 * @param y_dummy_high_low
 * @param y_dummy_high_high
 * @param coeff_low
 * @param coeff_high
 *
 */
void rdwt_free(double **x_dummy_low, double **x_dummy_high, double **y_dummy_low_low, double **y_dummy_low_high, 
  double **y_dummy_high_low, double **y_dummy_high_high, double **coeff_low, double **coeff_high) {
  rwt_free(*x_dummy_low);
  rwt_free(*x_dummy_high);
  rwt_free(*y_dummy_low_low);
  rwt_free(*y_dummy_low_high);
  rwt_free(*y_dummy_high_low);
  rwt_free(*y_dummy_high_high);
  rwt_free(*coeff_low);
  rwt_free(*coeff_high);
}


/*!
 * Put the scaling coeffients into a form ready for use in the convolution function
 *
 * @param ncoeff length of h / the number of scaling coefficients
 * @param h  the wavelet scaling coefficients
 * @param coeff_low the high pass coefficients - reversed h
 * @param coeff_high the high pass coefficients - forward h, even values are sign flipped
 *
 * The coefficients of our Quadrature Mirror Filter are described by
 * \f$ g\left[lh - 1 - n \right] = (-1)^n * h\left[n\right] \f$
 *
 * This is identical to dwt_coefficients() 
 *
 */
void rdwt_coefficients(int ncoeff, double *h, double **coeff_low, double **coeff_high) {
  int i;
  for (i=0; i<ncoeff; i++) {
    (*coeff_low)[i] = h[(ncoeff-i)-1];
    (*coeff_high)[i] = h[i];
  }
  for (i=0; i<ncoeff; i+=2)
    (*coeff_high)[i] = -((*coeff_high)[i]);
}


/*!
 * Perform the redundant discrete wavelet transform
 *
 * @param x      the input signal
 * @param nrows  number of rows in the input
 * @param ncols  number of columns in the input
 * @param h      wavelet scaling coefficients
 * @param ncoeff length of h / the number of scaling coefficients
 * @param levels the number of levels
 * @param yl
 * @param yh
 *
 */
void rdwt(double *x, size_t nrows, size_t ncols, double *h, int ncoeff, int levels, double *yl, double *yh) {
  double *coeff_low, *coeff_high, *y_dummy_low_low, *y_dummy_low_high, *y_dummy_high_low;
  double *y_dummy_high_high, *x_dummy_low, *x_dummy_high;
  long i;
  int current_level, three_n_L, sample_f;
  size_t current_rows, current_cols, idx_rows, idx_cols, n_c, n_cb, n_r, n_rb;
  size_t column_cursor, column_cursor_plus_n, column_cursor_plus_double_n;

  rdwt_allocate(nrows, ncols, ncoeff, &x_dummy_low, &x_dummy_high, &y_dummy_low_low, &y_dummy_low_high, 
    &y_dummy_high_low, &y_dummy_high_high, &coeff_low, &coeff_high);

  rdwt_coefficients(ncoeff, h, &coeff_low, &coeff_high);

  if (ncols==1) {
    ncols = nrows;
    nrows = 1;
  }  

  /* analysis lowpass and highpass */
 
  three_n_L = 3*ncols*levels;
  current_rows = 2*nrows;
  current_cols = 2*ncols;
  for (i=0; i<nrows*ncols; i++)
    yl[i] = x[i];
  
  /* main loop */
  sample_f = 1;
  for (current_level=1; current_level <= levels; current_level++) {
    current_rows = current_rows/2;
    current_cols = current_cols/2;
    /* actual (level dependent) column offset */
    if (nrows==1)
      column_cursor = ncols*(current_level-1);
    else
      column_cursor = 3*ncols*(current_level-1);
    column_cursor_plus_n = column_cursor + ncols;
    column_cursor_plus_double_n = column_cursor_plus_n + ncols;
    
    /* go by rows */
    n_cb = ncols/current_cols;         /* # of column blocks per row */
    for (idx_rows=0; idx_rows<nrows; idx_rows++) { /* loop over rows */
      for (n_c=0; n_c<n_cb; n_c++) {          /* loop within one row */      
	/* store in dummy variable */
	idx_cols = -sample_f + n_c;
	for (i=0; i<current_cols; i++) {
	  idx_cols = idx_cols + sample_f;
	  x_dummy_low[i] = mat(yl, idx_rows, idx_cols, nrows, ncols);  
	}
	/* perform filtering lowpass/highpass */
	rdwt_convolution(x_dummy_low, current_cols, coeff_low, coeff_high, ncoeff, y_dummy_low_low, y_dummy_high_high); 
	/* restore dummy variables in matridx_colses */
	idx_cols = -sample_f + n_c;
	for  (i=0; i<current_cols; i++) {
          idx_cols = idx_cols + sample_f;
          mat(yl, idx_rows, idx_cols,                 nrows, ncols)     = y_dummy_low_low[i];
          mat(yh, idx_rows, idx_cols + column_cursor, nrows, three_n_L) = y_dummy_high_high[i];  
	} 
      }
    }
      
    /* go by columns in case of a 2D signal*/
    if (nrows>1) {
      n_rb = nrows/current_rows;           /* # of row blocks per column */
      for (idx_cols=0; idx_cols<ncols; idx_cols++) { /* loop over column */
	for (n_r=0; n_r<n_rb; n_r++) {         /* loop within one column */
	  /* store in dummy variables */
	  idx_rows = -sample_f + n_r;
	  for (i=0; i<current_rows; i++) {    
	    idx_rows = idx_rows + sample_f;
	    x_dummy_low[i]  = mat(yl, idx_rows, idx_cols,                 nrows, ncols);
	    x_dummy_high[i] = mat(yh, idx_rows, idx_cols + column_cursor, nrows, three_n_L);
	  }
	  /* perform filtering: first LL/LH, then HL/HH */
	  rdwt_convolution(x_dummy_low,  current_rows, coeff_low, coeff_high, ncoeff, y_dummy_low_low,  y_dummy_low_high);
	  rdwt_convolution(x_dummy_high, current_rows, coeff_low, coeff_high, ncoeff, y_dummy_high_low, y_dummy_high_high);
	  /* restore dummy variables in matrices */
	  idx_rows = -sample_f + n_r;
	  for (i=0; i<current_rows; i++) {
	    idx_rows = idx_rows + sample_f;
	    mat(yl, idx_rows, idx_cols,                               nrows, ncols)     = y_dummy_low_low[i];
	    mat(yh, idx_rows, idx_cols + column_cursor,               nrows, three_n_L) = y_dummy_low_high[i];
	    mat(yh, idx_rows, idx_cols + column_cursor_plus_n,        nrows, three_n_L) = y_dummy_high_low[i];
	    mat(yh, idx_rows, idx_cols + column_cursor_plus_double_n, nrows, three_n_L) = y_dummy_high_high[i];
	  }
	}
      }
    }
    sample_f = sample_f*2;
  }
  rdwt_free(&x_dummy_low, &x_dummy_high, &y_dummy_low_low, &y_dummy_low_high, &y_dummy_high_low, &y_dummy_high_high, &coeff_low, &coeff_high);
}
