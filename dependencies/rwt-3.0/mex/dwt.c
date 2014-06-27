/*! \file dwt.c
    \brief Implementation of the discrete wavelet transform

*/

#include "rwt_platform.h"

/*!
 * Perform convolution for dwt
 *
 * @param x_in input signal values
 * @param lx the length of x_in
 * @param coeff_low the low pass coefficients
 * @param coeff_high the high pass coefficients
 * @param ncoeff_minus_one one less than the number of scaling coefficients
 * @param x_out_low low pass results
 * @param x_out_high high pass results
 * 
 * For the convolution we will calculate the output of the lowpass and highpass filters in parallel
 *
 * Normally we can describe the calculation of a convolution as
 * \f$ (\textbf{w} * \textbf{z})_k = \frac{1}{N} \sum\limits_{l=0}^{2N-1} w_{k-l} \cdot z_{l} \f$
 *
 * Our actual implementation resembles this
 *
 */
void dwt_convolution(double *x_in, size_t lx, double *coeff_low, double *coeff_high, int ncoeff_minus_one, double *x_out_low, double *x_out_high) {
  size_t i, j, ind;
  double x0, x1;
  for (i=lx; i<lx+ncoeff_minus_one; i++) { 
    x_in[i] = *(x_in+(i-lx)); /*! extend x_in by creating a small mirror at the end of length ncoeff_minus_one */
  }
  ind = 0;
  for (i=0; i<(lx); i+=2) {   /*! Step through the input values, moving right 2 values each loop */
    x0 = 0;
    x1 = 0;
    for (j=0; j<=ncoeff_minus_one; j++) {                   /*! Take the high and low filters in reverse order */
      x0 = x0 + x_in[i+j] * coeff_low[ncoeff_minus_one-j];  /*! Sum the product of the next ncoeff values of x_in with the filter coefficients */
      x1 = x1 + x_in[i+j] * coeff_high[ncoeff_minus_one-j];
    }
    x_out_low[ind] = x0; /*! Place these calculated sums in the next position of the output */
    x_out_high[ind++] = x1;
  }
}


/*!
 * Allocate memory for dwt
 *
 * @param m      the number of rows of the input matrix
 * @param n      the number of columns of the input matrix
 * @param ncoeff the number of scaling coefficients
 * @param x_dummy      storage space for input data being passed to the convolution
 * @param y_dummy_low  storage space for low pass convolution results
 * @param y_dummy_high storage space for high pass convolution results
 * @param coeff_low    storage space for the low pass coefficients
 * @param coeff_high   storage space for the high pass coefficients
 *
 * The low pass and high pass filter coefficients are the same size as the scaling coefficients
 * For the output storage area we will need as much space as the input: m*n
 * For the input storage area we will need the same plus one less than the length of the coeffiecients
 */
void dwt_allocate(size_t m, size_t n, int ncoeff, double **x_dummy, double **y_dummy_low, double **y_dummy_high, double **coeff_low, double **coeff_high) {
  *x_dummy      = (double *) rwt_calloc(max(m,n)+ncoeff-1, sizeof(double));
  *y_dummy_low  = (double *) rwt_calloc(max(m,n),          sizeof(double));
  *y_dummy_high = (double *) rwt_calloc(max(m,n),          sizeof(double));
  *coeff_low    = (double *) rwt_calloc(ncoeff,            sizeof(double));
  *coeff_high   = (double *) rwt_calloc(ncoeff,            sizeof(double));
}


/*!
 * Free memory that we allocated for dwt
 *
 * @param x_dummy      storage space for input data being passed to the convolution
 * @param y_dummy_low  storage space for low pass convolution results
 * @param y_dummy_high storage space for high pass convolution results
 * @param coeff_low    storage space for the low pass coefficients
 * @param coeff_high   storage space for the high pass coefficients
 *
 */
void dwt_free(double **x_dummy, double **y_dummy_low, double **y_dummy_high, double **coeff_low, double **coeff_high) {
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
 * @param coeff_low the low pass coefficients - reversed h
 * @param coeff_high the high pass coefficients - forward h, alternate values are sign flipped
 *
 * The coefficients of our Quadrature Mirror Filter are described by
 * \f$ g\left[lh - 1 - n \right] = (-1)^n * h\left[n\right] \f$
 *
 */
void dwt_coefficients(int ncoeff, double *h, double **coeff_low, double **coeff_high) {
  int i;
  for (i=0; i<ncoeff; i++) {
    (*coeff_low)[i] = h[(ncoeff-i)-1];
    (*coeff_high)[i] = h[i];
  }
  for (i=0; i<ncoeff; i+=2)
    (*coeff_high)[i] = -((*coeff_high)[i]);
}


/*!
 * Perform the discrete wavelet transform
 *
 * @param x      the input signal
 * @param nrows  number of rows in the input
 * @param ncols  number of columns in the input
 * @param h      wavelet scaling coefficients
 * @param ncoeff length of h / the number of scaling coefficients
 * @param levels the number of levels
 * @param y      the output signal with the wavelet transform applied
 *
 * The discrete wavelet transform begins with a set of samples of a signal whose length
 * is a power of 2. This exponent will be the maximum number of levels of the transform
 * that we can perform.
 *
 */
void dwt(double *x, size_t nrows, size_t ncols, double *h, int ncoeff, int levels, double *y) {
  double  *coeff_low, *coeff_high, *y_dummy_low, *y_dummy_high, *x_dummy;
  long i;
  int current_level, ncoeff_minus_one;
  size_t current_rows, current_cols, row_cursor, column_cursor, idx_rows, idx_columns;

  if (ncols==1) { /*! Accept either column vectors or row vectors. Store the length in the variable n */
    ncols = nrows;
    nrows = 1;
  }
  
  dwt_allocate(nrows, ncols, ncoeff, &x_dummy, &y_dummy_low, &y_dummy_high, &coeff_low, &coeff_high);
  dwt_coefficients(ncoeff, h, &coeff_low, &coeff_high); /*! For performance, calculate what we can outside the loops */
  ncoeff_minus_one = ncoeff - 1;
  current_rows = 2*nrows; /*! current_rows and current_cols start at 2x since we divide by 2 at the start of the loop */
  current_cols = 2*ncols;
 
  for (current_level=1; current_level<=levels; current_level++) {
    if (nrows==1)
      current_rows = 1;
    else{
      current_rows = current_rows/2;
      row_cursor = current_rows/2;     
    }
    current_cols = current_cols/2;
    column_cursor = current_cols/2;

    for (idx_rows=0; idx_rows<current_rows; idx_rows++) {
      for (i=0; i<current_cols; i++)
	if (current_level==1)  
	  x_dummy[i] = mat(x, idx_rows, i, nrows, ncols);  
	else 
	  x_dummy[i] = mat(y, idx_rows, i, nrows, ncols);  
      /*! Perform filtering lowpass and highpass*/
      dwt_convolution(x_dummy, current_cols, coeff_low, coeff_high, ncoeff_minus_one, y_dummy_low, y_dummy_high); 
      /*! Restore dummy variables in matrices */
      idx_columns = column_cursor;
      for (i=0; i<column_cursor; i++) {    
	mat(y, idx_rows, i,             nrows, ncols) = y_dummy_low[i];  
	mat(y, idx_rows, idx_columns++, nrows, ncols) = y_dummy_high[i];  
      } 
    }  
    
    /*! For the 2d transform, we go through each of the columns after having gone through the rows */
    if (nrows>1) {
      for (idx_columns=0; idx_columns<current_cols; idx_columns++) { /* loop over columns */
	/*! Store in dummy variables */
	for (i=0; i<current_rows; i++)
	  x_dummy[i] = mat(y, i, idx_columns, nrows, ncols);  
	/*! Perform filtering lowpass and highpass*/
	dwt_convolution(x_dummy, current_rows, coeff_low, coeff_high, ncoeff_minus_one, y_dummy_low, y_dummy_high); 
	/*! Restore dummy variables in matrix */
	idx_rows = row_cursor;
	for (i=0; i<row_cursor; i++) {
	  mat(y, i,          idx_columns, nrows, ncols) = y_dummy_low[i];  
	  mat(y, idx_rows++, idx_columns, nrows, ncols) = y_dummy_high[i];  
	}
      }
    }
  }
  dwt_free(&x_dummy, &y_dummy_low, &y_dummy_high, &coeff_low, &coeff_high);
}

