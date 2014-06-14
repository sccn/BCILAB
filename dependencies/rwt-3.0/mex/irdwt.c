/*! \file irdwt.c
    \brief Implementation of the inverse redundant discrete wavelet transform

*/

#include "rwt_platform.h"

void irdwt_convolution(double *x_out, size_t lx, double *coeff_low, double *coeff_high, int ncoeff, double *x_in_low, double *x_in_high) {
  int k;
  size_t i, j;
  double x0;

  for (k=ncoeff-2; k > -1; k--) {
    x_in_low[k] = x_in_low[lx+k];
    x_in_high[k] = x_in_high[lx+k];
  }
  for (i=0; i<lx; i++){
    x0 = 0;
    for (j=0; j<ncoeff; j++)
      x0 = x0 + (x_in_low[j+i] * coeff_low[ncoeff-1-j]) + (x_in_high[j+i] * coeff_high[ncoeff-1-j]);
	
    x_out[i] = x0;
  }
}


void irdwt_allocate(size_t m, size_t n, int ncoeff, double **x_high, double **x_dummy_low, double **x_dummy_high, double **y_dummy_low_low, 
  double **y_dummy_low_high, double **y_dummy_high_low, double **y_dummy_high_high, double **coeff_low, double **coeff_high) {
  *x_high            = (double *) rwt_calloc(m*n,               sizeof(double));
  *x_dummy_low       = (double *) rwt_calloc(max(m,n),          sizeof(double));
  *x_dummy_high      = (double *) rwt_calloc(max(m,n),          sizeof(double));
  *y_dummy_low_low   = (double *) rwt_calloc(max(m,n)+ncoeff-1, sizeof(double));
  *y_dummy_low_high  = (double *) rwt_calloc(max(m,n)+ncoeff-1, sizeof(double));
  *y_dummy_high_low  = (double *) rwt_calloc(max(m,n)+ncoeff-1, sizeof(double));
  *y_dummy_high_high = (double *) rwt_calloc(max(m,n)+ncoeff-1, sizeof(double));
  *coeff_low         = (double *) rwt_calloc(ncoeff,            sizeof(double));
  *coeff_high        = (double *) rwt_calloc(ncoeff,            sizeof(double));
}


void irdwt_free(double **x_dummy_low, double **x_dummy_high, double **y_dummy_low_low, double **y_dummy_low_high, 
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

/* not the same as idwt_coefficients */
void irdwt_coefficients(int ncoeff, double *h, double **coeff_low, double **coeff_high) {
  int i;
  for (i=0; i<ncoeff; i++) {
    (*coeff_low)[i] = h[i]/2;
    (*coeff_high)[i] = h[ncoeff-i-1]/2;
  }
  for (i=1; i<=ncoeff; i+=2)
    (*coeff_high)[i] = -((*coeff_high)[i]);
}


void irdwt(double *x, size_t nrows, size_t ncols, double *h, int ncoeff, int levels, double *y_low, double *y_high) {
  double  *coeff_low, *coeff_high, *y_dummy_low_low, *y_dummy_low_high, *y_dummy_high_low;
  double *y_dummy_high_high, *x_dummy_low , *x_dummy_high, *x_high;
  long i;
  int current_level, three_n_L, ncoeff_minus_one, sample_f;
  size_t current_rows, current_cols, column_cursor, column_blocks_per_row;
  size_t idx_rows, idx_cols, n_r, n_c;
  size_t row_blocks_per_column, column_cursor_plus_n, column_cursor_plus_double_n;

  irdwt_allocate(nrows, ncols, ncoeff, &x_high, &x_dummy_low, &x_dummy_high, &y_dummy_low_low, 
    &y_dummy_low_high, &y_dummy_high_low, &y_dummy_high_high, &coeff_low, &coeff_high);
  irdwt_coefficients(ncoeff, h, &coeff_low, &coeff_high);
 
  if (ncols==1) {
    ncols = nrows;
    nrows = 1;
  }
  /* analysis lowpass and highpass */
  
  three_n_L = 3*ncols*levels;
  ncoeff_minus_one = ncoeff - 1;
  /*! we calculate sample_f = 2^(levels - 1) with a loop since that is actually the recommended method for whole numbers */
  sample_f = 1;
  for (i=1; i<levels; i++)
    sample_f = sample_f*2;

  current_rows = nrows/sample_f;
  current_cols = ncols/sample_f;
  /* restore y_low in x */
  for (i=0; i<nrows*ncols; i++)
    x[i] = y_low[i];
  
  /* main loop */
  for (current_level=levels; current_level >= 1; current_level--) {
    /* actual (level dependent) column offset */
    if (nrows==1)
      column_cursor = ncols*(current_level-1);
    else
      column_cursor = 3*ncols*(current_level-1);
    column_cursor_plus_n = column_cursor + ncols;
    column_cursor_plus_double_n = column_cursor_plus_n + ncols;
    
    /* go by columns in case of a 2D signal*/
    if (nrows>1) {
      row_blocks_per_column = nrows/current_rows;   /* # of row blocks per column */
      for (idx_cols=0; idx_cols<ncols; idx_cols++) {          /* loop over column */
	for (n_r=0; n_r<row_blocks_per_column; n_r++) { /* loop within one column */
	  /* store in dummy variables */
	  idx_rows = -sample_f + n_r;
	  for (i=0; i<current_rows; i++) {    
	    idx_rows = idx_rows + sample_f;
	    y_dummy_low_low[i+ncoeff_minus_one]   = mat(x,      idx_rows, idx_cols,                               nrows, ncols);
	    y_dummy_low_high[i+ncoeff_minus_one]  = mat(y_high, idx_rows, idx_cols + column_cursor,               nrows, three_n_L);
	    y_dummy_high_low[i+ncoeff_minus_one]  = mat(y_high, idx_rows, idx_cols + column_cursor_plus_n,        nrows, three_n_L);
	    y_dummy_high_high[i+ncoeff_minus_one] = mat(y_high, idx_rows, idx_cols + column_cursor_plus_double_n, nrows, three_n_L);
	  }
	  /* perform filtering and adding: first LL/LH, then HL/HH */
	  irdwt_convolution(x_dummy_low,  current_rows, coeff_low, coeff_high, ncoeff, y_dummy_low_low,  y_dummy_low_high); 
	  irdwt_convolution(x_dummy_high, current_rows, coeff_low, coeff_high, ncoeff, y_dummy_high_low, y_dummy_high_high); 
	  /* store dummy variables in matrices */
	  idx_rows = -sample_f + n_r;
	  for (i=0; i<current_rows; i++) {
	    idx_rows = idx_rows + sample_f;
	    mat(x,      idx_rows, idx_cols, nrows, ncols) = x_dummy_low[i];
	    mat(x_high, idx_rows, idx_cols, nrows, ncols) = x_dummy_high[i];
	  }
	}
      }
    }
    
    /* go by rows */
    column_blocks_per_row = ncols/current_cols; /* # of column blocks per row */
    for (idx_rows=0; idx_rows<nrows; idx_rows++) {          /* loop over rows */
      for (n_c=0; n_c<column_blocks_per_row; n_c++) {  /* loop within one row */      
	/* store in dummy variable */
	idx_cols = -sample_f + n_c;
	for  (i=0; i<current_cols; i++) {    
	  idx_cols = idx_cols + sample_f;
	  y_dummy_low_low[i+ncoeff_minus_one] = mat(x, idx_rows, idx_cols, nrows, ncols);  
	  if (nrows>1)
	    y_dummy_high_high[i+ncoeff_minus_one] = mat(x_high, idx_rows, idx_cols,                 nrows, ncols);
	  else
            y_dummy_high_high[i+ncoeff_minus_one] = mat(y_high, idx_rows, idx_cols + column_cursor, nrows, three_n_L);
	} 
	/* perform filtering lowpass/highpass */
	irdwt_convolution(x_dummy_low, current_cols, coeff_low, coeff_high, ncoeff, y_dummy_low_low, y_dummy_high_high); 
	/* restore dummy variables in matrices */
	idx_cols = -sample_f + n_c;
	for (i=0; i<current_cols; i++) {    
	  idx_cols = idx_cols + sample_f;
	  mat(x, idx_rows, idx_cols, nrows, ncols) = x_dummy_low[i];  
	}
      }
    }
    sample_f = sample_f/2;
    current_rows = current_rows*2;
    current_cols = current_cols*2;
  }
  irdwt_free(&x_dummy_low, &x_dummy_high, &y_dummy_low_low, &y_dummy_low_high, &y_dummy_high_low, &y_dummy_high_high, &coeff_low, &coeff_high);
}

