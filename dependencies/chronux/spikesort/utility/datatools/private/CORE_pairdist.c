/*
 * CORE_pairdist.c
 * --------------------
 * The syntax for this function is described in the M file of
 * the same name.
 *
 * This is a MEX-file for MATLAB.  It was compiled using:
 *
 *      COREmex CORE_pairdist.c
 * 
 * samar@cs.stanford.edu, 07/12/05
 */

#include "mex.h"
#include "matrix.h"
#include "CORE_library.h"
#include "CORE_mextools.h"


/* Variables */
static mxArray *psst_dists = NULL;


/* Function: psst_cleanup()
 * ----------------------------------------------------------
 * Cleans up persistent memory.
 *
 */

void psst_cleanup(void) {
  psst_dists = safe_destroy_array(psst_dists);
}


/* ********************************************************** */
/* ********************* MEX GATEWAY  *********************** */
/* ********************************************************** */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  int mrows, nrows, dcols, m_psst, n_psst;
  double *vectsX, *normsqrX, *vectsY, *normsqrY, *dists;     //data
  double takeSqrt, reuseMem, useSafe;                        //flags

  /* ************* ACCESS MATLAB VARIABLES ************* */
  /* Get all of the flags first -- no error checking, so they better be ok */
  takeSqrt      = mxGetScalar(prhs[2]);
  reuseMem      = mxGetScalar(prhs[3]);
  useSafe       = mxGetScalar(prhs[4]);

  /* Next get pointer to vector lists */
  vectsX  = (double *) mxGetPr(prhs[0]);   vectsY = (double *)mxGetPr(prhs[1]);
  mrows = mxGetM(prhs[0]);   dcols = mxGetN(prhs[0]);    nrows = mxGetM(prhs[1]);

  if (reuseMem) {
    psst_dists = persistent_matrix(psst_dists, mrows, nrows);
    mexAtExit(psst_cleanup);
    plhs[0] = psst_dists;
  } else {
    psst_dists = safe_destroy_array(psst_dists);
    plhs[0] = mxCreateNumericMatrix(mrows, nrows, mxDOUBLE_CLASS, mxREAL);
  }
  dists = (double *) mxGetPr(plhs[0]);
  

  if (!useSafe) {
    /* Compute vector norms for fast algorithm */
    normsqrX = (double *) mxMalloc(mrows * sizeof(double));
    norm_squared_rows(vectsX, mrows, dcols, normsqrX);
    
    normsqrY = (double *) mxMalloc(nrows * sizeof(double));
    norm_squared_rows(vectsY, nrows, dcols, normsqrY);
  }
  
  /* ********************* COMPUTE ********************* */
  if (!useSafe) {
    pairdist_squared_rows(vectsX, normsqrX, vectsY, normsqrY, dists, 
			  mrows, nrows, dcols);
  } else {
    pairdist_squared_rows_safe(vectsX, vectsY, dists, mrows, nrows, dcols);
  }

  if (takeSqrt) {sqrt_matrix_inplace(dists, mrows, nrows);}


  /* ********************* CLEAN UP ******************** */
  if (!useSafe) {
    mxFree(normsqrX); 
    mxFree(normsqrY);
  }

  // Additional static allocated memory is freed when this MEX file
  // is cleared from memory in the function psst_cleanup above.
}
