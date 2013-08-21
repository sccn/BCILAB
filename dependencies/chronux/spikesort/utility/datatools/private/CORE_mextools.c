/*
 * FILE: CORE_mextools.c
 * --------------------
 * Common MEX utilities for Matlab CORE_function_ MEX files.  
 * Function descriptions can be found in CORE_mextools.h
 * 
 * samar@cs.stanford.edu, 09/18/05
 */

#include "CORE_mextools.h"

// ------------------------------------------------------------------------
// Function: persistent_matrix
// ------------------------------------------------------------------------
mxArray *safe_destroy_array(mxArray *ptrMatrix)
{
  if (ptrMatrix != NULL) {mxDestroyArray(ptrMatrix);}
  return NULL;
}
// ------------------------------------------------------------------------


// ------------------------------------------------------------------------
// Function: persistent_matrix
// ------------------------------------------------------------------------
mxArray *persistent_matrix(mxArray *ptrMatrix, int M, int N)
{
  int oldM, oldN;
  double *newspace;

  if (ptrMatrix != NULL) { // need to force new allocation if size larger
    oldM = mxGetM(ptrMatrix);  oldN = mxGetN(ptrMatrix);
    if ((oldM < M) || (oldN < N)) {   // kill the old array; its too small
      mxDestroyArray(ptrMatrix);    ptrMatrix = NULL;
    } else if ((oldM > M) || (oldN > N)) { // old array too big; try realloc
      newspace = (double *) mxRealloc(mxGetPr(ptrMatrix), M*N*sizeof(double));
      if (newspace == NULL) {
	mxDestroyArray(ptrMatrix);  ptrMatrix = NULL;  // realloc failed ...
      } else {
	mxSetM(ptrMatrix, M);  mxSetN(ptrMatrix, N);   // realloc succeeded ...
	mxSetPr(ptrMatrix, newspace);
      }
    }
  }

  if (ptrMatrix == NULL) { // if still no memory allocated, allocate it
    ptrMatrix = mxCreateNumericMatrix(M, N, mxDOUBLE_CLASS, mxREAL);
    mexMakeArrayPersistent(ptrMatrix);
  }

  return ptrMatrix;
}
// ------------------------------------------------------------------------
