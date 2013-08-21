/*
 * FILE: CORE_mextools.h
 * ---------------------
 * Common MEX utilities for Matlab CORE_function_ MEX files.
 * 
 * samar@cs.stanford.edu, 09/18/05
 */

#ifndef __CORE_MEXTOOLS__
#define __CORE_MEXTOOLS__

#include "mex.h"
#include "matrix.h"

// ------------------------ UTILITY FUNCTIONS -------------------------------
/* SUMMARY:
 *   safe_destroy_array    calls mxDestroyArray if != NULL
 *   persistent_matrix     manages a persistent mxArray
 */


/* Function: safe_destroy_array
 * ----------------------------
 * Calls mxDestroyArray only if the input point is not NULL.
 * As a convenience, returns NULL so that it can be used in
 * e.g.,
 *    void cleanup(void) {ptr = safe_destroy_array(ptr);}
 * To guarantee that ptr is deallocated and set to NULL.
 */ 
mxArray *safe_destroy_array(mxArray *ptrMatrix);


/* Function: persistent_matrix
 * ---------------------------
 * This function assists in allocating persistent memory for
 * repeated MEX function calls.  This is useful when allocating
 * large amounts of memory if the size of the memory block is
 * unlikely to change between calls to the MEX function.
 * NOTE: This function requires that ptrMatrix point to a Matlab
 * matrix of type DOUBLE.
 * 
 * In general follow these steps:
 *   1) Declare a static mxArray * in the MEX .c file,
 *      and initialize it to NULL.
 *   2) Add a line to the cleanup function for the MEX
 *      file (or create one if necessary) calling
 *            ptr = safe_destroy_array(ptr);
 *      where ptr is your static mxArray.
 *   3) Call this function with desired number of rows (M) and
 *      columns (N) as:
 *            ptr = persistent_matrix(ptr, M, N);
 *      This will only allocate memory if ptr is NULL or if the
 *      requested M,N are different from the values currently
 *      allocated for the mxArray *ptr.
 *   4) Be sure to register your clean up function using:
 *            mexAtExit(cleanup);
 *   5) You can explicitly call safe_destroy_array(ptr) to
 *      de-allocate previously persistently allocated memory.
 */
mxArray *persistent_matrix(mxArray *ptrMatrix, int M, int N);

#endif
