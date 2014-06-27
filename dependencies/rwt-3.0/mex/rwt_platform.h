/*! \file rwt_platform.h
    \brief Abstract away environment differences and provide some common macros
*/
#ifndef RWT_PLATFORM_H
#define RWT_PLATFORM_H

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

/*! For MATLAB we address 2d inputs and outputs in column-major order */
/*! For Python we address 2d inputs and outputs in row-major order */
/*! The offset macros are for debugging */
/*! The parameters for the mat() macro are:
 *    a - the base pointer to the matrix of values
 *    i - index of the target row
 *    j - index of the target column
 *    m - the number of rows
 *    n - the number of columns
 */
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
  #define COLUMN_MAJOR_ORDER 1
  #include "mex.h"
  #ifndef OCTAVE_MEX_FILE
    #include "matrix.h"
  #endif
  #define mat(a, i, j, m, n) (*(a + (m*(j)+i)))
  #define mat_offset(a, i, j, m, n) (m*(j)+i)
  #define offset_row(offset, m, n) (offset % m)
  #define offset_col(offset, m, n) ((offset - (offset % m)) / m)
  #define rwt_printf(fmt, ...) mexPrintf(fmt, ##__VA_ARGS__)
  #define rwt_errormsg(msg) mexErrMsgTxt(msg)
#else
  #define ROW_MAJOR_ORDER 1
  #define mat(a, i, j, m, n) (*(a + (n*(i)+j)))
  #define mat_offset(a, i, j, m, n) (n*(i)+j)
  #define offset_row(offset, m, n) ((offset - (offset % n)) / n)
  #define offset_col(offset, m, n) (offset % n)
  #define rwt_printf(fmt, ...) printf(fmt, ##__VA_ARGS__)
  #define rwt_errormsg(msg) printf("\033[91m%s\033[0m\n", msg);
#endif

#define max(A,B) (A > B ? A : B)
#define min(A,B) (A < B ? A : B)
#define even(x)  ((x & 1) ? 0 : 1)

#ifdef __cplusplus
extern "C" {
#endif

void *rwt_malloc(size_t size);
void *rwt_calloc(size_t num, size_t size);
void rwt_free(void *ptr);

#ifdef __cplusplus
}
#endif

#endif
