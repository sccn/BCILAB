/*
 * FILE: CORE_library.h
 * --------------------
 * Common computational routines for Matlab CORE_function_ MEX files.
 * 
 * samar@cs.stanford.edu, 01/18/06
 */

#ifndef __CORE_LIBRARY__
#define __CORE_LIBRARY__

/* --------------------------  EXTERNAL FUNCTIONS --------------------------- */
#ifdef matrix_h   /* If linking to a file using MEX, translate to Matlab */
/* Needed for: find_minmax */
unsigned char IsNaN(double value)  {return mxIsNaN(value);}
double GetNaN(void)  {return mxGetNaN();}
#else             /* If no MEX, prototypes only; define depending on compiler, etc. */
/* Needed for: find_minmax */
unsigned char IsNaN(double value);
double GetNaN(void);
#endif

/* ---------------------------- PUBLIC FUNCTIONS ---------------------------- */
/* SUMMARY:
 *   count_2D                    2D histogram
 *   count_diffs                 Pairwise time differences
 *   count_in_cols               Column-by-column histogram
 *   find_minmax                 Simult. min/max, values & indices (ignore NaN)
 ****k_centers_approx            K clusters of min (within x2) possible radius
 ****k_means                     K clusters with min (local opt) mse
 *   minumum_column_indices      same as 2nd output of Matlab's min(A,[],2)
 *   minimum_distances           Closest vectors from two lists
 *   norm_squared_rows           same as Matlab's sum(A.^2,2)
 *   pairdist_squared_rows       Pairwise vector distances
 *   pairdist_squared_rows_safe  Pairwise vector distances (alternative version)
 *   sqrt_matrix_inplace         Replace elements with their square roots
 */


/*
 * Function:  count_2D
 *    Input:  arrays rowinds and colinds of length N with integer-valued
 *            entries; all entries in rowinds must be between 1 and maxrow,
 *            all entries in colinds must be between 1 and maxcol
 *    Output: column-major integer-valued counts with maxrow rows and
 *            maxcol columns, where, for the (jth row,kth col) of counts,
 *               counts[j][k] = #(rowinds[i]==j && colinds[i]==k), 0<=i<N
 *    Memory: (maxrow * maxcol)*sizeof(double) bytes must be pre-allocated
 *            for counts
 */
void count_2D(double *rowinds, double *colinds, int N, 
              double *counts, int maxrow, int maxcol);


/*
 * Function:  count_diffs
 *    Input:  arrays x of length Nx and y of length Ny with integer-valued
 *            entries
 *                 >>> Both x and y must be sorted in ascending order <<<
 *    Output: diffs[i] =  #(x[j]-y[k] = i), for |i| <= maxdiff
 *    Memory: (2*maxdiff + 1)*sizeof(double) bytes must be pre-allocated
 *            for diffs
 */
void count_diffs(double *x, int Nx, double *y, int Ny, double *diffs, int maxdiff);


/*
 * Function:  count_in_cols
 *    Input:  column-major matrix with M rows and N columns with integer-
 *            valued entries between 1 and maxval
 *    Output: column-major integer-valued counts with maxval rows and N
 *            columns, where the jth column of counts is such that 
 *            counts[k][j] = #(matrix[i][j] = k), 0 <= i < M
 *    Memory: (maxval * N)*sizeof(double) bytes must be pre-allocated for
 *            counts
 */
void count_in_cols(double *matrix, int M, int N, double *counts, int maxval);


/*
 * Function:  find_minmax
 *    Input:  array of length M
 *    Output: minimum (minval) and maximum (maxval) values of input and
 *            their integer-valued, 0-indexed indices (minind, maxind) in 
 *            the input array -- ties are broken in favor of the earlier
 *            entry.
 *    Memory: sizeof(double) bytes must be pre-allocated for minval, maxval
 *            sizeof(int) bytes must be pre-allocated for minind, maxind
 *    Special: IEC60559/IEEE754 +Inf/-Inf values behave as greater/less
 *             than all other finite values respectively.
 *             NaN values are more difficult because of compiler variation.
 *             The desired behavior is to ignore NaN values unless the
 *             input is all NaN, in which case min/max values are NaN and
 *             min/max indices are 0.  To accomplish this, the caller must
 *             define two functions: IsNaN, which returns 1 or 0, and
 *             GetNaN, which returns the value to be treated as NaN.
 *             The following should be used if NaN will never appear:
 *                unsigned char IsNaN(double value) {return 0;}
 *                double GetNaN(void) {return 0;}
 *             These should be used with Matlab/MEX:
 *                unsigned char IsNaN(double value) {return mxIsNaN(value);}
 *                double GetNaN(void) {return mxGetNaN();}
 */
void find_minmax(double *input, int M, 
                 double *minval, int *minind, double *maxval, int *maxind);


/*
 * Function:  k_centers_approx
 *    Input:  column-major matrix vects with M rows and N columns
 *            Target number of centers K
 *    Output: column major matrix centers with K rows and N columns
 *    Memory: K*N*sizeof(double) bytes must be pre-allocated for centers
 */
void k_centers_approx(double *vects, int M, int N, int K, double *centers);


/*
 * Function:  k_means
 *    Input:  column-major matrix vects with M rows and N columns
 *            Target number of centroids K
 *    Output: column major matrix centroids with K rows and N columns
 *    Memory: K*N*sizeof(double) bytes must be pre-allocated for centroids
 */
void k_means(double *vects, int M, int N, int K, double *centroids);


/*
 * Function:  minimum_column_indices
 *    Input:  column-major matrix with M rows and N columns
 *    Output: integer-valued length M array mininds such that
 *                mininds[j] = 1 + argmin_i (matrix[j][i]),  0 <= i < N
 *            Ties are broken in favor of the lowest-indexed column.
 *            NOTE that the indices in mininds are 1-indexed.
 *    Memory: M*sizeof(double) bytes must be pre-allocated for mininds
 */
void minimum_column_indices(double *matrix, int M, int N, double *mininds);


/*
 * Function:  minimum_distances
 *    Input:  column-major matrix vectsX with M rows and D columns
 *            column-major matrix vectsY with N rows and D columns
 *            normsqrX must be an array of length M where
 *                    normsqrX[j] = sum(vectsX[j][k]^2), 0 <= k < D.
 *            Similarly, normsqrY must be an array of length N where
 *                    normsqrY[j] = sum(vectsY[j][k]^2), 0 <= k < D.
 *    Output: vectors Yinds and dists, each with M rows, such that,
 *            for 0<=j<N,
 *                dists[i] = min_j (|| vectsX[i,:]-vectsY[j,:] ||^2)
 *                yinds[i] = 1 + argmin_j (|| vectsX[i,:]-vectsY[j,:] ||)
 *            All values yinds[i] will be integer-valued doubles in [1,N].
 *            Ties are broken in favor of the lowest-indexed column.
 *            NOTE that the indices in yinds are 1-indexed.
 *            NOTE that the returned distances are actually squared dists.
 *    Memory: M*sizeof(double) bytes must be pre-allocated for yinds
 *            M*sizeof(double) bytes must be pre-allocated for dists
 */
void minimum_distances(double *vectsX, double *normsqrX, double *vectsY, double *normsqrY,
                           int M, int N, int D, double *Yinds, double *distsqr);


/*
 * Function:  norm_squared_rows
 *    Input:  column-major matrix with M rows and N columns
 *    Output: length M array normsqr such that
 *                normsqr[j] = sum(matrix[j][k]^2), 0 <= k < N
 *    Memory: M*sizeof(double) bytes must be pre-allocated for normsqr
 */
void norm_squared_rows(double *matrix, int M, int N, double *normsqr);


/*
 * Function:  pairdist_squared_rows
 *    Input:  column-major matrix vectsX with M rows and D columns
 *            column-major matrix vectsY with N rows and D columns
 *            normsqrX must be an array of length M where
 *                    normsqrX[j] = sum(vectsX[j][k]^2), 0 <= k < D.
 *            Similarly, normsqrY must be an array of length N where
 *                    normsqrY[j] = sum(vectsY[j][k]^2), 0 <= k < D.
 *    Output: column-major dists with M rows and N columns such that
 *                dists[i][j] = sum((vectsX[i][k]-vectsY[j][k])^2)
 *            If the pointers vectsX and vectsY are the same, and
 *              M == N, all dists along the diagonal are forced to
 *              zero to correct for numerical errors.  These errors
 *              may remain, however, for any other case in which some
 *              row of vectsX is identical to a row of vectsY. If this
 *              is a possibility, the slower pairdist_squared_rows_safe
 *              (below) should be used instead.
 *    Memory: M*N*sizeof(double) bytes must be pre-allocated for dists
 */
void pairdist_squared_rows(double *vectsX, double *normsqrX, double *vectsY, double *normsqrY,
                           double *dists, int M, int N, int D);


/*
 * Function:  pairdist_squared_rows_safe
 *    Input:  column-major matrix vectsX with M rows and D columns
 *            column-major matrix vectsY with N rows and D columns
 *    Output: column-major dists with M rows and N columns such that
 *                dists[i][j] = sum((vectsX[i][k]-vectsY[j][k])^2)
 *            This function is slower than pairdist_squared_rows above,
 *                but correctly handles cases where dists[i][j] == 0.
 *    Memory: M*N*sizeof(double) bytes must be pre-allocated for dists
 */
void pairdist_squared_rows_safe(double *vectsX, double *vectsY, double *dists, 
                                int M, int N, int D);

/*
 * Function:  sqrt_matrix_inplace
 *    Input:  column-major matrix with M rows and N columns
 *    Output: elements of matrix are overwritten with their square roots
 *    Memory: (performed in-place)
 */
void sqrt_matrix_inplace(double *matrix, int M, int N);


/*---------------------------------------------------------------------------*/


#endif
