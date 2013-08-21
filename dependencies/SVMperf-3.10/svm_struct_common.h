/***********************************************************************/
/*                                                                     */
/*   svm_struct_common.h                                               */
/*                                                                     */
/*   Functions and types used by multiple components of SVM-struct.    */
/*                                                                     */
/*   Author: Thorsten Joachims                                         */
/*   Date: 31.10.05                                                    */
/*                                                                     */
/*   Copyright (c) 2005  Thorsten Joachims - All rights reserved       */
/*                                                                     */
/*   This software is available for non-commercial use only. It must   */
/*   not be modified and distributed without prior permission of the   */
/*   author. The author is not responsible for implications from the   */
/*   use of this software.                                             */
/*                                                                     */
/***********************************************************************/

#ifndef svm_struct_common
#define svm_struct_common

# define STRUCT_VERSION       "V3.10"
# define STRUCT_VERSION_DATE  "14.08.08"

#ifdef __cplusplus
extern "C" {
#endif
#include "svm_common.h"
#ifdef __cplusplus
}
#endif
#include "svm_struct_api_types.h"

typedef struct example {  /* an example is a pair of pattern and label */
	PATTERN x;
	LABEL y;
} EXAMPLE;

typedef struct sample { /* a sample is a set of examples */
	int     n;            /* n is the total number of examples */
	EXAMPLE *examples;
} SAMPLE;

typedef struct constset { /* a set of linear inequality constrains of
 for lhs[i]*w >= rhs[i] */
	int     m;            /* m is the total number of constrains */
	DOC     **lhs;
	double  *rhs;
} CONSTSET;


/**** print methods ****/
void printIntArray(int*,int);
void printDoubleArray(double*,int);
void printWordArray(WORD*);
void printModel(MODEL *);
void printW(double *, long, long, double);

extern long   struct_verbosity;              /* verbosity level (0-4) */

#endif
