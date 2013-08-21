#ifndef _GLOBAL_H
#define _GLOBAL_H

#include "mem_clean.h"


/****************************************
 * global malloc array - used to track 
 * blocks that are allocated my my_malloc 
 * and can be used to prevent double free()'ing
 * and to release memory at the end of execution.
 */
extern HASH_ARRAY *malloc_array;

/* from svm_hideo.c */
# define DEF_PRECISION          1E-5
# define DEF_MAX_ITERATIONS     200
# define DEF_LINDEP_SENSITIVITY 1E-8
# define EPSILON_HIDEO          1E-20
# define EPSILON_EQ             1E-5

/******************************************
 * from svm_hideo.c - made global 
 * to allow a single global initialization
 */
extern double *primal;
extern double *dual;
extern long   precision_violations;
extern double opt_precision;
extern long   maxiter;
extern double lindep_sensitivity;
extern double *buffer;
extern long   *nonoptimal;

extern long  smallroundcount;
extern long  roundnumber;

extern long   kernel_cache_statistic; 
extern long   verbosity;

/****************************************
 * declare a global initialization routine
 **************************************** 
*/

void global_init( );
void global_destroy( );

void show_model(MODEL *model);
void show_doc(DOC *doc);
void show_kparm(KERNEL_PARM *parm);

#endif
