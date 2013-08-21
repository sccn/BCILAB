/***********************************************************************/
/*								       */
/*   global.c							       */
/*								       */
/*   Global declarations used throughout the svm package	       */
/*								       */
/*   Author: Tom Briggs						       */
/*   For use with SVM Lite by Thorsten Joachims			       */
/*   Date: January 1, 2005					       */
/*								       */
/*   Copyright (c) 2002	 Thorsten Joachims - All rights reserved       */
/*								       */
/*   This software is available for non-commercial use only. It must   */
/*   not be modified and distributed without prior permission of the   */
/*   author. The author is not responsible for implications from the   */
/*   use of this software.					       */
/*								       */
/***********************************************************************/

#include "svm_common.h"
#include "global.h"
#include "mem_clean.h"


long   verbosity;	       /* verbosity level (0-4)  */
long   kernel_cache_statistic; 

/* these are necessary to support the memory 
 * cleaner - and my be defined outside of the
 * my_malloc() functions */
#ifdef MATLAB_MEX

#include <mex.h>
#include <matrix.h>

#define MALLOC mxMalloc
#define FREE mxFree

#else

#include <malloc.h>

#define MALLOC malloc
#define FREE free
#endif

/* from svm_hideo.c */
double *primal=0,*dual=0; 
long   precision_violations=0; 
double opt_precision=DEF_PRECISION; 
long   maxiter=DEF_MAX_ITERATIONS;
double lindep_sensitivity=DEF_LINDEP_SENSITIVITY;
double *buffer; 
long   *nonoptimal; 

long  smallroundcount=0;
long  roundnumber=0;

/* malloc array */
HASH_ARRAY *malloc_array;

/* ***********************************
 * global initialization routine - set all
 * globals to known values 
 ************************************ */

void global_init( )
{
  verbosity = 0;
  kernel_cache_statistic = 0;

  /* global variables from svm_hideo.c */
  primal=0;
  dual=0;
  precision_violations=0;
  opt_precision=DEF_PRECISION;
  maxiter=DEF_MAX_ITERATIONS;
  lindep_sensitivity=DEF_LINDEP_SENSITIVITY;
  smallroundcount=0;
  roundnumber=0;
  
  /* ********************************* */
  /* build the malloc array */
  malloc_array = MALLOC(sizeof(HASH_ARRAY));
  hash_init_array(malloc_array);
}


/*********************************
 * global destructor - called at program
 * exit to clear out any allocated memory
 * and purge any global values.
 ********************************** */
void global_destroy( )
{
  hash_destroy_array(malloc_array);
  FREE(malloc_array);
}


/*********************************
 * For debugging - display the parameters
 * of a given model
 ********************************* */
void show_model(MODEL *model)
{
  int i;

  printf("------------------------------\n");
  printf("sv_num %ld\n", model->sv_num);
  printf("at_ub %ld\n", model->at_upper_bound);
  printf("b %f", model->b);

  printf("alphas\n");
  if (model->alpha != NULL)
    {
      for (i = 0; i < model->totdoc; i++)
	printf("%d  %f\n",i,  model->alpha[i]);
    }

  printf("indices:\n");
  for (i = 0; i < model->totdoc; i++)
    printf("%d %ld\n", i, model->index[i]);

  show_kparm(&model->kernel_parm);

  for (i = 0; i < model->sv_num; i++)
    show_doc(model->supvec[i]);

  printf("------------------------------\n");


}

/* ***********************************************
 * Display the given document 
 * *********************************************** */
void show_doc(DOC *doc)
{
  if (doc) {
    printf("   docnum: %ld  queryid: %ld  costfactor: %f slackid: %ld\n", 
	   doc->docnum, doc->queryid, doc->costfactor, doc->slackid);
  }

}


/* **********************************************
 * Display the kernel parameters   
 * ********************************************** */
void show_kparm(KERNEL_PARM *parm)
{
  printf("----------| START Kernel Parm |--------------\n");
  printf("  type: %ld\n", parm->kernel_type);
  printf("  degree: %ld\n", parm->poly_degree);
  printf("  rbf_gamma: %f\n", parm->rbf_gamma);
  printf("  coef lin: %f\n", parm->coef_lin);
  printf("  coef_const: %f\n", parm->coef_const);
  printf("  custom: %s\n", parm->custom);
  printf("----------| END Kernel Parm |--------------\n");


}
