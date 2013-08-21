/*
 *  mex_interface.h
 *  mex_svmperf
 *
 *  Created by Oscar Luaces on 22/02/08.
 *  Copyright 2008 AIC-University of Oviedo at Gijón. All rights reserved.
 *
 */

#ifndef _MEX_INTERFACE_
#define _MEX_INTERFACE_

#include <mex.h>

#define printf	mexPrintf
#define perror	mexErrMsgTxt
#define exit(n) mexErrMsgTxt("EXITING...")
#define malloc	mxMalloc
#define realloc my_realloc
#define free	my_free
#define fflush(x)	nop()
#define fprintf	my_fprintf

#define MAX_ARGVS	128

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void create_argc_argv(const mxArray *prhs,int *argc,char **argv[]);
void my_fprintf(FILE *, const char *, ...); 
void my_free(void *ptr); 
void *my_malloc(size_t size);
void *my_realloc(void *ptr, size_t size); 

void nop(); 

#endif