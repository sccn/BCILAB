#ifndef _MEXSVMLEARN
#define _MEXSVMLEARN



#include "mex.h"
#include "matrix.h"
#include "math.h"

char **make_argv(mxArray *, int *);

void read_input_parameters(int,char **,char *,char *, char *,long *,
			   LEARN_PARM *,KERNEL_PARM *);

void mexToDOC(mxArray *, mxArray *, DOC ***, 
	      double **, long int *, long int *);

int parse_mxEntry(int, mxArray *, mxArray *, WORD *);
void   read_input_parameters(int, char **, char *, char *, char *, long *, 
			     LEARN_PARM *, KERNEL_PARM *);
void wait_any_key();
void store_model(struct model *model, mxArray *mxOut[] );
void myBzero(char *ptr, int numbytes);
void extract_user_opts(mxArray *options, KERNEL_PARM *kernel_parm);
void freeArgv(int argc, char **argv);
void freeModelMem(MODEL *model);
#endif
