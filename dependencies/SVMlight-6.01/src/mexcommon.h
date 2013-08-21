#ifndef _MEXCOMMON
#define _MEXCOMMON

#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "svm_common.h"
#include "svm_learn.h"
#include "mexsvmlearn.h"
#include <stdio.h>


#define ERR001 "ERR-001 - Malformed arguments: call as model = mexsvmlearn(samples,labels,options)"
#define ERR002 "ERR-002 - samples must be an mxn matrix, lables must be an mx1, The m dimensions don't match."
#define ERR003 "ERR-003 - labels must be an mx1 matrix"
#define ERR004 "ERR-004 - unknown learning parameter type"
#define ERR005 "ERR-005 - unknown kernel function"
#define ERR006 "ERR-006 - model length is zero."
#define ERR007 "ERR-007 - Malloc failed, out of memory"
#define ERR008 "ERR-008 - SVM Regression not implemented for MEX/MATLAB"
#define ERR009 "ERR-009 - SVM Ranking not implemented for MEX/MATLAB"
#define ERR010 "ERR-010 - Error NULL Pointer given to my_realloc"
#define ERR011 "ERR-011 - Error Invalid point given to my_realloc"

void storeValue(mxArray *mxStruct, char *fieldName, double value);
void storeString(mxArray *mxStruct, char *fieldName, char *value);
void storeArray(mxArray *mxStruct, char *name, double *srcArray, int n);
void storeArrayLong(mxArray *mxStruct, char *name, long *srcArray, int n);
void storeDocs(mxArray *mxStruct, char *name, DOC **docs, int ndocs, int nwords);
void store_kern_parms(mxArray *mxStruct, char *name,KERNEL_PARM *kernel_parm);
void store_model(struct model *model, mxArray *mxOut[] );
void mexToDOC(mxArray *mxData, mxArray *mxLabels, DOC ***docs, double **label,  long int *totwords, long int *totdoc);
int parse_mxEntry(int row, mxArray *mxData, mxArray *mxLabels, WORD *words);
int computeOffset(int numRows, int numCols, int i, int j);
void myBzero(char *ptr, int numbytes);
MODEL *restore_model(const mxArray *);

#endif
