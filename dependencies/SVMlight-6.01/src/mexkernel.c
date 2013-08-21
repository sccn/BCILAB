/*
 * =============================================================
 * mexsvmlearn - link MATLAB with SVM-Lite by Thorstenn Joachims.
 * This is a MEX-file for MATLAB.   
 *
 * Compute the entire kernel matrix for the given vectors and kernel
 *
 * MEX-Interface:
 *    Thomas Briggs (c) 2004
 * SVM-Lite:
 *    Thorsten Joachims
 * =============================================================
 */
#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "svm_common.h"
#include "svm_learn.h"
#include "mexsvmlearn.h"
#include "mexcommon.h"
#include "global.h"



void checkParameters(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

/* call as  model = mexsvmlearn(data,labels,options) */
void mexFunction(int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[])
{
  int rows, cols,i,j,offset;
  double *kern,k;
  DOC **docs;
  SVECTOR *a, *b;
  MODEL *model;

  global_init();

  /* load model parameters from the "model" parameter */
  model = restore_model((mxArray *)prhs[1]);

  rows = mxGetM(prhs[0]);
  cols = mxGetN(prhs[0]);

  /* load the testing arrays into docs */
  mexToDOC((mxArray *)prhs[0], NULL, &docs, NULL, NULL, NULL);
  
  /* setup output environment */
  plhs[0] = mxCreateDoubleMatrix(rows,rows,mxREAL);
  kern = mxGetPr(plhs[0]);

  for (i = 0; i < rows; i++) {
    a = docs[i]->fvec;
    for (j = 0; j < rows; j++) {
      b = docs[j]->fvec;
      k = single_kernel(&(model->kernel_parm), a, b);

      offset = computeOffset(rows, rows, i, j);      
      kern[offset] = k;
    }
  }

  global_destroy();
}


void checkParameters(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  if (nlhs != 1)
    mexErrMsgTxt("Error: call only with two parameters: [ err, predictions ]");

  if (nrhs != 2)
    mexErrMsgTxt("Error: call with three input parameters mexsvmclassify(x,y,model)");

  if (!mxIsDouble(prhs[0]))
    mexErrMsgTxt("Error: x and y must be double arrays ");

  if (!mxIsStruct(prhs[1]))
    mexErrMsgTxt("Error: model is not a struct (output from mexsvmlearn)");


}
