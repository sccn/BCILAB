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

void checkParameters(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

/* call as  model = mexsvmlearn(data,labels,options) */
void mexFunction(int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[])
{
  int rows, cols,i,j,offset;
  double *kern,k;
  DOC **docsA, **docsB;
  SVECTOR *a, *b;
  MODEL *model;

  global_init( );

  /* load model parameters from the "model" parameter */
  model = restore_model((mxArray *)prhs[2]);

  rows = mxGetM(prhs[0]);
  cols = mxGetN(prhs[0]);

  /* load the testing arrays into docs */
  mexToDOC((mxArray *)prhs[0], NULL, &docsA, NULL, NULL, NULL);
  mexToDOC((mxArray *)prhs[1], NULL, &docsB, NULL, NULL, NULL);
  
  /* setup output environment */
  plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
  kern = mxGetPr(plhs[0]);

  a = docsA[0]->fvec;
  b = docsB[0]->fvec;

  kern[0] = single_kernel(&(model->kernel_parm), a, b);


  global_destroy();
}


void checkParameters(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  if (nlhs != 1)
    mexErrMsgTxt("Error: call only with two parameters: [ err, predictions ]");

  if (nrhs != 3)
    mexErrMsgTxt("Error: call with three input parameters mexsvmclassify(x,y,model)");

  if ((!mxIsDouble(prhs[0])) && (!mxIsDouble(prhs[1])))
    mexErrMsgTxt("Error: x and y must be double arrays ");

  if ((mxGetM(prhs[0]) != 1) || (mxGetM(prhs[1])!=1))
    mexErrMsgTxt("Error: x and y must be row vectors.");

  if (!mxIsStruct(prhs[2]))
    mexErrMsgTxt("Error: model is not a struct (output from mexsvmlearn)");


}
