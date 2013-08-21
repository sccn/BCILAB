/*
 * =============================================================
 * mexsvmlearn - link MATLAB with SVM-Lite by Thorstenn Joachims.
 * This is a MEX-file for MATLAB.  
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

  
  DOC **docs; /* hold a test example */
  double *target; /* hold labels */
  WORD *words;  /* the words read from the example */
  long rows, cols; /* the number of rows and cols in the test data */
  double dist, doc_label, costfactor;
  double *err,*pred;
  long int correct=0, incorrect=0, none=0,i;
  MODEL *model;
  checkParameters(nlhs, plhs, nrhs, prhs);

  global_init( );

  /* load model parameters from the "model" parameter */
  model = restore_model((mxArray *)prhs[2]);

  rows = mxGetM(prhs[0]);
  cols = mxGetN(prhs[0]);

  /* load the testing arrays into docs */

  mexToDOC((mxArray *)prhs[0], (mxArray *)prhs[1], &docs, &target, NULL, NULL);
  
  /* setup output environment */
  plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
  plhs[1] = mxCreateDoubleMatrix(rows,1,mxREAL);

  err = mxGetPr(plhs[0]);
  pred = mxGetPr(plhs[1]);

  /* classify examples */
  for (i = 0; i < rows; i++) {

    dist = classify_example(model, docs[i]);
    pred[i] = dist;

    if (dist > 0) {
      if (target[i] > 0) correct++;
      else incorrect++;
    } else {
      if (target[i] < 0) correct++;
      else incorrect++;
    }

    if ((int)(0.1+(target[i] * target[i]))!=1)
      none++;

  }

  err[0] = incorrect / (double) rows;

  
  global_destroy( );

}


void checkParameters(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  if (nlhs > 2)
    mexErrMsgTxt("Error: call only with two parameters: [ err, predictions ]");

  if (nrhs != 3)
    mexErrMsgTxt("Error: call with three input parameters mexsvmclassify(x,y,model)");

  if ((!mxIsDouble(prhs[0])) && (mxIsDouble(prhs[1])))
    mexErrMsgTxt("Error: x and y must be double arrays ");

  if (mxGetM(prhs[0]) != mxGetM(prhs[1]))
    mexErrMsgTxt("Error: x and y must have the same number of rows");

  if (!mxIsStruct(prhs[2]))
    mexErrMsgTxt("Error: model is not a struct (output from mexsvmlearn)");


}
