/*
 * Wrapper for the LBGFS library
 */

#include <stdio.h>
#include <lbfgs.h>
#include "mex.h"

mxArray *rhs[2];

static lbfgsfloatval_t evaluate(
void *instance,
const lbfgsfloatval_t *x,
lbfgsfloatval_t *g,
const int n,
const lbfgsfloatval_t step
) {
    lbfgsfloatval_t fx = 0.0;
    mxArray *lhs[2];
        
    memcpy(mxGetPr(rhs[1]),x,n*sizeof(double));

    /* Call matlab function and return function value and gradient */
    /* rhs contains: (function handle,x) */
    /* lhs contains: (function value,gradient) */
    if (mexCallMATLAB(2,lhs,2,rhs,"feval"))
        mexErrMsgTxt("Error while evaluating objective function.");        
    
    /* Retrieve function value and gradient */
    memcpy(g,mxGetPr(lhs[1]),n*sizeof(double));
    fx = mxGetScalar(lhs[0]); 
    mxDestroyArray(lhs[0]); mxDestroyArray(lhs[1]);
    
    return fx;
}


/*
 * This functions display progress
 */
static int progress(
void *instance,
const lbfgsfloatval_t *x,
const lbfgsfloatval_t *g,
const lbfgsfloatval_t fx,
const lbfgsfloatval_t xnorm,
const lbfgsfloatval_t gnorm,
const lbfgsfloatval_t step,
int n,
int k,
int ls
) {
    mexPrintf("Iteration %d:\n", k);
    mexPrintf("  fx = %f, x[0] = %f\n", fx, x[0]);
    mexPrintf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
    mexPrintf("\n");
    /* Make sure info is printed in "real time" */
    mexEvalString("drawnow;");  
    return 0;
}


/* Gateway routine
 *
 * Inputs: 1. f - function handle
 *         2. x0 - initial guess
 *         3. param - parameter struct
 *
 * Outputs: 1. x   - solution vector
 *          2. msg - output message 
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSize N,ret;
    mxArray *display;
    char *buf;
    mwSize buflen;
    lbfgs_parameter_t param;
    lbfgsfloatval_t *x, fx;
    
    if (nrhs<3)
        mexErrMsgTxt("lbfgs_ requires exactly 3 input arguments.");
    if( mxGetClassID(prhs[0]) != mxFUNCTION_CLASS )    
        mexErrMsgTxt("The first argument must be a function handle.");
    if (!mxIsDouble(prhs[1]))
        mexErrMsgTxt("The second argument must be a column vector.");
    if (!mxIsStruct(prhs[2]))
        mexErrMsgTxt("The third argument must be a struct.");        

    N = mxGetM(prhs[1]);
    
    /* Get function handle for use with mexCallMATLAB*/
    rhs[0] = mxDuplicateArray(prhs[0]);    
    rhs[1] = mxCreateDoubleMatrix(N,1,mxREAL);
    
    /* Allocate memory for solution vector and initialize it to x0 */
    x = lbfgs_malloc(N);
    memcpy(x,mxGetPr(prhs[1]),N*sizeof(double));
    
    /* Get parameters for optimization */
    param.m = (int) mxGetScalar(mxGetField(prhs[2],0,"m"));
    param.epsilon = (lbfgsfloatval_t) mxGetScalar(mxGetField(prhs[2],0,"epsilon"));
    param.past = (int) mxGetScalar(mxGetField(prhs[2],0,"past"));
    param.delta = (lbfgsfloatval_t) mxGetScalar(mxGetField(prhs[2],0,"delta"));
    param.max_iterations = (int) mxGetScalar(mxGetField(prhs[2],0,"MaxIter"));    
    param.linesearch = (int) mxGetScalar(mxGetField(prhs[2],0,"linesearch"));
    param.max_linesearch = (int) mxGetScalar(mxGetField(prhs[2],0,"max_linesearch"));
    param.min_step = (lbfgsfloatval_t) mxGetScalar(mxGetField(prhs[2],0,"min_step"));    
    param.max_step = (lbfgsfloatval_t) mxGetScalar(mxGetField(prhs[2],0,"max_step"));    
    param.ftol = (lbfgsfloatval_t) mxGetScalar(mxGetField(prhs[2],0,"ftol"));        
    param.wolfe = (lbfgsfloatval_t) mxGetScalar(mxGetField(prhs[2],0,"wolfe"));            
    param.gtol = (lbfgsfloatval_t) mxGetScalar(mxGetField(prhs[2],0,"gtol"));                
    param.xtol = (lbfgsfloatval_t) mxGetScalar(mxGetField(prhs[2],0,"xtol")); 
    param.orthantwise_c = (lbfgsfloatval_t) mxGetScalar(mxGetField(prhs[2],0,"orthantwise_c"));                    
    param.orthantwise_start = (int) mxGetScalar(mxGetField(prhs[2],0,"orthantwise_start"));                        
    param.orthantwise_end = (int) mxGetScalar(mxGetField(prhs[2],0,"orthantwise_end"));
    
    display = mxGetField(prhs[2],0,"Display");
    buflen = mxGetN(display)*sizeof(mxChar)+1;
    buf = mxMalloc(buflen);
    if (mxGetString(display, buf, buflen))
        mexErrMsgTxt("Problem with the field 'Display'.");
    
    /* Call solver */
    if (!strcmp(buf,"iter")) /* print info at each iterations */
        ret = lbfgs(N, x, &fx, evaluate, progress, NULL, &param);   
    else if (!strcmp(buf,"final") || !strcmp(buf,"none")) /* no printout */
        ret = lbfgs(N, x, &fx, evaluate, NULL, NULL, &param);
    else
        mexErrMsgTxt("Unknown value for 'options.Display'. Available options are 'iter', 'final' or 'none'.\n");    
    
    /* Allocate outputs */
    plhs[0] = mxCreateDoubleMatrix(N,1,mxREAL);    
    plhs[1] = mxCreateDoubleScalar(fx);
    plhs[2] = mxCreateDoubleScalar(ret);   
    
    /* Copy current iterate to output */
    memcpy(mxGetPr(plhs[0]),x,N*sizeof(double));
    
    /* Release allocated memory */
    lbfgs_free(x);
    mxDestroyArray(rhs[0]);
    mxDestroyArray(rhs[1]);
}
