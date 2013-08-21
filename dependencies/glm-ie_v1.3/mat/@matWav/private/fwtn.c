/*

 FWTN.C	- Periodized Wavelet Transform and its Inverse Transform
   The wavelet transform is an orthonormal transform i.e. norm(x(:))==norm(y(:)) 
   with linear computational complexity O(numel(x)).

   The calling syntax is:

			y = FWTN(x,l,f{,transp});
  
   x        data whose size is l times divisable by 2 i.e. size(x) = sz*2^l
   l        number of levels in the pyramid, 0 means FWTN reduces to y=x;
   f        1d quadrature mirror filter coefficients
   transp   0 for Fast Wavelet Transform and 
            1 for Inverse Fast Wavelet Transform, default is 0 
   y        output data of same size (and norm) as input data x
  
  
  Examples for f corresponding to well-known kinds of wavelets
      [-0.075765714789341, -0.029635527645954,                       (Symmlet,4)
        0.497618667632458,  0.803738751805216, 0.297857795605542,
       -0.099219543576935, -0.012603967262261, 0.032223100604071]
      [ 1+sqrt(3), 3+sqrt(3), 3-sqrt(3), 1-sqrt(3) ]/sqrt(32)    (Daubechies, 4)
      [ 1, 1 ]/sqrt(2)                                                    (Haar)

  Written by Hannes Nickisch, 17 Jan 2008

*/

#include "wavelet.h"
#include "mex.h"
#include <math.h>
#include <string.h>

/* Input Arguments */
#define	SIG_IN	prhs[0] /* X */
#define	LEV_IN	prhs[1]
#define FIL_IN  prhs[2]
#define TRP_IN  prhs[3]

/* Output Arguments */
#define	SIG_OUT	plhs[0] /* Y */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    double *in, *out, *f, *tmp, *dnl;
    int i, j, k, nf, d, *n, nl, max;
    mxArray *temp;
    int transpose=0;

    /* Check for proper number of arguments */
    if (nrhs < 3) {
        mexErrMsgTxt("FWTN requires at least 3 input arguments.");
    } else if (nlhs != 1) {
        mexErrMsgTxt("FWTN requires one output argument.");
    }
    if (nrhs > 3)
        if(mxGetM(TRP_IN)>0)
            transpose = (int)mxGetPr(TRP_IN)[0];
    
    /* Analyze input */
    d  = mxGetNumberOfDimensions(SIG_IN);
    n  = mxMalloc(d*sizeof(int));
    memcpy(n, mxGetDimensions(SIG_IN), d*sizeof(int));
    if (mxGetM(LEV_IN)*mxGetN(LEV_IN)!=1)
        mexErrMsgTxt("FWTN requires a suitable numbers of layers");
    nl = floor((mxGetPr(LEV_IN))[0]+.5);
    if(nl<0) mexErrMsgTxt("FWTN requires a suitable numbers of layers");
    for(i=0; i<d; i++) {
        for(j=1,k=0; j<n[i]; j*=2,k++);
        if(k<nl)    mexErrMsgTxt("FWTN requires a suitable numbers of layers");       
        for(j=n[i],k=0; k<nl; k++,j/=2);
        for(       k=0; k<nl; k++,j*=2);
        if(j!=n[i]) mexErrMsgTxt("FWTN requires appropriate sizes");
    }  
    f   = mxGetPr(FIL_IN);
    nf  = (int) (mxGetM(FIL_IN)*mxGetN(FIL_IN));
    
    /* Create a matrix for the return argument */
    SIG_OUT = mxCreateNumericArray(d, n, mxDOUBLE_CLASS, mxREAL);
    max = n[0]; for(i=1; i<d; i++) if(max<n[i]) max=n[i];
    temp = mxCreateDoubleMatrix(max, 3,  mxREAL);

    /* Assign pointers to the various parameters */
    out = mxGetPr(SIG_OUT);
    tmp = mxGetPr(temp);
    in  = mxGetPr(SIG_IN);
    
    /* Do the actual computations depending on the flag */
    if (transpose)
        iwtn(in, n, d, nl, f, nf, out, tmp);
    else
         wtn(in, n, d, nl, f, nf, out, tmp);
    
    mxDestroyArray(temp);
}
