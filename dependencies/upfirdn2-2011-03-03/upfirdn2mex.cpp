#include "mex.h"
#include "Resampler.h"


/* Upsample by P (inserting zeroes), apply given FIR filter H, and downsample by Q.
 * [Y,Zf,Tf,Of] = upfirdn2mex(X,H,P,Q,Zi,Ti,Oi)
 *
 * Similar to upfirdn, except that final filter conditions are being returned, which can
 * be used as initial conditions of a subsequent run of upfirdn2mex. 
 * Also, X must be real and a vector.
 *
 * In:
 *  X : the signal to be resampled; real-valued column vector
 *
 *  H : the filter kernel to be used; column vector
 *
 *  P : upsampling factor, must be a positive integer
 *
 *  Q : downsampling factor, must be a positive integer
 *
 *  Zi : optional initial filter conditions
 *
 *  Ti : optional initial fractional time
 *
 *  Oi : optional initial fractional offset
 *
 * Out:
 *  Y : resampled singnal; column vector (note: Y may be delayed depending on the
 *      delay of the FIR filter H, and the length of Y may depend the initial conditions)
 *
 *  Zf : final conditions of the filter
 *
 *  Tf : final filter fractional time
 *
 *  Of : final filter fractional offset
 *
 *                          Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
 *                          2011-03-03
 */

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]) {
    int Lx, Lh, Ly, Ls, k; /* length of X, length of H, length of Y, length of Zi/Zf */
    double p, q, Ti, Oi;
    void *X,*H,*Zi;    
    double *ptrd;
    float *ptrf;
    
    /* check input argument count */
    if (!((nrhs == 4) || (nrhs == 7)))
        mexErrMsgTxt("Either 4 or 7 input arguments may be specified.");
    
    /* check output argument count */
    if (!((nlhs == 1) || (nlhs == 4)))
        mexErrMsgTxt("Either 1 or 4 output arguments may be specified.");
    
    /* Parse X */
    Lx = mxGetM(prhs[0]);
    if (mxGetN(prhs[0]) > 1)
        mexErrMsgTxt("X must be a column vector.");
    if (!(mxIsDouble(prhs[0]) || mxIsSingle(prhs[0])))
        mexErrMsgTxt("X must be either of type double or single.");
    if (mxIsComplex(prhs[0]))
        mexErrMsgTxt("X must not be complex-valued.");
    X = mxGetData(prhs[0]);
    
    /* Parse H */
    Lh = mxGetM(prhs[1]);
    if (mxGetN(prhs[1]) > 1)
        mexErrMsgTxt("H must be a column vector.");
    if (!(mxIsDouble(prhs[1]) || mxIsSingle(prhs[1])))
        mexErrMsgTxt("H must be either of type double or single.");
    if (mxIsComplex(prhs[1]))
        mexErrMsgTxt("H must not be complex-valued.");
    H = mxGetData(prhs[1]);
    
    /* Parse P */
    if (mxGetNumberOfElements(prhs[2]) > 1)
        mexErrMsgTxt("P must be a scalar.");
    if (!mxIsDouble(prhs[2]))
        mexErrMsgTxt("P must be of type double.");
    p = mxGetScalar(prhs[2]);
    if ((p - ceil(p)) > 0 || p < 1)
        mexErrMsgTxt("P must a positive integer.");
    
    /* Parse Q */
    if (mxGetNumberOfElements(prhs[3]) > 1)
        mexErrMsgTxt("Q must be a scalar.");
    if (!mxIsDouble(prhs[3]))
        mexErrMsgTxt("Q must be of type double.");
    q = mxGetScalar(prhs[3]);
    if ((q - ceil(q)) > 0 || q < 1)
        mexErrMsgTxt("Q must a positive integer.");
    
    /* Parse Zi, Ti, Oi, if applicable */
    if (nrhs == 7) {
        if (mxGetClassID(prhs[4]) != mxGetClassID(prhs[0]))
            mexErrMsgTxt("Zi must be of the same type as X.");
        if (mxIsComplex(prhs[4]))
            mexErrMsgTxt("Zi cannot be complex-valued.");
        if (mxGetN(prhs[4]) > 1)
            mexErrMsgTxt("Zi must be a column vector.");
        Zi = mxGetData(prhs[4]);
        
        if (mxGetNumberOfElements(prhs[5]) > 1)
            mexErrMsgTxt("Ti must be a scalar.");
        if (!mxIsDouble(prhs[5]))
            mexErrMsgTxt("Ti must be of type double.");
        Ti = mxGetScalar(prhs[5]);
        if ((Ti - ceil(Ti)) != 0)
            mexErrMsgTxt("Ti must an integer.");
        
        if (mxGetNumberOfElements(prhs[6]) > 1)
            mexErrMsgTxt("Oi must be a scalar.");
        if (!mxIsDouble(prhs[6]))
            mexErrMsgTxt("Oi must be of type double.");
        Oi = mxGetScalar(prhs[6]);
        if ((Oi - ceil(Oi)) != 0)
            mexErrMsgTxt("Oi must an integer.");
    } else {
        /* otherwise all null's (or null ptrs)*/
        Zi = 0;
        Ti = 0;
        Oi = 0;
    }
    
    if (mxIsDouble(prhs[0])) {
        if (mxIsDouble(prhs[1])) {           
            /* double X, double H */
            Resampler<double,double,double> res((int)p, (int)q, (double*)H, Lh, (double*)Zi, (int)Ti, (int)Oi);
            /* create output array... */
            Ly = res.neededOutCount(Lx);
            plhs[0] = mxCreateDoubleMatrix(Ly, 1, mxREAL);
            /* apply resampler */
            res.apply((double*)X,Lx,(double*)mxGetData(plhs[0]),Ly);
            //mexPrintf("sdfdfsdf %i %i\n",res._t,res._xOffset);
            /* create final conditions array, if applicable */
            if (nlhs == 4) {
                Ls = res.coefsPerPhase();
                plhs[1] = mxCreateDoubleMatrix(Ls, 1, mxREAL);
                ptrd = (double*)mxGetData(plhs[1]);
                for (k=0;k<Ls;k++)
                    ptrd[k] = res.finalConds()[k];
                plhs[2] = mxCreateDoubleScalar(res.get_t());
                plhs[3] = mxCreateDoubleScalar(res.get_xoff());
            }
        } else {
            /* double X, single H */
            Resampler<double,double,float> res((int)p, (int)q, (float*)H, Lh, (double*)Zi, (int)Ti, (int)Oi);
            /* create output array... */
            Ly = res.neededOutCount(Lx);
            plhs[0] = mxCreateDoubleMatrix(Ly, 1, mxREAL);
            /* apply resampler */
            res.apply((double*)X,Lx,(double*)mxGetData(plhs[0]),Ly);
            /* create final conditions array, if applicable */
            if (nlhs == 4) {
                Ls = res.coefsPerPhase();
                plhs[1] = mxCreateDoubleMatrix(Ls, 1, mxREAL);
                ptrd = (double*)mxGetData(plhs[1]);
                for (k=0;k<Ls;k++)
                    ptrd[k] = res.finalConds()[k];
                plhs[2] = mxCreateDoubleScalar(res.get_t());
                plhs[3] = mxCreateDoubleScalar(res.get_xoff());
            }            
        }
    } else {
        if (mxIsDouble(prhs[1])) {
            /* single X, double H */
            Resampler<float,float,double> res((int)p, (int)q, (double*)H, Lh, (float*)Zi, (int)Ti, (int)Oi);
            /* create output array... */
            Ly = res.neededOutCount(Lx);
            plhs[0] = mxCreateNumericMatrix(Ly, 1, mxSINGLE_CLASS, mxREAL);
            /* apply resampler */
            res.apply((float*)X,Lx,(float*)mxGetData(plhs[0]),Ly);
            /* create final conditions array, if applicable */
            if (nlhs == 4) {
                Ls = res.coefsPerPhase();
                plhs[1] = mxCreateNumericMatrix(Ls, 1, mxSINGLE_CLASS, mxREAL);
                ptrf = (float*)mxGetData(plhs[1]);
                for (k=0;k<Ls;k++)
                    ptrf[k] = res.finalConds()[k];
                plhs[2] = mxCreateDoubleScalar(res.get_t());
                plhs[3] = mxCreateDoubleScalar(res.get_xoff());
            }            
        } else {
            /* single X, single H */
            Resampler<float,float,float> res((int)p, (int)q, (float*)H, Lh, (float*)Zi, (int)Ti, (int)Oi);
            /* create output array... */
            Ly = res.neededOutCount(Lx);
            plhs[0] = mxCreateNumericMatrix(Ly, 1, mxSINGLE_CLASS, mxREAL);
            /* apply resampler */
            res.apply((float*)X,Lx,(float*)mxGetData(plhs[0]),Ly);            
            /* create final conditions array, if applicable */
            if (nlhs == 4) {
                Ls = res.coefsPerPhase();
                plhs[1] = mxCreateNumericMatrix(Ls, 1, mxSINGLE_CLASS, mxREAL);
                ptrf = (float*)mxGetData(plhs[1]);
                for (k=0;k<Ls;k++)
                    ptrf[k] = res.finalConds()[k];
                plhs[2] = mxCreateDoubleScalar(res.get_t());
                plhs[3] = mxCreateDoubleScalar(res.get_xoff());
            }            
        }
    }
}
