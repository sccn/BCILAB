/*=================================================================
 * mexlocfit.c 
 *
 * starting a locfit interface.
 *
/* $Revision: 1.5 $ */
#include "mex.h"
#include "lfev.h"

design des;
lfit lf;
int lf_error;

extern void lfmxdata(), lfmxsp(), lfmxevs();

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    int  d, i, nvc[5], nvm, vc, nv, nc, mvc, lw[7];
    double *y;
    mxArray *iwkc, *temparray;
    const char * fpcnames[] = {"evaluation_points","fitted_values","evaluation_vectors","fit_limits","family_link","kappa"};
    const char * evsnames[] = {"cell","splitvar","lo","hi"};

    mut_redirect(mexPrintf);
    if (nrhs != 3) mexErrMsgTxt("mexlf requires 3 inputs.");
    lf_error = 0;

    lfit_alloc(&lf);
    lfit_init(&lf);

    lfmxdata(&lf.lfd,prhs[0]);
    d = lf.lfd.d;

    lfmxsp(&lf.sp,prhs[2],d);
    lfmxevs(&lf,prhs[1]);

    guessnv(&lf.evs,&lf.sp,&lf.mdl,lf.lfd.n,d,lw,nvc);
    nvm = lf.fp.nvm = nvc[0];
    vc = nvc[2];
    if (ev(&lf.evs) != EPRES) lf.fp.xev = mxCalloc(nvm*d,sizeof(double));
    lf.fp.lev = nvm*d;
    lf.fp.wk = lf.fp.coef = mxCalloc(lw[1],sizeof(double));
    lf.fp.lwk = lw[1];
    lf.evs.iwk = mxCalloc(lw[2],sizeof(int));
    lf.evs.liw = lw[2];
    plhs[1] = mxCreateDoubleMatrix(1,lw[3],mxREAL);
    lf.pc.wk = mxGetPr(plhs[1]);
    lf.pc.lwk = lw[3];
    lf.fp.kap = mxCalloc(lw[5],sizeof(double));
/* should also allocate design here */
    
    if (!lf_error) startmodule(&lf,&des,NULL,NULL);

/* now, store the results: 
   plhs[0] stores informtion about fit points and evaluation structure.
    it is now a matlab structure, not a cell   
   fit_points.evaluation_points - matrix of evaluation points.
   fit_points.fitted_values - matrix fitted values etc.
   fit_points.evaluation_vectors - structure of 'integer' vectors {cell,splitvar,lo,hi}
   fit_points.fit_limits - fit limit (matrix with 2 cols).
   fit_points.family_link - [familt link] numeric vector.
   fit_points.kappa - kap vector.
*/
    plhs[0] = mxCreateStructMatrix(1,1,6,fpcnames);
    if ( plhs[0] == NULL ) {
       printf("Problem with CreateStructMatrix for plhs[0]\n");fflush(stdout);
    }

    mxSetField(plhs[0],0,"evaluation_points",mxCreateDoubleMatrix(d,lf.fp.nv,mxREAL));
    memcpy(mxGetPr(mxGetField(plhs[0],0,"evaluation_points")), lf.fp.xev, d*lf.fp.nv*sizeof(double));

    mxSetField(plhs[0],0,"fitted_values",mxCreateDoubleMatrix(lf.fp.nv,lf.mdl.keepv,mxREAL));
    for (i=0; i<lf.mdl.keepv; i++)
      memcpy(&mxGetPr(mxGetField(plhs[0],0,"fitted_values"))[i*lf.fp.nv], &lf.fp.coef[i*nvm], lf.fp.nv*sizeof(double));
    /* another bit to save here? -- split values, kdtree */
   
    temparray = mxCreateStructMatrix(1,1,4,evsnames);
    if ( temparray == NULL ) {
      printf("Problem with CreateStructMatrix for temparray\n");fflush(stdout);
    }
    mxSetField(plhs[0],0,"evaluation_vectors",temparray);
    iwkc = mxGetField(plhs[0],0,"evaluation_vectors");
    nv = lf.fp.nv;
    nc = lf.evs.nce;
    mvc = (nv>nc) ? nv : nc;
    mxSetField(iwkc,0,"cell",mxCreateDoubleMatrix(vc,nc,mxREAL));  /* ce */
    mxSetField(iwkc,0,"splitvar",mxCreateDoubleMatrix(1,mvc,mxREAL));  /* s  */
    mxSetField(iwkc,0,"lo",mxCreateDoubleMatrix(1,mvc,mxREAL));  /* lo */
    mxSetField(iwkc,0,"hi",mxCreateDoubleMatrix(1,mvc,mxREAL));  /* hi */
    y = mxGetPr(mxGetField(iwkc,0,"cell"));
    for (i=0; i<vc*nc; i++) y[i] = lf.evs.ce[i];
    y = mxGetPr(mxGetField(iwkc,0,"splitvar"));
    for (i=0; i<mvc; i++) y[i] = lf.evs.s[i];
    y = mxGetPr(mxGetField(iwkc,0,"lo"));
    for (i=0; i<mvc; i++) y[i] = lf.evs.lo[i];
    y = mxGetPr(mxGetField(iwkc,0,"hi"));
    for (i=0; i<mvc; i++) y[i] = lf.evs.hi[i];

    mxSetField(plhs[0],0,"fit_limits",mxCreateDoubleMatrix(d,2,mxREAL));
    memcpy(mxGetPr(mxGetField(plhs[0],0,"fit_limits")), lf.evs.fl, 2*d*sizeof(double));
    
    mxSetField(plhs[0],0,"family_link",mxCreateDoubleMatrix(1,2,mxREAL));
    y = mxGetPr(mxGetField(plhs[0],0,"family_link"));
    y[0] = fam(&lf.sp);
    y[1] = link(&lf.sp);

    mxSetField(plhs[0],0,"kappa",mxCreateDoubleMatrix(1,lf.mdl.keepc,mxREAL));
    memcpy(mxGetPr(mxGetField(plhs[0],0,"kappa")),lf.fp.kap,lf.mdl.keepc*sizeof(double));
}
