/*=================================================================
 * mexpp.c 
 *
 * starting a locfit interface.
 *
 * $Revision: 1.5 $ */
#include "mex.h"
#include "lfev.h"

extern void lfmxdata(), lfmxevs();

design des;
lfit lf;
int lf_error;

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{ int i, vc, nc, mvc, mg[MXDIM], where, m, wh, dir, nkap;
  const mxArray *dat, *fpc, *sp, *pc, *evs, *iwkc;
  double *y, *x, *fy, *fx[MXDIM], *kap, cv, *fh, *se, *lo, *hi, lev;
  char band[16], wher[16], what[16], rest[16];

  if (nrhs<6) mexErrMsgTxt("ppmex requires 5 inputs.");
  lf_error = 0;

  mut_redirect(mexPrintf);
  dat= mxGetField(prhs[1],0,"data");
  evs= mxGetField(prhs[1],0,"evaluation_structure");
  sp = mxGetField(prhs[1],0,"smoothing_parameters");
  fpc= mxGetField(prhs[1],0,"fit_points");
  pc = mxGetField(prhs[1],0,"parametric_component");
  mxGetString(prhs[2],band,16);
  mxGetString(prhs[3],what,16);
  mxGetString(prhs[4],rest,16);
  dir = mxGetPr(prhs[5])[0];
  kap = mxGetPr(prhs[6]);
  nkap = mxGetN(prhs[6]);
  lev = mxGetPr(prhs[7])[0];
  lfmxdata(&lf.lfd,dat);
  lfmxsp(&lf.sp,sp,lf.lfd.d);
  lfmxevs(&lf,evs); /* this has initialized module */

  if (rest[0]=='n')
    wh = ppwhat(what);
  else
    wh = 64+restyp(rest);
   
    lf.fp.xev = mxGetPr(mxGetField(fpc,0,"evaluation_points"));
    lf.lfd.d = lf.fp.d = mxGetM(mxGetField(fpc,0,"evaluation_points"));
    lf.fp.nv = lf.fp.nvm = mxGetN(mxGetField(fpc,0,"evaluation_points"));
    lf.fp.wk = mxGetPr(mxGetField(fpc,0,"fitted_values"));
    lf.fp.h = NULL;
    if (lf.mdl.alloc!=NULL) lf.mdl.alloc(&lf);
        
    iwkc = mxGetField(fpc,0,"evaluation_vectors");
    vc = mxGetM(mxGetField(iwkc,0,"cell"));
    nc = mxGetN(mxGetField(iwkc,0,"cell"));
    mvc = mxGetN(mxGetField(iwkc,0,"splitvar"));
    lf.evs.iwk = mxCalloc(vc*nc+3*mvc,sizeof(int));
    lf.evs.nce = lf.evs.ncm = nc;
    y = mxGetPr(mxGetField(iwkc,0,"cell"));
    lf.evs.ce = lf.evs.iwk;
    for (i=0; i<vc*nc; i++) lf.evs.ce[i] = y[i];
    y = mxGetPr(mxGetField(iwkc,0,"splitvar"));
    lf.evs.s = &lf.evs.iwk[vc*nc];
    for (i=0; i<mvc; i++) lf.evs.s[i] = y[i];
    y = mxGetPr(mxGetField(iwkc,0,"lo"));
    lf.evs.lo = &lf.evs.s[mvc];
    for (i=0; i<mvc; i++) lf.evs.lo[i] = y[i];
    y = mxGetPr(mxGetField(iwkc,0,"hi"));
    lf.evs.hi = &lf.evs.lo[mvc];
    for (i=0; i<mvc; i++) lf.evs.hi[i] = y[i];
    lf.fp.hasd = deg(&lf.sp)>0;

    lf.fp.kap = mxGetPr(mxGetField(fpc,0,"kappa"));

    lf.pc.wk = mxGetPr(pc);
    lf.pc.lwk = mxGetN(pc);
    pcchk(&lf.pc,lf.fp.d,npar(&lf.sp),1);
    haspc(&lf.pc) = !noparcomp(&lf.sp);
    lf.pc.xtwx.st = JAC_EIGD;
    
  where = 0;
  if (mxIsCell(prhs[0])) /* interpret as grid margins */
  { m = 1;
    for (i=0; i<lf.fp.d; i++)
    { fx[i] = mxGetPr(mxGetCell(prhs[0],i));
      mg[i] = mxGetM(mxGetCell(prhs[0],i));
      m *= mg[i];
    }
    where = 2;
  }
  if (mxIsChar(prhs[0]))
  { mxGetString(prhs[0],wher,16);
    where = 3; /* data */
    m = mg[0] = lf.lfd.n;
    if (wher[0] == 'f') /* or fit points */
    { where = 4;
      m = mg[0] = lf.fp.nv;
    }
  }
  if (where==0)                   /* interpret as numeric vector */
  { x = mxGetPr(prhs[0]);
    m = mg[0] = mxGetM(prhs[0]); /* should check mxGetN == d */
    for (i=0; i<lf.lfd.d; i++)
      fx[i] = &x[i*m];
    where=1;
  }

  plhs[0] = mxCreateDoubleMatrix(m,1,mxREAL); /* fitted values */
  plhs[1] = mxCreateDoubleMatrix(m,1,mxREAL); /* std errors */
  fh = mxGetPr(plhs[0]);
  se = mxGetPr(plhs[1]);

  preplot(&lf,fx,fh,se,band[0],mg,where,wh,dir);

  cv = 0.0;
  if (band[0] != 'n')
  { cv = critval(1-lev,kap,nkap,lf.lfd.d,TWO_SIDED,0.0,GAUSS);
    plhs[2] = mxCreateDoubleMatrix(m,2,mxREAL);
    lo = mxGetPr(plhs[2]);
    hi = &lo[m];
    for (i=0; i<m; i++)
    { lo[i] = fh[i]-cv*se[i]; 
      hi[i] = fh[i]+cv*se[i]; 
    }
  }
  else plhs[2] = mxCreateDoubleScalar(0.0);
}
