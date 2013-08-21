#include "mex.h"
#include "lfev.h"

void lfmxdata(lfdata *lfd, const mxArray *data)
{ const mxArray *tmp;
  double *x, *sc, *xl;
  int n, d, i, k;
  char c, ststr[16];
  
  n = lfd->n = mxGetM(mxGetField(data,0,"x"));
  d = lfd->d = mxGetN(mxGetField(data,0,"x"));
  x = mxGetPr(mxGetField(data,0,"x"));
  for (i=0; i<d; i++) lfd->x[i] = &x[i*n];
  
  lfd->y = mxGetPr(mxGetField(data,0,"y"));
  
  tmp = mxGetField(data,0,"weights");
  k = mxGetM(tmp);
  lfd->w = (k>1) ? mxGetPr(tmp) : NULL;

  tmp = mxGetField(data,0,"censor");
  k = mxGetM(tmp); 
  lfd->c = (k>1) ? mxGetPr(tmp) : NULL;

  tmp = mxGetField(data,0,"baseline");
  k = mxGetM(tmp);
  lfd->b = (k>1) ? mxGetPr(tmp) : NULL;

  mxGetString(mxGetField(data,0,"style"),ststr,16);
  k = strlen(ststr);
  for (i=0; i<d; i++)
  { c = (k==1) ? ststr[0] : ststr[i];
    switch(c)
    { case 'a': lfd->sty[i] = STANGL; break;
      case 'l': lfd->sty[i] = STLEFT; break;
      case 'r': lfd->sty[i] = STRIGH; break;
      case 'c': lfd->sty[i] = STCPAR; break;
      default:  lfd->sty[i] = 1; break;
    }
  }

  sc = mxGetPr(mxGetField(data,0,"scales"));
  k = mxGetN(mxGetField(data,0,"scales"));
  for (i=0; i<d; i++) lfd->sca[i] = (k==1) ? sc[0] : sc[i];

  xl = mxGetPr(mxGetField(data,0,"xlim"));
  for (i=0; i<d; i++)
  { lfd->xl[i] = xl[2*i];
    lfd->xl[i+d] = xl[2*i+1];
  }
}

void lfmxsp(smpar *sp, mxArray *mcell, int d)
{ double *alpha;
  char str[16];

  alpha = mxGetPr(mxGetField(mcell,0,"alpha"));
  nn(sp)  = alpha[0];
  fixh(sp)= alpha[1];
  pen(sp) = alpha[2];

  mxGetString(mxGetField(mcell,0,"adaptive_criterion"),str,16);
  acri(sp) = lfacri(str);

  deg(sp) = mxGetPr(mxGetField(mcell,0,"degree"))[0];
  deg0(sp) = -1;

  mxGetString(mxGetField(mcell,0,"family"),str,16);
  fam(sp) = lffamily(str);
  mxGetString(mxGetField(mcell,0,"link"),str,16);
  link(sp) = lflink(str);
  setfamily(sp);

  mxGetString(mxGetField(mcell,0,"kernel"),str,16);
  ker(sp) = lfkernel(str);
  mxGetString(mxGetField(mcell,0,"kernel_type"),str,16);
  kt(sp) = lfketype(str);
  npar(sp) = calcp(sp,d);

  de_renorm = (int)(mxGetPr(mxGetField(mcell,0,"deren"))[0]);
  mxGetString(mxGetField(mcell,0,"deit"),str,16);
  de_itype = deitype(str);
  de_mint = (int)(mxGetPr(mxGetField(mcell,0,"demint"))[0]);
  lf_debug = (int)(mxGetPr(mxGetField(mcell,0,"debug"))[0]);
}

void lfmxevs(lfit *lf, mxArray *mcell)
{ int d, i, j;
  double *ll, *ur, *mg, *drv;
  char evstr[16], mod[16], mdir[256];
  evstruc *evs;
  fitpt *fp;
  deriv *dv;

  evs = &lf->evs;
  fp  = &lf->fp;
  dv  = &lf->dv;
  d   = lf->lfd.d;

  if (mxIsChar(mxGetField(mcell,0,"type")))
  { mxGetString(mxGetField(mcell,0,"type"),evstr,16);
    ev(evs) = lfevstr(evstr);
  }
  else
  { ev(evs) = EPRES;
    evs->mg[0] = mxGetN(mxGetField(mcell,0,"type"));
    fp->xev = mxGetPr(mxGetField(mcell,0,"type"));
  }


  mxGetString(mxGetField(mxGetField(mcell,0,"module"),0,"name"),mod,16);
  mxGetString(mxGetField(mxGetField(mcell,0,"module"),0,"directory"),mdir,256);
  MODPARAMS(lf) = mxGetPr(mxGetField(mxGetField(mcell,0,"module"),0,"parameters"));
  MODNPARAMS(lf) = mxGetN(mxGetField(mxGetField(mcell,0,"module"),0,"parameters"));
  initmodule(&lf->mdl,mod,mdir,lf);


  ll = mxGetPr(mxGetField(mcell,0,"lower_left"));
  ur = mxGetPr(mxGetField(mcell,0,"upper_right"));
  mg = mxGetPr(mxGetField(mcell,0,"grid"));
  j =  mxGetN(mxGetField(mcell,0,"grid"));
  cut(evs) = mxGetPr(mxGetField(mcell,0,"cut"))[0];
  for (i=0; i<d; i++)
  { evs->fl[i] = ll[i];
    evs->fl[i+d] = ur[i];
    if (ev(evs) != EPRES) evs->mg[i] = (j==1) ? mg[0] : mg[i];
  }
  mk(evs) = mxGetPr(mxGetField(mcell,0,"maxk"))[0];

  drv = mxGetPr(mxGetField(mcell,0,"derivative")); 
  j = mxGetN(mxGetField(mcell,0,"derivative")); 
  for (i=0; i<j; i++) dv->deriv[i] = drv[i]-1;
  dv->nd = (drv[0]>0) ? j : 0;
}
