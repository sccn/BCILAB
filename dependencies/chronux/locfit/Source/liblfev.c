/*
 * Copyright 1996-2006 Catherine Loader.
 */

#include "mex.h"
/*
 * Copyright 1996-2006 Catherine Loader.
 */
#include "lfev.h"

extern void fitoptions();

static double hmin, gmin, sig2, pen, vr, tb;
static lfit *blf;
static design *bdes;

int procvbind(des,lf,v)
design *des;
lfit *lf;
int v;
{ double s0, s1, bi;
  int i, ii, k;
  k = procv_var(des,lf,v);
  wdiag(&lf->lfd, &lf->sp, des,des->wd,&lf->dv,0,1,0);
  s0 = s1 = 0.0;
  for (i=0; i<des->n; i++)
  { ii = des->ind[i];
    s0+= prwt(&lf->lfd,ii)*des->wd[i]*des->wd[i];
    bi = prwt(&lf->lfd,ii)*fabs(des->wd[i]*ipower(dist(des,ii),deg(&lf->sp)+1));
    s1+= bi*bi;
  }
  vr += s0;
  tb += s1;
  return(k);
}

double bcri(h,c,cri)
double h;
int c, cri;
{ double num, den, res[10];
  int (*pv)();
  if (c==DALP)
    blf->sp.nn = h;
  else
    blf->sp.fixh = h;
  if ((cri&63)==BIND)
  { pv = procvbind;
    vr = tb = 0.0;
  }
  else pv = procvstd;
  if (cri<64) startlf(bdes,blf,pv,0);
  switch(cri&63)
  { case BGCV:
      ressumm(blf,bdes,res);
      num = -2*blf->lfd.n*res[0];
      den = blf->lfd.n-res[1];
      return(num/(den*den));
    case BCP:
      ressumm(blf,bdes,res);
      return(-2*res[0]/sig2-blf->lfd.n+pen*res[1]);
    case BIND:
      return(vr+pen*pen*tb);
  } 
  LERR(("bcri: unknown criterion"));
  return(0.0);
}

void bsel2(h0,g0,ifact,c,cri)
double h0, g0, ifact;
int c, cri;
{ int done, inc;
  double h1, g1;
  h1 = h0; g1 = g0;
  done = inc = 0;
  while (!done)
  { h1 *= 1+ifact;
    g0 = g1;
    g1 = bcri(h1,c,cri);
    if (g1<gmin) { hmin = h1; gmin = g1; }
    if (g1>g0) inc++; else inc = 0;
    switch(cri)
    { case BIND:
        done = (inc>=4) & (vr<blf->fp.nv);
        break;
      default:
        done = (inc>=4);
    }
  }
}

void bsel3(h0,g0,ifact,c,cri)
double h0, g0, ifact;
int c, cri;
{ double h1, g1;
  int i;
  hmin = h0; gmin = g0;
  for (i=-1; i<=1; i++) if (i!=0)
  { h1 = h0*(1+i*ifact);
    g1 = bcri(h1,c,cri);
    if (g1<gmin) { hmin = h1; gmin = g1; }
  }
  return;
}

void bselect(lf,des,c,cri,pn)
lfit *lf;
design *des;
int c, cri;
double pn;
{ double h0, g0, ifact;
  int i;
  pen = pn;
  blf = lf;
  bdes = des;
  if (cri==BIND) pen /= factorial(deg(&lf->sp)+1);
  hmin = h0 = (c==DFXH) ? fixh(&lf->sp) : nn(&lf->sp);
  if (h0==0) LERR(("bselect: initial bandwidth is 0"));
  if (lf_error) return;
  sig2 = 1.0;

  gmin = g0 = bcri(h0,c,cri);
  if (cri==BCP)
  { sig2 = rv(&lf->fp);
    g0 = gmin = bcri(h0,c,cri+64);
  }
  
  ifact = 0.3;
  bsel2(h0,g0,ifact,c,cri);

  for (i=0; i<5; i++)
  { ifact = ifact/2;
    bsel3(hmin,gmin,ifact,c,cri);
  }
  if (c==DFXH)
    fixh(&lf->sp) = hmin;
  else
    nn(&lf->sp) = hmin;
  startlf(des,lf,procvstd,0);
  ressumm(lf,des,lf->fp.kap);
}

double compsda(x,h,n)
double *x, h;
int n;
/* n/(n-1) * int( fhat''(x)^2 dx ); bandwidth h */
{ int i, j;
  double ik, sd, z;
  ik = wint(1,NULL,0,WGAUS);
  sd = 0;

  for (i=0; i<n; i++)
    for (j=i; j<n; j++)
    { z = (x[i]-x[j])/h;
      sd += (2-(i==j))*Wconv4(z,WGAUS)/(ik*ik);
    }
  sd = sd/(n*(n-1)*h*h*h*h*h);
  return(sd);
}

double widthsj(x,lambda,n)
double *x, lambda;
int n;
{ double ik, a, b, td, sa, z, c, c1, c2, c3;
  int i, j;
  a = GFACT*0.920*lambda*exp(-log((double)n)/7)/SQRT2;
  b = GFACT*0.912*lambda*exp(-log((double)n)/9)/SQRT2;
  ik = wint(1,NULL,0,WGAUS);

  td = 0;
  for (i=0; i<n; i++)
    for (j=i; j<n; j++)
    { z = (x[i]-x[j])/b;
      td += (2-(i==j))*Wconv6(z,WGAUS)/(ik*ik);
    }

  td = -td/(n*(n-1));
  j = 2.0;
  c1 = Wconv4(0.0,WGAUS);
  c2 = wint(1,&j,1,WGAUS);
  c3 = Wconv(0.0,WGAUS);  /* (2*c1/(c2*c3))^(1/7)=1.357 */
  sa = compsda(x,a,n);
  c = b*exp(log(c1*ik/(c2*c3*GFACT*GFACT*GFACT*GFACT)*sa/td)/7)*SQRT2;
  return(c);
}

void kdecri(x,h,res,c,k,ker,n)
double *x, h, *res, c;
int k, ker, n;
{ int i, j;
  double degfree, dfd, pen, s, r0, r1, d0, d1, ik, wij;

  if (h<=0) WARN(("kdecri, h = %6.4f",h));

  res[0] = res[1] = 0.0;
  ik = wint(1,NULL,0,ker);
  switch(k)
  { case 1: /* aic */
      pen = 2.0;
      for (i=0; i<n; i++)
      { r0 = d0 = 0.0;
        for (j=0; j<n; j++)
        { s = (x[i]-x[j])/h;
          r0 += W(s,ker);
          d0 += s*s*Wd(s,ker);
        }
        d0 = -(d0+r0)/(n*h*h*ik);  /* d0 = d/dh fhat(xi) */
        r0 /= n*h*ik;              /* r0 = fhat(xi) */
        res[0] += -2*log(r0)+pen*W(0.0,ker)/(n*h*ik*r0);
        res[1] += -2*d0/r0-pen*W(0.0,ker)/(n*h*ik*r0)*(d0/r0+1.0/h);
      }
      return;
    case 2: /* ocv */
      for (i=0; i<n; i++)
      { r0 = 0.0; d0 = 0.0;
        for (j=0; j<n; j++) if (i!=j)
        { s = (x[i]-x[j])/h;
          r0 += W(s,ker);
          d0 += s*s*Wd(s,ker);
        }
        d0 = -(d0+r0)/((n-1)*h*h);
        r0 = r0/((n-1)*h);
        res[0] -= log(r0);
        res[1] -= d0/r0;
      }
      return;
    case 3: /* lscv */
      r0 = r1 = d0 = d1 = degfree = 0.0;
      for (i=0; i<n; i++)
      { dfd = 0.0;
        for (j=0; j<n; j++)
        { s = (x[i]-x[j])/h;
          wij = W(s,ker);
          dfd += wij;
/* 
 *  r0 = \int fhat * fhat = sum_{i,j} W*W( (Xi-Xj)/h ) / n^2 h.
 *  d0 is it's derivative wrt h.
 *
 *  r1 = 1/n sum( f_{-i}(X_i) ).
 *  d1 is  it's derivative wrt h.
 *
 *  degfree = sum_i ( W_0 / sum_j W( (Xi-Xj)/h ) ) is fitted d.f.
 */
          r0 += Wconv(s,ker);
          d0 += -s*s*Wconv1(s,ker);
          if (i != j)
          { r1 += wij;
            d1 += -s*s*Wd(s,ker);
          }
        }
        degfree += 1.0/dfd;
      }
      d1 -= r1;
      d0 -= r0;
      res[0] = r0/(n*n*h*ik*ik)   - 2*r1/(n*(n-1)*h*ik);
      res[1] = d0/(n*n*h*h*ik*ik) - 2*d1/(n*(n-1)*h*h*ik);
      res[2] = degfree;
      return;
    case 4: /* bcv */
      r0 = d0 = 0.0;
      for (i=0; i<n; i++)
        for (j=i+1; j<n; j++)
        { s = (x[i]-x[j])/h;
          r0 += 2*Wconv4(s,ker);
          d0 += 2*s*Wconv5(s,ker);
        }
      d0 = (-d0-r0)/(n*n*h*h*ik*ik);
      r0 = r0/(n*n*h*ik*ik);
      j = 2.0;
      d1 = wint(1,&j,1,ker);
      r1 = Wconv(0.0,ker);
      res[0] = (d1*d1*r0/4+r1/(n*h))/(ik*ik);
      res[1] = (d1*d1*d0/4-r1/(n*h*h))/(ik*ik);
      return;
    case 5: /* sjpi */
      s = c*exp(5*log(h)/7)/SQRT2;
      d0 = compsda(x,s,n);
      res[0] = d0; /* this is S(alpha) in SJ */
      res[1] = exp(log(Wikk(WGAUS,0)/(d0*n))/5)-h;
      return;
    case 6: /* gas-k-k */
      s = exp(log(1.0*n)/10)*h;
      d0 = compsda(x,s,n);
      res[0] = d0;
      res[1] = exp(log(Wikk(WGAUS,0)/(d0*n))/5)-h;
      return;
  }
  LERR(("kdecri: what???"));
  return;
}

double esolve(x,j,h0,h1,k,c,ker,n)
double *x, h0, h1, c;
int j, k, ker, n;
{ double h[7], d[7], r[7], res[4], min, minh, fact;
  int i, nc;
  min = 1.0e30; minh = 0.0;
  fact = 1.00001;
  h[6] = h0; kdecri(x,h[6],res,c,j,ker,n);
  r[6] = res[0]; d[6] = res[1];
  if (lf_error) return(0.0);
  nc = 0;
  for (i=0; i<k; i++)
  { h[5] = h[6]; r[5] = r[6]; d[5] = d[6];
    h[6] = h0*exp((i+1)*log(h1/h0)/k);
    kdecri(x,h[6],res,c,j,ker,n);
    r[6] = res[0]; d[6] = res[1];
    if (lf_error) return(0.0);
    if (d[5]*d[6]<0)
    { h[2] = h[0] = h[5]; d[2] = d[0] = d[5]; r[2] = r[0] = r[5];
      h[3] = h[1] = h[6]; d[3] = d[1] = d[6]; r[3] = r[1] = r[6];
      while ((h[3]>fact*h[2])|(h[2]>fact*h[3]))
      { h[4] = h[3]-d[3]*(h[3]-h[2])/(d[3]-d[2]);
        if ((h[4]<h[0]) | (h[4]>h[1])) h[4] = (h[0]+h[1])/2;
        kdecri(x,h[4],res,c,j,ker,n);
        r[4] = res[0]; d[4] = res[1];
        if (lf_error) return(0.0);
        h[2] = h[3]; h[3] = h[4];
        d[2] = d[3]; d[3] = d[4];
        r[2] = r[3]; r[3] = r[4];
        if (d[4]*d[0]>0) { h[0] = h[4]; d[0] = d[4]; r[0] = r[4]; }
                    else { h[1] = h[4]; d[1] = d[4]; r[1] = r[4]; }
      }
      if (j>=4) return(h[4]); /* first min for BCV etc */
      if (r[4]<=min) { min = r[4]; minh = h[4]; }
      nc++;
    }
  }
  if (nc==0) minh = (r[5]<r[6]) ? h0 : h1;
  return(minh);
}

void kdeselect(band,x,ind,h0,h1,meth,nm,ker,n)
double h0, h1, *band, *x;
int *ind, nm, ker, n, *meth;
{ double scale, c;
  int i, k;
  k = n/4;
  for (i=0; i<n; i++) ind[i] = i;
  scale = kordstat(x,n+1-k,n,ind) - kordstat(x,k,n,ind);
  c = widthsj(x,scale,n);
  for (i=0; i<nm; i++)
    band[i] = esolve(x,meth[i],h0,h1,10,c,ker,n);
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
/*
 *   The function dens_integrate(lf,des,z) is used to integrate a density
 *   estimate (z=1) or the density squared (z=2). This is used to renormalize
 *   the estimate (function dens_renorm) or in the computation of LSCV
 *   (in modlscv.c). The implementation is presently for d=1.
 *
 *   The computation orders the fit points selected by locfit, and
 *   integrates analytically over each interval. For the log-link,
 *   the interpolant used is peicewise quadratic (with one knot in
 *   the middle of each interval); this differs from the cubic interpolant
 *   used elsewhere in Locfit.
 *
 *   TODO: allow for xlim. What can be done simply in >=2 dimensions?
 */

#include "lfev.h"

/*
 * Finds the order of observations in the array x, and
 * stores in integer array ind.
 * At input, lset l=0 and r=length(x)-1.
 * At output, x[ind[0]] <= x[ind[1]] <= ...
 */
void lforder(ind,x,l,r)
int *ind, l, r;
double *x;
{ double piv;
  int i, i0, i1;
  piv = (x[ind[l]]+x[ind[r]])/2;
  i0 = l; i1 = r;
  while (i0<=i1)
  { while ((i0<=i1) && (x[ind[i0]]<=piv)) i0++;
    while ((i0<=i1) && (x[ind[i1]]>piv))  i1--;
    if (i0<i1)
    { ISWAP(ind[i0],ind[i1]);
      i0++; i1--;
    }
  }
  /* now, x[ind[l..i1]] <= piv < x[ind[i0..r]].
     put the ties in the middle */
  while ((i1>=l) && (x[ind[i1]]==piv)) i1--;
  for (i=l; i<=i1; i++)
    if (x[ind[i]]==piv)
    { ISWAP(ind[i],ind[i1]);
      while (x[ind[i1]]==piv) i1--;
    }

  if (l<i1) lforder(ind,x,l,i1);
  if (i0<r) lforder(ind,x,i0,r);
}

/*
 *  estdiv integrates the density between fit points x0 and x1.
 *  f0, f1 are function values, d0, d1 are derivatives.
 */
double estdiv(x0,x1,f0,f1,d0,d1,lin)
double x0, x1, f0, f1, d0, d1;
int lin;
{ double cf[4], I[2], dlt, e0, e1;

  if (x0==x1) return(0.0);

  if (lin==LIDENT)
  {
/* cf are integrals of hermite polynomials.
 * Then adjust for x1-x0.
 */
    cf[0] = cf[1] = 0.5;
    cf[2] = 1.0/12.0; cf[3] = -cf[2];
    return( (cf[0]*f0+cf[1]*f1)*(x1-x0)
          + (cf[2]*d0+cf[3]*d1)*(x1-x0)*(x1-x0) );
  }

/*
 * this is for LLOG
 */

  dlt = (x1-x0)/2;
  cf[0] = f0;
  cf[1] = d0;
  cf[2] = ( 2*(f1-f0) - dlt*(d1+3*d0) )/(4*dlt*dlt);
  recurint(0.0,dlt,cf,I,0,WRECT);
  e0 = I[0];

  cf[0] = f1;
  cf[1] = -d1;
  cf[2] = ( 2*(f0-f1) + dlt*(d0+3*d1) )/( 4*dlt*dlt );
  recurint(0.0,dlt,cf,I,0,WRECT);
  e1 = I[0];

  return(e0+e1);
}

/*
 *   Evaluate the integral of the density estimate to the power z.
 *   This would be severely messed up, if I ever implement parcomp
 *   for densities.
 */
double dens_integrate(lf,des,z)
lfit *lf;
design *des;
int z;
{ int has_deriv, i, i0, i1, nv, *ind;
  double *xev, *fit, *deriv, sum, term;
  double d0, d1, f0, f1;
  fitpt *fp;

  fp = &lf->fp;

  if (fp->d >= 2)
  { WARN(("dens_integrate requires d=1"));
    return(0.0);
  }

  has_deriv = (deg(&lf->sp) > 0); /* not right? */
  fit = fp->coef;
  if (has_deriv)
    deriv = &fit[fp->nvm];
  xev = evp(fp);

  /*
   * order the vertices
   */
  nv = fp->nv;
  if (lf->lfd.n<nv) return(0.0);
  ind = des->ind;
  for (i=0; i<nv; i++) ind[i] = i;
  lforder(ind,xev,0,nv-1);
  sum = 0.0;

  /*
   * Estimate the contribution of the boundaries.
   * should really check flim here.
   */
  i0 = ind[0]; i1 = ind[1];
  f1 = fit[i0];
  d1 = (has_deriv) ? deriv[i0] :
         (fit[i1]-fit[i0])/(xev[i1]-xev[i0]);
  if (d1 <= 0) WARN(("dens_integrate - ouch!"));
  if (z==2)
  { if (link(&lf->sp)==LLOG)
    { f1 *= 2; d1 *= 2; }
    else
    { d1 = 2*d1*f1; f1 = f1*f1; }
  }
  term = (link(&lf->sp)==LIDENT) ? f1*f1/(2*d1) : exp(f1)/d1;
  sum += term;

  i0 = ind[nv-2]; i1 = ind[nv-1];
  f0 = fit[i1];
  d0 = (has_deriv) ? deriv[i1] :
         (fit[i1]-fit[i0])/(xev[i1]-xev[i0]);
  if (d0 >= 0) WARN(("dens_integrate - ouch!"));
  if (z==2)
  { if (link(&lf->sp)==LLOG)
    { f0 *= 2; d0 *= 2; }
    else
    { d0 = 2*d0*f0; f0 = f0*f0; }
  }
  term = (link(&lf->sp)==LIDENT) ? -f0*f0/(2*d0) : exp(f0)/d0;
  sum += term;
  
  for (i=1; i<nv; i++)
  { i0 = ind[i-1]; i1 = ind[i];
    f0 = fit[i0]; f1 = fit[i1];
    d0 = (has_deriv) ? deriv[i0] :
              (f1-f0)/(xev[i1]-xev[i0]);
    d1 = (has_deriv) ? deriv[i1] : d0;
    if (z==2)
    { if (link(&lf->sp)==LLOG)
      { f0 *= 2; f1 *= 2; d0 *= 2; d1 *= 2; }
      else
      { d0 *= 2*f0; d1 *= 2*f1; f0 = f0*f0; f1 = f1*f1; }
    }
    term = estdiv(xev[i0],xev[i1],f0,f1,d0,d1,link(&lf->sp));
    sum += term;
  }

  return(sum);
}

void dens_renorm(lf,des)
lfit *lf;
design *des;
{ int i;
  double sum;
  sum = dens_integrate(lf,des,1);
  if (sum==0.0) return;
  sum = log(sum);
  for (i=0; i<lf->fp.nv; i++) lf->fp.coef[i] -= sum;
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
/*
 *   This file contains functions for constructing and
 *   interpolating the adaptive tree structure. This is
 *   the default evaluation structure used by Locfit.
 */

#include "lfev.h"

/*
  Guess the number of fitting points.
  Needs improving!
*/
void atree_guessnv(evs,nvm,ncm,vc,d,alp)
evstruc *evs;
double alp;
int *nvm, *ncm, *vc, d;
{ double a0, cu, ifl;
  int i, nv, nc;

  *ncm = 1<<30; *nvm = 1<<30;
  *vc = 1 << d;

  if (alp>0)
  { a0 = (alp > 1) ? 1 : 1/alp;
    if (cut(evs)<0.01)
    { WARN(("guessnv: cut too small."));
      cut(evs) = 0.01;
    }
    cu = 1;
    for (i=0; i<d; i++) cu *= MIN(1.0,cut(evs));
    nv = (int)((5*a0/cu)**vc);   /* this allows 10*a0/cu splits */
    nc = (int)(10*a0/cu+1);      /* and 10*a0/cu cells */
    if (nv<*nvm) *nvm = nv;
    if (nc<*ncm) *ncm = nc;
  }

  if (*nvm == 1<<30) /* by default, allow 100 splits */
  { *nvm = 102**vc;
    *ncm = 201;
  }

  /* inflation based on mk */
  ifl = mk(evs)/100.0;
  *nvm = *vc+(int)(ifl**nvm);
  *ncm = (int)(ifl**ncm);
  
}

/*
  Determine whether a cell in the tree needs splitting.
  If so, return the split variable (0..d-1).
  Otherwise, return -1.
*/
int atree_split(lf,ce,le,ll,ur)
lfit *lf;
int *ce;
double *le, *ll, *ur;
{ int d, vc, i, is;
  double h, hmin, score[MXDIM];
  d = lf->fp.d; vc = 1<<d;

  hmin = 0.0;
  for (i=0; i<vc; i++)
  { h = lf->fp.h[ce[i]];
    if ((h>0) && ((hmin==0)|(h<hmin))) hmin = h;
  }

  is = 0;
  for (i=0; i<d; i++)
  { le[i] = (ur[i]-ll[i])/lf->lfd.sca[i];
    if ((lf->lfd.sty[i]==STCPAR) || (hmin==0))
      score[i] = 2*(ur[i]-ll[i])/(lf->evs.fl[i+d]-lf->evs.fl[i]);
    else
      score[i] = le[i]/hmin;
    if (score[i]>score[is]) is = i;
  }
  if (cut(&lf->evs)<score[is]) return(is);
  return(-1);
}

/*
  recursively grow the tree structure, begining with the parent cell.
*/
void atree_grow(des,lf,ce,ct,term,ll,ur)
design *des;
lfit *lf;
int *ce, *ct, *term;
double *ll, *ur;
{ int nce[1<<MXDIM];
  int i, i0, i1, d, vc, pv, tk, ns;
  double le[MXDIM], z;
  d = lf->fp.d; vc = 1<<d;

  /* does this cell need splitting?
     If not, wrap up and return.
  */
  ns = atree_split(lf,ce,le,ll,ur);
  if (ns==-1)
  { if (ct != NULL) /* reconstructing terminal cells */
    { for (i=0; i<vc; i++) term[*ct*vc+i] = ce[i];
      (*ct)++;
    }
    return;
  }

  /* split the cell at the midpoint on side ns */
  tk = 1<<ns;
  for (i=0; i<vc; i++)
  { if ((i&tk)==0) nce[i] = ce[i];
    else
    { i0 = ce[i];
      i1 = ce[i-tk];
      pv = (lf->lfd.sty[i]!=STCPAR) &&
           (le[ns] < (cut(&lf->evs)*MIN(lf->fp.h[i0],lf->fp.h[i1])));
      nce[i] = newsplit(des,lf,i0,i1,pv);
      if (lf_error) return;
    }
  }
  z = ur[ns]; ur[ns] = (z+ll[ns])/2;
  atree_grow(des,lf,nce,ct,term,ll,ur);
  if (lf_error) return;
  ur[ns] = z;
  for (i=0; i<vc; i++)
    nce[i] = ((i&tk)== 0) ? nce[i+tk] : ce[i];
  z = ll[ns]; ll[ns] = (z+ur[ns])/2;
  atree_grow(des,lf,nce,ct,term,ll,ur);
  if (lf_error) return;
  ll[ns] = z;
}

void atree_start(des,lf)
design *des;
lfit *lf;
{ int d, i, j, k, vc, ncm, nvm;
  double ll[MXDIM], ur[MXDIM];

  if (lf_debug>1) mut_printf(" In atree_start\n");
  d = lf->fp.d;
  atree_guessnv(&lf->evs,&nvm,&ncm,&vc,d,nn(&lf->sp));
  if (lf_debug>2) mut_printf(" atree_start: nvm %d ncm %d\n",nvm,ncm);
  trchck(lf,nvm,ncm,vc);

  /* Set the lower left, upper right limits. */
  for (j=0; j<d; j++)
  { ll[j] = lf->evs.fl[j];
    ur[j] = lf->evs.fl[j+d];
  }

  /* Set the initial cell; fit at the vertices. */
  for (i=0; i<vc; i++)
  { j = i;
    for (k=0; k<d; ++k)
    { evptx(&lf->fp,i,k) = (j%2) ? ur[k] : ll[k];
      j >>= 1;
    }
    lf->evs.ce[i] = i;
    PROC_VERTEX(des,lf,i);
    if (lf_error) return;
    lf->evs.s[i] = 0;
  }
  lf->fp.nv = vc;

  /* build the tree */
  atree_grow(des,lf,lf->evs.ce,NULL,NULL,ll,ur);
  lf->evs.nce = 1;
}

double atree_int(lf,x,what)
lfit *lf;
double *x;
int what;
{ double vv[64][64], *ll, *ur, h, xx[MXDIM];
  int lo, tk, ns, nv, nc, d, i, vc;
  int ce[64];
  fitpt *fp;
  evstruc *evs;

  fp = &lf->fp;
  evs= &lf->evs;
  d = fp->d;
  vc = 1<<d;

  for (i=0; i<vc; i++)
  { setzero(vv[i],vc);
    nc = exvval(fp,vv[i],i,d,what,1);
    ce[i] = evs->ce[i];
  }
  ns = 0;
  while(ns!=-1)
  { ll = evpt(fp,ce[0]); ur = evpt(fp,ce[vc-1]);
    ns = atree_split(lf,ce,xx,ll,ur);
    if (ns!=-1)
    { tk = 1<<ns;
      h = ur[ns]-ll[ns];
      lo = (2*(x[ns]-ll[ns])) < h;
      for (i=0; i<vc; i++) if ((tk&i)==0)
      { nv = findpt(fp,evs,(int)ce[i],(int)ce[i+tk]);
        if (nv==-1) LERR(("Descend tree problem"));
        if (lf_error) return(0.0);
        if (lo)
        { ce[i+tk] = nv;
          if (evs->s[nv]) exvvalpv(vv[i+tk],vv[i],vv[i+tk],d,ns,h,nc);
                    else exvval(fp,vv[i+tk],nv,d,what,1);
        }
        else
        { ce[i] = nv;
          if (evs->s[nv]) exvvalpv(vv[i],vv[i],vv[i+tk],d,ns,h,nc);
                    else exvval(fp,vv[i],nv,d,what,1);
      } }
  } }
  ll = evpt(fp,ce[0]); ur = evpt(fp,ce[vc-1]);
  return(rectcell_interp(x,vv,ll,ur,d,nc));
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
#include "lfev.h"

double linear_interp(h,d,f0,f1)
double h, d, f0, f1;
{ if (d==0) return(f0);
  return( ( (d-h)*f0 + h*f1 ) / d );
}

void hermite2(x,z,phi)
double x, z, *phi;
{ double h;
  if (z==0)
  { phi[0] = 1.0; phi[1] = phi[2] = phi[3] = 0.0;
    return;
  }
  h = x/z;
  if (h<0)
  { phi[0] = 1; phi[1] = 0;
    phi[2] = h; phi[3] = 0;
    return;
  }
  if (h>1)
  { phi[0] = 0; phi[1] = 1;
    phi[2] = 0; phi[3] = h-1;
    return;
  }
  phi[1] = h*h*(3-2*h);
  phi[0] = 1-phi[1];
  phi[2] = h*(1-h)*(1-h);
  phi[3] = h*h*(h - 1);
}

double cubic_interp(h,f0,f1,d0,d1)
double h, f0, f1, d0, d1;
{ double phi[4];
  hermite2(h,1.0,phi);
  return(phi[0]*f0+phi[1]*f1+phi[2]*d0+phi[3]*d1);
}

double cubintd(h,f0,f1,d0,d1)
double h, f0, f1, d0, d1;
{ double phi[4];
  phi[1] = 6*h*(1-h);
  phi[0] = -phi[1];
  phi[2] = (1-h)*(1-3*h);
  phi[3] = h*(3*h-2);
  return(phi[0]*f0+phi[1]*f1+phi[2]*d0+phi[3]*d1);
}

/*
  interpolate over a rectangular cell.
    x = interpolation point. 
    vv = array of vertex values.
    ll = lower left corner.
    ur = upper right corner.
    d = dimension.
    nc = no of coefficients.
*/
double rectcell_interp(x,vv,ll,ur,d,nc)
double *x, vv[64][64], *ll, *ur;
int d, nc;
{ double phi[4];
  int i, j, k, tk;

  tk = 1<<d;
  for (i=0; i<tk; i++) if (vv[i][0]==NOSLN) return(NOSLN);

  /* no derivatives - use multilinear interpolation */
  if (nc==1)
  { for (i=d-1; i>=0; i--)
    { tk = 1<<i;
      for (j=0; j<tk; j++)
        vv[j][0] = linear_interp(x[i]-ll[i],ur[i]-ll[i],vv[j][0],vv[j+tk][0]);
    }
    return(vv[0][0]);
  }

  /* with derivatives -- use cubic */
  if (nc==d+1)
  { for (i=d-1; i>=0; i--)
    { hermite2(x[i]-ll[i],ur[i]-ll[i],phi);
      tk = 1<<i;
      phi[2] *= ur[i]-ll[i];
      phi[3] *= ur[i]-ll[i];
      for (j=0; j<tk; j++)
      { vv[j][0] = phi[0]*vv[j][0] + phi[1]*vv[j+tk][0]
                 + phi[2]*vv[j][i+1] + phi[3]*vv[j+tk][i+1];
        for (k=1; k<=i; k++)
          vv[j][k] = phi[0]*vv[j][k] + phi[1]*vv[j+tk][k];
      }
    }
    return(vv[0][0]); 
  }

  /* with all coefs -- use multicubic */
  for (i=d-1; i>=0; i--)
  { hermite2(x[i]-ll[i],ur[i]-ll[i],phi);
    tk = 1<<i;
    phi[2] *= ur[i]-ll[i];
    phi[3] *= ur[i]-ll[i];
    for (j=0; j<tk; j++)
      for (k=0; k<tk; k++)
        vv[j][k] = phi[0]*vv[j][k] + phi[1]*vv[j+tk][k]
                 + phi[2]*vv[j][k+tk] + phi[3]*vv[j+tk][k+tk];
  }
  return(vv[0][0]);
}

int exvval(fp,vv,nv,d,what,z)
fitpt *fp;
double *vv;
int nv, d, z, what;
{ int i, k;
  double *values;

  k = (z) ? 1<<d : d+1;
  for (i=1; i<k; i++) vv[i] = 0.0;
  switch(what)
  { case PCOEF:
      values = fp->coef;
      break;
    case PVARI:
    case PNLX:
      values = fp->nlx;
      break;
    case PT0:
      values = fp->t0;
      break;
    case PBAND:
      vv[0] = fp->h[nv];
      return(1);
    case PDEGR:
      vv[0] = fp->deg[nv];
      return(1);
    case PLIK:
      vv[0] = fp->lik[nv];
      return(1);
    case PRDF:
      vv[0] = fp->lik[2*fp->nvm+nv];
      return(1);
    default:
      LERR(("Invalid what in exvval"));
      return(0);
  }
  vv[0] = values[nv];
  if (!fp->hasd) return(1);
  if (z)
  { for (i=0; i<d; i++) vv[1<<i] = values[(i+1)*fp->nvm+nv];
    return(1<<d);
  }
  else
  { for (i=1; i<=d; i++) vv[i] = values[i*fp->nvm+nv];
    return(d+1);
  }
}

void exvvalpv(vv,vl,vr,d,k,dl,nc)
double *vv, *vl, *vr, dl;
int d, k, nc;
{ int i, tk, td;
  double f0, f1;
  if (nc==1)
  { vv[0] = (vl[0]+vr[0])/2;
    return;
  }
  tk = 1<<k;
  td = 1<<d;
  for (i=0; i<td; i++) if ((i&tk)==0)
  { f0 = (vl[i]+vr[i])/2 + dl*(vl[i+tk]-vr[i+tk])/8;
    f1 = 1.5*(vr[i]-vl[i])/dl - (vl[i+tk]+vr[i+tk])/4;
    vv[i] = f0;
    vv[i+tk] = f1;
  }
} 

double grid_int(fp,evs,x,what)
fitpt *fp;
evstruc *evs;
double *x;
int what;
{ int d, i, j, jj, nc, sk, v[MXDIM], vc, z0, nce[1<<MXDIM], *mg;
  double *ll, *ur, vv[64][64], z;

  d = fp->d;
  ll = evpt(fp,0); ur = evpt(fp,fp->nv-1);
  mg = mg(evs);

  z0 = 0; vc = 1<<d;
  for (j=d-1; j>=0; j--)
  { v[j] = (int)((mg[j]-1)*(x[j]-ll[j])/(ur[j]-ll[j]));
    if (v[j]<0) v[j]=0;
    if (v[j]>=mg[j]-1) v[j] = mg[j]-2;
    z0 = z0*mg[j]+v[j];
  }
  nce[0] = z0; nce[1] = z0+1; sk = jj = 1; 
  for (i=1; i<d; i++)
  { sk *= mg[i-1];
    jj<<=1;
    for (j=0; j<jj; j++)
      nce[j+jj] = nce[j]+sk;
  }
  for (i=0; i<vc; i++)
    nc = exvval(fp,vv[i],nce[i],d,what,1);
  ll = evpt(fp,nce[0]);
  ur = evpt(fp,nce[vc-1]);
  z = rectcell_interp(x,vv,ll,ur,d,nc);
  return(z);
}

double fitp_int(fp,x,what,i)
fitpt *fp;
double *x;
int what, i;
{ double vv[1+MXDIM];
  exvval(fp,vv,i,fp->d,what,0);
  return(vv[0]);
}

double xbar_int(fp,x,what)
fitpt *fp;
double *x;
int what;
{ int i, nc;
  double vv[1+MXDIM], f;
  nc = exvval(fp,vv,0,fp->d,what,0);
  f = vv[0];
  if (nc>1)
    for (i=0; i<fp->d; i++)
      f += vv[i+1]*(x[i]-evptx(fp,0,i));
  return(f);
}

double dointpoint(lf,x,what,ev,j)
lfit *lf;
double *x;
int what, ev, j;
{ double xf, f, link[LLEN];
  int i, rt;
  fitpt *fp;
  evstruc *evs;

  fp = &lf->fp; evs = &lf->evs;
  for (i=0; i<fp->d; i++) if (lf->lfd.sty[i]==STANGL)
  { xf = floor(x[i]/(2*PI*lf->lfd.sca[i]));
    x[i] -= xf*2*PI*lf->lfd.sca[i];
  }
  if (what > 64)
  { rt = what-64;
    what = PCOEF;
  }
  else rt = 0;

  switch(ev)
  { case EGRID: f = grid_int(fp,evs,x,what); break;
    case EKDTR: f = kdtre_int(fp,evs,x,what); break;
    case ETREE: f = atree_int(lf,x,what); break;
    case EPHULL: f = triang_int(lf,x,what); break;
    case EFITP: f = fitp_int(fp,x,what,j); break;
    case EXBAR: f = xbar_int(fp,x,what); break;
    case ENONE: f = 0; break;
    case ESPHR: f = sphere_int(lf,x,what); break;
    default: LERR(("dointpoint: cannot interpolate structure %d",ev));
  }
  if (((what==PT0)|(what==PNLX)) && (f<0)) f = 0.0;
  f += addparcomp(lf,x,what);

  if (rt>0)
  {
    stdlinks(link,&lf->lfd,&lf->sp,j,f,rsc(fp));
    f = resid(resp(&lf->lfd,j),prwt(&lf->lfd,j),f,fam(&lf->sp),rt,link);
  }

  return(f);
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
/*
 *   Routines for building and interpolating the kd tree.
 *   Initially, this started from the loess code.
 *
 *   Todo: EKDCE isn't working.
 */

#include "lfev.h"

void newcell();
static int nterm;

void kdtre_guessnv(evs,nvm,ncm,vc,n,d,alp)
evstruc *evs;
double alp;
int *nvm, *ncm, *vc, n, d;
{ int k;
  if (ev(evs) == EKDTR)
  { nterm = (int)(cut(evs)/4 * n * MIN(alp,1.0) );
    k = 2*n/nterm;
    *vc = 1<<d;
    *ncm = 2*k+1;
    *nvm = (k+2)**vc/2;
    return;
  }
  if (ev(evs) == EKDCE)
  { nterm = (int)(n * alp);
    *vc = 1;
    *nvm = 1+(int)(2*n/nterm);
    *ncm = 2**nvm+1;
    return;
  }
  *nvm = *ncm = *vc = 0;
  return;
}

/*
  Split x[pi[l..r]] into two equal sized sets.

  Let m=(l+r)/2.
  At return,
    x[pi[l..m]] < x[pi[m+1..r]].
    Assuming no ties:
      If l+r is odd, the sets have the same size.
      If l+r is even, the low set is larger by 1.
    If there are ties, all ties go in the low set.
*/      
int ksmall(l, r, m, x, pi)
int l, r, m, *pi;
double *x;
{
  int il, ir, jl, jr;
  double t;


  while(l<r)
  { t = x[pi[m]];

    /*
      permute the observations so that
        x[pi[l..il]] < t <= x[pi[ir..r]].
    */
    ir = l; il = r;
    while (ir<il)
    { while ((ir<=r) && (x[pi[ir]] < t)) ir++;
      while ((il>=l) && (x[pi[il]]>= t)) il--;
      if (ir<il) ISWAP(pi[ir],pi[il]);
    }

    /*
      move  = t to the middle:
        x[pi[l..il]]  < t
        x[pi[jl..jr]] = t
        x[pi[ir..r]] > t
    */
    jl = ir; jr = r;
    while (ir<jr)
    { while ((ir<=r)  && (x[pi[ir]]== t)) ir++;
      while ((jr>=jl) && (x[pi[jr]] > t)) jr--;
      if (ir<jr) ISWAP(pi[ir],pi[jr]);
    }

    /*
      we're done if m is in the middle, jl <= m <= jr.
    */
    if ((jl<=m) & (jr>=m)) return(jr);

    /*
      update l or r.
    */
    if (m>=ir) l = ir;
    if (m<=il) r = il;
  }
  if (l==r) return(l);
  LERR(("ksmall failure"));
  return(0);
}

int terminal(lf,p,pi,fc,d,m,split_val)
lfit *lf;
int p, d, fc, *m, *pi;
double *split_val;
{ int i, k, lo, hi, split_var;
  double max, min, score, max_score, t;

  /*
    if there are fewer than fc points in the cell, this cell
    is terminal.
  */
  lo = lf->evs.lo[p]; hi = lf->evs.hi[p];
  if (hi-lo < fc) return(-1);

  /* determine the split variable */
  max_score = 0.0; split_var = 0;
  for (k=0; k<d; k++)
  { max = min = datum(&lf->lfd, k, pi[lo]);
    for (i=lo+1; i<=hi; i++)
    { t = datum(&lf->lfd,k,pi[i]);
      if (t<min) min = t;
      if (t>max) max = t;
    }
    score = (max-min) / lf->lfd.sca[k];
    if (score > max_score)
    { max_score = score;
      split_var = k;
    }
  }
  if (max_score==0) /* all points in the cell are equal */
    return(-1);

  *m = ksmall(lo,hi,(lo+hi)/2, dvari(&lf->lfd,split_var), pi);
  *split_val = datum(&lf->lfd, split_var, pi[*m]);

  if (*m==hi) /* all observations go lo */
    return(-1);
  return(split_var);
}

void kdtre_start(des,lf)
design *des;
lfit *lf;
{
  int i, j, vc, d, nc, nv, ncm, nvm, k, m, n, p, *pi;
  double sv;
  d = lf->lfd.d; n = lf->lfd.n; pi = des->ind;
  kdtre_guessnv(&lf->evs,&nvm,&ncm,&vc,n,d,nn(&lf->sp));
  trchck(lf,nvm,ncm,vc);

  nv = 0;
  if (ev(&lf->evs) != EKDCE)
  { for (i=0; i<vc; i++)
    { j = i;
      for (k=0; k<d; ++k)
      { evptx(&lf->fp,i,k) = lf->evs.fl[d*(j%2)+k];
        j >>= 1;
      }
    }
    nv = vc;
    for (j=0; j<vc; j++) lf->evs.ce[j] = j;
  }

  for (i=0; i<n; i++) pi[i] = i;
  p = 0; nc = 1;
  lf->evs.lo[p] = 0; lf->evs.hi[p] = n-1;
  lf->evs.s[p] = -1;
  while (p<nc)
  { k = terminal(lf,p,pi,nterm,d,&m,&sv);
    if (k>=0)
    {
      if ((ncm<nc+2) | (2*nvm<2*nv+vc))
      { WARN(("Insufficient space for full tree"));
        lf->evs.nce = nc; lf->fp.nv = nv;
        return;
      }

      /* new lo cell has obsn's lo[p]..m */
      lf->evs.lo[nc] = lf->evs.lo[p];
      lf->evs.hi[nc] = m;
      lf->evs.s[nc] = -1;

      /* new hi cell has obsn's m+1..hi[p] */
      lf->evs.lo[nc+1] = m+1;
      lf->evs.hi[nc+1] = lf->evs.hi[p];
      lf->evs.s[nc+1] = -1;

      /* cell p is split on variable k, value sv */
      lf->evs.s[p] = k;
      lf->evs.sv[p] = sv;
      lf->evs.lo[p] = nc; lf->evs.hi[p] = nc+1;

      nc=nc+2; i = nv;

      /* now compute the new vertices. */
      if (ev(&lf->evs) != EKDCE)
        newcell(&nv,vc,evp(&lf->fp), d, k, sv,
             &lf->evs.ce[p*vc], &lf->evs.ce[(nc-2)*vc], &lf->evs.ce[(nc-1)*vc]);

    }
    else if (ev(&lf->evs)==EKDCE) /* new vertex at cell center */
    { sv = 0;
      for (i=0; i<d; i++) evptx(&lf->fp,nv,i) = 0;
      for (j=lf->evs.lo[p]; j<=lf->evs.hi[p]; j++)
      { sv += prwt(&lf->lfd,(int)pi[j]);
        for (i=0; i<d; i++)
          evptx(&lf->fp,nv,i) += datum(&lf->lfd,i,pi[j])*prwt(&lf->lfd,(int)pi[j]);
      }
      for (i=0; i<d; i++) evptx(&lf->fp,nv,i) /= sv;
      lf->lfd.n = lf->evs.hi[p] - lf->evs.lo[p] + 1;
      des->ind = &pi[lf->evs.lo[p]]; /* why? */
      PROC_VERTEX(des,lf,nv);
      lf->lfd.n = n; des->ind = pi;
      nv++;
    }
    p++;
  }

  /* We've built the tree. Now do the fitting. */
  if (ev(&lf->evs)==EKDTR)
    for (i=0; i<nv; i++) PROC_VERTEX(des,lf,i);

  lf->evs.nce = nc; lf->fp.nv = nv;
  return;
}

void newcell(nv,vc,xev, d, k, split_val, cpar, clef, crig)
double *xev, split_val;
int *cpar, *clef, *crig;
int *nv, vc, d, k;
{ int i, ii, j, j2, tk, match;
  tk = 1<<k;
  for (i=0; i<vc; i++)
  { if ((i&tk) == 0)
    { for (j=0; j<d; j++) xev[*nv*d+j] = xev[d*cpar[i]+j];
      xev[*nv*d+k] = split_val;
      match = 0; j = vc; /* no matches in first vc points */
      while ((j<*nv) && (!match))
      { j2 = 0;
        while ((j2<d) && (xev[*nv*d+j2] == xev[j*d+j2])) j2++;
        match = (j2==d);
        if (!match) j++;
      }
      ii = i+tk;
      clef[i] = cpar[i];
      clef[ii]= crig[i] = j;
      crig[ii]= cpar[ii];
      if (!match) (*nv)++;
    }
  }
  return;
}

extern void hermite2();

double blend(fp,evs,s,x,ll,ur,j,nt,t,what)
fitpt *fp;
evstruc *evs;
double s, *x, *ll, *ur;
int j, nt, *t, what;
{
  int *ce, k, k1, m, nc, j0, j1;
  double v0, v1, xibar, g0[3], g1[3], gg[4], gp[4], phi[4];
  ce = evs->ce;
  for (k=0; k<4; k++)  /* North South East West */
  { k1 = (k>1);
    v0 = ll[k1]; v1 = ur[k1];
    j0 = ce[j+2*(k==0)+(k==2)];
    j1 = ce[j+3-2*(k==1)-(k==3)];
    xibar = (k%2==0) ? ur[k<2] : ll[k<2];
    m = nt;
    while ((m>=0) && ((evs->s[t[m]] != (k<=1)) | (evs->sv[t[m]] != xibar))) m--;
    if (m >= 0)
    { m = (k%2==1) ? evs->lo[t[m]] : evs->hi[t[m]];
      while (evs->s[m] != -1)
        m = (x[evs->s[m]] < evs->sv[m]) ? evs->lo[m] : evs->hi[m];
      if (v0 < evptx(fp,ce[4*m+2*(k==1)+(k==3)],k1))
      { j0 = ce[4*m+2*(k==1)+(k==3)];
        v0 = evptx(fp,j0,k1);
      }
      if (evptx(fp,ce[4*m+3-2*(k==0)-(k==2)],k1) < v1)
      { j1 = ce[4*m+3-2*(k==0)-(k==2)];
        v1 = evptx(fp,j1,k1);
      }
    }
    nc = exvval(fp,g0,j0,2,what,0);
    nc = exvval(fp,g1,j1,2,what,0);
    if (nc==1)
      gg[k] = linear_interp((x[(k>1)]-v0),v1-v0,g0[0],g1[0]);
    else
    { hermite2(x[(k>1)]-v0,v1-v0,phi);
      gg[k] = phi[0]*g0[0]+phi[1]*g1[0]+(phi[2]*g0[1+k1]+phi[3]*g1[1+k1])*(v1-v0);
      gp[k] = phi[0]*g0[2-k1] + phi[1]*g1[2-k1];
    }
  }
  s = -s;
  if (nc==1)
    for (k=0; k<2; k++)
      s += linear_interp(x[k]-ll[k],ur[k]-ll[k],gg[3-2*k],gg[2-2*k]);
    else
    for (k=0; k<2; k++) /* EW NS */
    { hermite2(x[k]-ll[k],ur[k]-ll[k],phi);
      s += phi[0]*gg[3-2*k] + phi[1]*gg[2-2*k]
          +(phi[2]*gp[3-2*k] + phi[3]*gp[2-2*k]) * (ur[k]-ll[k]);
    }
  return(s);
}

double kdtre_int(fp,evs,x,what)
fitpt *fp;
evstruc *evs;
double *x;
int what;
{
  int *ce, k, vc, t[20], nt, nc, j, d;
  double *ll, *ur, ff, vv[64][64];
  d = fp->d;
  vc = 1<<d;
  if (d > 6) { LERR(("d too large in kdint")); return(NOSLN); }

  /* descend the tree to find the terminal cell */
  nt = 0; t[nt] = 0; k = 0;
  while (evs->s[k] != -1)
  { nt++;
    if (nt>=20) { LERR(("Too many levels in kdint")); return(NOSLN); }
    k = t[nt] = (x[evs->s[k]] < evs->sv[k]) ? evs->lo[k] : evs->hi[k];
  }

  ce = &evs->ce[k*vc];
  ll = evpt(fp,ce[0]);
  ur = evpt(fp,ce[vc-1]);
  nc = 0;
  for (j=0; j<vc; j++) nc = exvval(fp,vv[j],(int)ce[j],d,what,0);
  ff = rectcell_interp(x,vv,ll,ur,d,nc);

  if (d==2) ff = blend(fp,evs,ff,x,ll,ur,k*vc,nt,t,what);
  return(ff);
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
#include "lfev.h"

/*
 * convert eval. structure string to numeric code.
 */
#define NETYPE 11
static char *etype[NETYPE]= { "tree",     "phull", "data", "grid", "kdtree",
                          "kdcenter", "cross", "preset", "xbar", "none",
                          "sphere" };
static int   evals[NETYPE]= { ETREE, EPHULL, EDATA, EGRID, EKDTR,
                          EKDCE, ECROS,  EPRES, EXBAR, ENONE, ESPHR };
int lfevstr(char *z)
{ return(pmatch(z, etype, evals, NETYPE, -1));
}

void evstruc_init(evs)
evstruc *evs;
{ int i;
  ev(evs) = ETREE;
  mk(evs) = 100;
  cut(evs) = 0.8;
  for (i=0; i<MXDIM; i++)
  { evs->fl[i] = evs->fl[i+MXDIM] = 0.0;
    evs->mg[i] = 10;
  }
  evs->nce = evs->ncm = 0;
}

int evstruc_reqi(nvm,ncm,vc)
int nvm, ncm, vc;
{ return(ncm*vc+3*MAX(ncm,nvm));
}

/* al=1: allows dynamic allocation.
 * al=0: inhibit. use with caution.
 */
void evstruc_alloc(evs,nvm,ncm,vc,al)
evstruc *evs;
int nvm, ncm, vc, al;
{ int rw, *k;

  if (al)
  { rw = evstruc_reqi(nvm,ncm,vc);
    if (evs->liw<rw)
    { evs->iwk = (int *)calloc(rw,sizeof(int));
    if ( evs->iwk == NULL ) {
      printf("Problem allocating memory for evs->iwk\n");fflush(stdout);
    }
      evs->liw = rw;
    }
  }
  k = evs->iwk;
  evs->ce = k; k += vc*ncm;
  evs->s  = k; k += MAX(ncm,nvm);
  evs->lo = k; k += MAX(ncm,nvm);
  evs->hi = k; k += MAX(ncm,nvm);
}

void guessnv(evs,sp,mdl,n,d,lw,nvc)
evstruc *evs;
smpar *sp;
module *mdl;
int n, d, *lw, *nvc;
{ int i, nvm, ncm, vc;

  npar(sp) = calcp(sp,d);
  switch(ev(evs))
  { case ETREE:
      atree_guessnv(evs,&nvm,&ncm,&vc,d,nn(sp));
      break;
    case EPHULL:
      nvm = ncm = mk(evs)*d;
      vc = d+1;
      break;
    case EDATA:
    case ECROS:
      nvm = n;
      ncm = vc = 0;
      break;
    case EKDTR:
    case EKDCE:
      kdtre_guessnv(evs,&nvm,&ncm,&vc,n,d,nn(sp));
      break;
    case EGRID:
      nvm = 1;
      for (i=0; i<d; i++) nvm *= evs->mg[i];
      ncm = 0;
      vc = 1<<d;
      break;
    case EXBAR:
    case ENONE:
      nvm = 1;
      ncm = vc = 0;
      break;
    case EPRES:
      nvm = evs->mg[0];
      ncm = vc = 0;
      break;
    default:
      LERR(("guessnv: I don't know this evaluation structure."));
      nvm = ncm = vc = 0;
  }

  lw[0] = des_reqd(n,npar(sp));
  lw[1] = lfit_reqd(d,nvm,ncm,mdl);
  lw[2] = evstruc_reqi(nvm,ncm,vc);
  lw[6] = des_reqi(n,npar(sp));
  lw[3] = pc_reqd(d,npar(sp));
  lw[4] = mdl->keepv;
  lw[5] = mdl->keepc;

  if (nvc==NULL) return;
  nvc[0] = nvm;
  nvc[1] = ncm;
  nvc[2] = vc;
  nvc[3] = nvc[4] = 0;
}

/*
 * trchck checks the working space on the lfit structure 
 * has space for nvm vertices and ncm cells.
 */
void lfit_alloc(lf)
lfit *lf;
{ lf->fp.lwk = lf->fp.lev = lf->evs.liw = lf->pc.lwk = 0;
  lf->lf_init_id = LF_INIT_ID;
}
int lfit_reqd(d,nvm,ncm,mdl)
int d, nvm, ncm;
module *mdl;
{ int z;
  z = mdl->keepv;
  return(nvm*z+ncm);
}

void trchck(lf,nvm,ncm,vc)
lfit *lf;
int nvm, ncm, vc;
{ int rw, d, *k;
  double *z;

  if (lf->lf_init_id != LF_INIT_ID) lfit_alloc(lf);

  lf->fp.nvm = nvm; lf->evs.ncm = ncm;
  d = lf->lfd.d;

  if (lf->fp.lev < d*nvm)
  { lf->fp.xev = (double *)calloc(d*nvm,sizeof(double));
    if ( lf->fp.xev == NULL ) {
      printf("Problem allocating memory for lf->fp.xev\n");fflush(stdout);
    }
    lf->fp.lev = d*nvm;
  }

  rw = lfit_reqd(d,nvm,ncm,&lf->mdl);
  if (lf->fp.lwk < rw)
  {
    lf->fp.coef = (double *)calloc(rw,sizeof(double));
    if ( lf->fp.coef == NULL ) {
      printf("Problem allocating memory for lf->fp.coef\n");fflush(stdout);
    }
    lf->fp.lwk = rw;
  }
  z = lf->fp.wk = lf->fp.coef;

  lf->fp.h = NULL;
  if (!lf->mdl.isset) mut_printf("module not set.\n");
  else
  { if (lf->mdl.alloc!=NULL) lf->mdl.alloc(lf);
    z += KEEPV(lf) * nvm;
  }
  lf->evs.sv = z; z += ncm;

  evstruc_alloc(&lf->evs,nvm,ncm,vc,1);
}

void data_guessnv(nvm,ncm,vc,n)
int *nvm, *ncm, *vc, n;
{ *nvm = n;
  *ncm = *vc = 0;
}

void dataf(des,lf)
design *des;
lfit *lf;
{
  int d, i, j, ncm, nv, vc;

  d = lf->lfd.d;
  data_guessnv(&nv,&ncm,&vc,lf->lfd.n);
  trchck(lf,nv,ncm,vc);

  for (i=0; i<nv; i++)
    for (j=0; j<d; j++) evptx(&lf->fp,i,j) = datum(&lf->lfd,j,i);
  for (i=0; i<nv; i++)
  {
    PROC_VERTEX(des,lf,i);
    lf->evs.s[i] = 0;
  }
  lf->fp.nv = lf->fp.nvm = nv; lf->evs.nce = 0;
}

void xbar_guessnv(nvm,ncm,vc)
int *nvm, *ncm, *vc;
{ *nvm = 1;
  *ncm = *vc = 0;
  return;
}

void xbarf(des,lf)
design *des;
lfit *lf;
{ int i, d, nvm, ncm, vc;
  d = lf->lfd.d;
  xbar_guessnv(&nvm,&ncm,&vc);
  trchck(lf,1,0,0);
  for (i=0; i<d; i++) evptx(&lf->fp,0,i) = lf->pc.xbar[i];
  PROC_VERTEX(des,lf,0);
  lf->evs.s[0] = 0;
  lf->fp.nv = 1; lf->evs.nce = 0;
}

void preset(des,lf)
design *des;
lfit *lf;
{ int i, nv;

  nv = lf->fp.nvm;
  trchck(lf,nv,0,0);
  for (i=0; i<nv; i++)
  { 
    PROC_VERTEX(des,lf,i);
    lf->evs.s[i] = 0;
  }
  lf->fp.nv = nv; lf->evs.nce = 0;
}

void crossf(des,lf)
design *des;
lfit *lf;
{ int d, i, j, n, nv, ncm, vc;
  double w;

  n = lf->lfd.n; d = lf->lfd.d;
  data_guessnv(&nv,&ncm,&vc,n);
  trchck(lf,nv,ncm,vc);

  if (lf->lfd.w==NULL) LERR(("crossf() needs prior weights"));
  for (i=0; i<n; i++)
    for (j=0; j<d; j++) evptx(&lf->fp,i,j) = datum(&lf->lfd,j,i);
  for (i=0; i<n; i++)
  { lf->evs.s[i] = 0;
    w = prwt(&lf->lfd,i);
    lf->lfd.w[i] = 0;
    PROC_VERTEX(des,lf,i);
    lf->lfd.w[i] = w;
  }
  lf->fp.nv = n; lf->evs.nce = 0;
}

void gridf(des,lf)
design *des;
lfit *lf;
{ int d, i, j, nv, u0, u1, z;
  nv = 1; d = lf->lfd.d;
  for (i=0; i<d; i++)
  { if (lf->evs.mg[i]==0)
      lf->evs.mg[i] = 2+(int)((lf->evs.fl[i+d]-lf->evs.fl[i])/(lf->lfd.sca[i]*cut(&lf->evs)));
    nv *= lf->evs.mg[i];
  }
  trchck(lf,nv,0,1<<d);
  for (i=0; i<nv; i++)
  { z = i;
    for (j=0; j<d; j++)
    { u0 = z%lf->evs.mg[j];
      u1 = lf->evs.mg[j]-1-u0;
      evptx(&lf->fp,i,j) = (lf->evs.mg[j]==1) ? lf->evs.fl[j] :
                      (u1*lf->evs.fl[j]+u0*lf->evs.fl[j+d])/(lf->evs.mg[j]-1);
      z = z/lf->evs.mg[j];
    }
    lf->evs.s[i] = 0;
    PROC_VERTEX(des,lf,i);
  }
  lf->fp.nv = nv; lf->evs.nce = 0;
}

int findpt(fp,evs,i0,i1)
fitpt *fp;
evstruc *evs;
int i0, i1;
{ int i;
  if (i0>i1) ISWAP(i0,i1);
  for (i=i1+1; i<fp->nv; i++)
    if ((evs->lo[i]==i0) && (evs->hi[i]==i1)) return(i);
  return(-1);
}

/*
  add a new vertex at the midpoint of (x[i0],x[i1]).
  return the vertex number.
*/
int newsplit(des,lf,i0,i1,pv)
design *des;
lfit *lf;
int i0, i1, pv;
{ int i, nv;

  i = findpt(&lf->fp,&lf->evs,i0,i1);
  if (i>=0) return(i);

  if (i0>i1) ISWAP(i0,i1);
  nv = lf->fp.nv;
  
  /* the point is new. Now check we have space for the new point. */
  if (nv>=lf->fp.nvm)
  {
    LERR(("newsplit: out of vertex space"));
    return(-1);
  }

  /* compute the new point, and evaluate the fit */
  lf->evs.lo[nv] = i0;
  lf->evs.hi[nv] = i1;
  for (i=0; i<lf->fp.d; i++)
    evptx(&lf->fp,nv,i) = (evptx(&lf->fp,i0,i)+evptx(&lf->fp,i1,i))/2;
  if (pv) /* pseudo vertex */
  { lf->fp.h[nv] = (lf->fp.h[i0]+lf->fp.h[i1])/2;
    lf->evs.s[nv] = 1; /* pseudo-vertex */
  }
  else /* real vertex */
  {
    PROC_VERTEX(des,lf,nv);
    lf->evs.s[nv] = 0;
  }
  lf->fp.nv++;

  return(nv);
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
/*
 *   Functions for constructing the fit and
 *   interpolating on the circle/sphere. d=2 only.
 */

#include "lfev.h"

/*
 * Guess the number of fitting points.
 */
void sphere_guessnv(nvm,ncm,vc,mg)
int *nvm, *ncm, *vc, *mg;
{ *nvm = mg[1]*(mg[0]+1);
  *ncm = 0;
  *vc = 0;
}

void sphere_start(des,lf)
design *des;
lfit *lf;
{ int d, i, j, ct, nv, ncm, vc, *mg;
  double rmin, rmax, *orig, r, th, c, s;

  mg = mg(&lf->evs);
  sphere_guessnv(&nv,&ncm,&vc,mg);
  trchck(lf,nv,0,0);
  d = lf->lfd.d;

  rmin = lf->evs.fl[0];
  rmax = lf->evs.fl[1];
  orig = &lf->evs.fl[2];
rmin = 0; rmax = 1; orig[0] = orig[1] = 0.0;

  ct = 0;
  for (i=0; i<mg[1]; i++)
  { th = 2*PI*i/mg[1];
    c = cos(th);
    s = sin(th);
    for (j=0; j<=mg[0]; j++)
    { r = rmin + (rmax-rmin)*j/mg[0];
      evptx(&lf->fp,ct,0) = orig[0] + r*c;
      evptx(&lf->fp,ct,1) = orig[1] + r*s;
      PROC_VERTEX(des,lf,ct);
      ct++;
    }
  }
  lf->fp.nv = ct;
  lf->evs.nce = 0;
}

double sphere_int(lf,x,what)
lfit *lf;
double *x;
int what;
{ double rmin, rmax, *orig, dx, dy, r, th, th0, th1;
  double v[64][64], c0, c1, s0, s1, r0, r1, d0, d1;
  double ll[2], ur[2], xx[2];
  int i0, j0, i1, j1, *mg, nc, ce[4];

  rmin = lf->evs.fl[0];
  rmax = lf->evs.fl[1];
  orig = &lf->evs.fl[2];
rmin = 0; rmax = 1; orig[0] = orig[1] = 0.0;
  mg = mg(&lf->evs);

  dx = x[0] - orig[0];
  dy = x[1] - orig[1];
  r = sqrt(dx*dx+dy*dy);
  th = atan2(dy,dx); /* between -pi and pi */

  i0 = (int)floor(mg[1]*th/(2*PI)) % mg[1];
  j0 = (int)(mg[0]*(r-rmin)/(rmax-rmin));

  i1 = (i0+1) % mg[1];
  j1 = j0+1; if (j1>mg[0]) { j0 = mg[0]-1; j1 = mg[0]; }

  ce[0] = i0*(mg[0]+1)+j0;
  ce[1] = i0*(mg[0]+1)+j1;
  ce[2] = i1*(mg[0]+1)+j0;
  ce[3] = i1*(mg[0]+1)+j1;
  nc = exvval(&lf->fp,v[0],ce[0],2,what,1);
  nc = exvval(&lf->fp,v[1],ce[1],2,what,1);
  nc = exvval(&lf->fp,v[2],ce[2],2,what,1);
  nc = exvval(&lf->fp,v[3],ce[3],2,what,1);

  th0 = 2*PI*i0/mg[1]; c0 = cos(th0); s0 = sin(th0);
  th1 = 2*PI*i1/mg[1]; c1 = cos(th1); s1 = sin(th1);
  r0 = rmin + j0*(rmax-rmin)/mg[0];
  r1 = rmin + j1*(rmax-rmin)/mg[0];
  
  d0 = c0*v[0][1] + s0*v[0][2];
  d1 = r0*(c0*v[0][2]-s0*v[0][1]);
  v[0][1] = d0; v[0][2] = d1;

  d0 = c0*v[1][1] + s0*v[1][2];
  d1 = r1*(c0*v[1][2]-s0*v[1][1]);
  v[1][1] = d0; v[1][2] = d1;

  d0 = c1*v[2][1] + s1*v[2][2];
  d1 = r0*(c1*v[2][2]-s1*v[2][1]);
  v[2][1] = d0; v[2][2] = d1;

  d0 = c1*v[3][1] + s1*v[3][2];
  d1 = r1*(c1*v[3][2]-s1*v[3][1]);
  v[3][1] = d0; v[3][2] = d1;

  xx[0] = r; xx[1] = th;
  ll[0] = r0; ll[1] = th0;
  ur[0] = r1; ur[1] = th1;
  return(rectcell_interp(xx,v,ll,ur,2,nc));
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
#include "lfev.h"

void solve(A,b,d) /* this is crude! A organized by column. */
double *A, *b;
int d;
{ int i, j, k;
  double piv;
  for (i=0; i<d; i++)
  { piv = A[(d+1)*i];
    for (j=i; j<d; j++) A[j*d+i] /= piv;
    b[i] /= piv;
    for (j=0; j<d; j++) if (j != i)
    { piv = A[i*d+j];
      A[i*d+j] = 0;
      for (k=i+1; k<d; k++)
        A[k*d+j] -= piv*A[k*d+i];
      b[j] -= piv*b[i];
    }
  }
}

void triang_guessnv(nvm,ncm,vc,d,mk)
int *nvm, *ncm, *vc, d, mk;
{ *nvm = *ncm = mk*d;
  *vc = d+1;
  return;         
}

int triang_split(lf,ce,le)
lfit *lf;
double *le;
int *ce;
{ int d, i, j, k, nts, vc;
  double di, dfx[MXDIM];
  nts = 0; d = lf->fp.d; vc = d+1;
  for (i=0; i<d; i++)
    for (j=i+1; j<=d; j++)
    { for (k=0; k<d; k++)
        dfx[k] = evptx(&lf->fp,ce[i],k)-evptx(&lf->fp,ce[j],k);
      di = rho(dfx,lf->lfd.sca,d,KSPH,NULL);
      le[i*vc+j] = le[j*vc+i] = di/MIN(lf->fp.h[ce[i]],lf->fp.h[ce[j]]);
      nts = nts || le[i*vc+j]>cut(&lf->evs);
    }
  return(nts);
}

void resort(pv,xev,dig)
double *xev;
int *pv, *dig;
{ double d0, d1, d2;
  int i;
  d0 = d1 = d2 = 0;
  for (i=0; i<3; i++)
  { d0 += (xev[3*pv[11]+i]-xev[3*pv[1]+i])*(xev[3*pv[11]+i]-xev[3*pv[1]+i]);
    d1 += (xev[3*pv[ 7]+i]-xev[3*pv[2]+i])*(xev[3*pv[ 7]+i]-xev[3*pv[2]+i]);
    d2 += (xev[3*pv[ 6]+i]-xev[3*pv[3]+i])*(xev[3*pv[ 6]+i]-xev[3*pv[3]+i]);
  }
  if ((d0<=d1) & (d0<=d2))
  { dig[0] = pv[1]; dig[1] = pv[11];
    dig[2] = pv[2]; dig[3] = pv[7];
    dig[4] = pv[3]; dig[5] = pv[6];
  }
  else if (d1<=d2)
  { dig[0] = pv[2]; dig[1] = pv[7];
    dig[2] = pv[1]; dig[3] = pv[11];
    dig[4] = pv[3]; dig[5] = pv[6];
  }
  else
  { dig[0] = pv[3]; dig[1] = pv[6];
    dig[2] = pv[2]; dig[3] = pv[7];
    dig[4] = pv[1]; dig[5] = pv[11];
  }
}

void triang_grow(des,lf,ce,ct,term)
design *des;
lfit *lf;
int *ce, *ct, *term;
{ double le[(1+MXDIM)*(1+MXDIM)], ml;
  int d, i, j, im, jm, vc, pv[(1+MXDIM)*(1+MXDIM)], dig[6];
  int nce[1+MXDIM];
  if (lf_error) return;
  d = lf->fp.d; vc = d+1;
  if (!triang_split(lf,ce,le))
  { if (ct != NULL)
    { for (i=0; i<vc; i++) term[*ct*vc+i] = ce[i];
      (*ct)++;
    }
    return;
  }
  if (d>3)
  { ml = 0;
    for (i=0; i<d; i++)
      for (j=i+1; j<vc; j++)
        if (le[i*vc+j]>ml) { ml = le[i*vc+j]; im = i; jm = j; }
    pv[0] = newsplit(des,lf,(int)ce[im],(int)ce[jm],0);
    for (i=0; i<vc; i++) nce[i] = ce[i];
    nce[im] = pv[0]; triang_grow(des,lf,nce,ct,term); nce[im] = ce[im];
    nce[jm] = pv[0]; triang_grow(des,lf,nce,ct,term);
    return;
  }

  for (i=0; i<d; i++)
    for (j=i+1; j<=d; j++)
      pv[i*vc+j] = pv[j*vc+i]
        = newsplit(des,lf,(int)ce[i],(int)ce[j],le[i*vc+j]<=cut(&lf->evs));
  for (i=0; i<=d; i++) /* corners */
  { for (j=0; j<=d; j++) nce[j] = (j==i) ? ce[i] : pv[i*vc+j];
    triang_grow(des,lf,nce,ct,term);
  }
  
  if (d==2) /* center for d=2 */
  { nce[0] = pv[5]; nce[1] = pv[2]; nce[2] = pv[1];
    triang_grow(des,lf,nce,ct,term);
  }
  if (d==3) /* center for d=3 */
  { resort(pv,evp(&lf->fp),dig);
    nce[0] = dig[0]; nce[1] = dig[1];
    nce[2] = dig[2]; nce[3] = dig[4]; triang_grow(des,lf,nce,ct,term);
    nce[2] = dig[5]; nce[3] = dig[3]; triang_grow(des,lf,nce,ct,term);
    nce[2] = dig[2]; nce[3] = dig[5]; triang_grow(des,lf,nce,ct,term);
    nce[2] = dig[4]; nce[3] = dig[3]; triang_grow(des,lf,nce,ct,term);
  }
}

void triang_descend(tr,xa,ce)
lfit *tr;
double *xa;
int *ce;
{ double le[(1+MXDIM)*(1+MXDIM)], ml;
  int d, vc, i, j, im, jm, pv[(1+MXDIM)*(1+MXDIM)];
  design *des;
  des = NULL;
  if (!triang_split(tr,ce,le)) return;
  d = tr->fp.d; vc = d+1;

  if (d>3) /* split longest edge */
  { ml = 0;
    for (i=0; i<d; i++)
      for (j=i+1; j<vc; j++)
        if (le[i*vc+j]>ml) { ml = le[i*vc+j]; im = i; jm = j; }
    pv[0] = newsplit(des,tr,(int)ce[im],(int)ce[jm],0);
    if (xa[im]>xa[jm])
    { xa[im] -= xa[jm]; xa[jm] *= 2; ce[jm] = pv[0]; }
    else
    { xa[jm] -= xa[im]; xa[im] *= 2; ce[im] = pv[0]; }
    triang_descend(tr,xa,ce);
    return;
  }

  for (i=0; i<d; i++)
    for (j=i+1; j<=d; j++)
      pv[i*vc+j] = pv[j*vc+i]
        = newsplit(des,tr,(int)ce[i],(int)ce[j],le[i*d+j]<=cut(&tr->evs));
  for (i=0; i<=d; i++) if (xa[i]>=0.5) /* in corner */
  { for (j=0; j<=d; j++)
    { if (i!=j) ce[j] = pv[i*vc+j];
      xa[j] = 2*xa[j];
    }
    xa[i] -= 1;
    triang_descend(tr,xa,ce);
    return;
  }
  if (d==1) { LERR(("weights sum to < 1")); }
  if (d==2) /* center */
  { ce[0] = pv[5]; xa[0] = 1-2*xa[0];
    ce[1] = pv[2]; xa[1] = 1-2*xa[1];
    ce[2] = pv[1]; xa[2] = 1-2*xa[2];
    triang_descend(tr,xa,ce);
  }
  if (d==3) /* center */
  { double z; int dig[6];
    resort(pv,evp(&tr->fp),dig);
    ce[0] = dig[0]; ce[1] = dig[1];
    xa[0] *= 2; xa[1] *= 2; xa[2] *= 2; xa[3] *= 2;
    if (xa[0]+xa[2]>=1)
    { if (xa[0]+xa[3]>=1)
      { ce[2] = dig[2]; ce[3] = dig[4];
        z = xa[0];
        xa[3] += z-1; xa[2] += z-1; xa[0] = xa[1]; xa[1] = 1-z;
      }
      else
      { ce[2] = dig[2]; ce[3] = dig[5];
        z = xa[3]; xa[3] = xa[1]+xa[2]-1; xa[1] = z;
        z = xa[2]; xa[2] += xa[0]-1; xa[0] = 1-z;
    } }
    else
    { if (xa[1]+xa[2]>=1)
      { ce[2] = dig[5]; ce[3] = dig[3];
        xa[1] = 1-xa[1]; xa[2] -= xa[1]; xa[3] -= xa[1];
      }
      else
      { ce[2] = dig[4]; ce[3] = dig[3];
        z = xa[3]; xa[3] += xa[1]-1; xa[1] = xa[2];
        xa[2] = z+xa[0]-1; xa[0] = 1-z;
    } }
    triang_descend(tr,xa,ce);
} }

void covrofdata(lfd,V,mn) /* covar of data; mean in mn */
lfdata *lfd;
double *V, *mn;
{ int d, i, j, k;
  double s;
  s = 0; d = lfd->d;
  for (i=0; i<d*d; i++) V[i] = 0;
  for (i=0; i<lfd->n; i++)
  { s += prwt(lfd,i);
    for (j=0; j<d; j++)
      for (k=0; k<d; k++)
        V[j*d+k] += prwt(lfd,i)*(datum(lfd,j,i)-mn[j])*(datum(lfd,k,i)-mn[k]);
  }
  for (i=0; i<d*d; i++) V[i] /= s;
}

int intri(x,w,xev,xa,d) /* is x in triangle bounded by xd[0..d-1]? */
double *x, *xev, *xa;
int *w, d;
{ int i, j;
  double eps, *r, xd[MXDIM*MXDIM];
  eps = 1.0e-10;
  r = &xev[w[d]*d];
  for (i=0; i<d; i++)
  { xa[i] = x[i]-r[i];
    for (j=0; j<d; j++) xd[i*d+j] = xev[w[i]*d+j]-r[j];
  }
  solve(xd,xa,d);
  xa[d] = 1.0;
  for (i=0; i<d; i++) xa[d] -= xa[i];
  for (i=0; i<=d; i++) if ((xa[i]<-eps) | (xa[i]>1+eps)) return(0);
  return(1);
}

void triang_start(des,lf) /* Triangulation with polyhedral start */
design *des;
lfit *lf;
{
  int i, j, k, n, d, nc, nvm, ncm, vc;
  int *ce, ed[1+MXDIM];
  double V[MXDIM*MXDIM], P[MXDIM*MXDIM], sigma, z[MXDIM], xa[1+MXDIM], *xev;
  xev = evp(&lf->fp);
  d = lf->lfd.d; n = lf->lfd.n;
  lf->fp.nv = nc = 0;

  triang_guessnv(&nvm,&ncm,&vc,d,mk(&lf->evs));
  trchck(lf,nvm,ncm,vc);

  ce = lf->evs.ce;
  for (j=0; j<d; j++) xev[j] = lf->pc.xbar[j];
  lf->fp.nv = 1;
  covrofdata(&lf->lfd,V,xev); /* fix this with scaling */
  eig_dec(V,P,d);

  for (i=0; i<d; i++) /* add vertices +- 2sigma*eigenvect */
  { sigma = sqrt(V[i*(d+1)]);
    for (j=0; j<d; j++)
      xev[lf->fp.nv*d+j] = xev[j]-2*sigma*P[j*d+i];
    lf->fp.nv++;
    for (j=0; j<d; j++)
      xev[lf->fp.nv*d+j] = xev[j]+2*sigma*P[j*d+i];
    lf->fp.nv++;
  }

  for (i=0; i<n; i++) /* is point i inside? */
  { ed[0] = 0;
    for (j=0; j<d; j++)
    { z[j] = 0;
      for (k=0; k<d; k++) z[j] += P[k*d+j]*(datum(&lf->lfd,k,i)-xev[k]);
      ed[j+1] = 2*j+1+(z[j]>0);
      for (k=0; k<d; k++) z[j] = datum(&lf->lfd,j,i);
    }
    k = intri(z,ed,xev,xa,d);
    if (xa[0]<0)
    { for (j=1; j<=d; j++)
        for (k=0; k<d; k++)
          xev[ed[j]*d+k] = xa[0]*xev[k]+(1-xa[0])*xev[ed[j]*d+k];
    }
  }

  nc = 1<<d; /* create initial cells */
  for (i=0; i<nc; i++)
  { ce[i*vc] = 0; k = i;
    for (j=0; j<d; j++)
    { ce[i*vc+j+1] = 2*j+(k%2)+1;
      k>>=1;
    }
  }

  for (i=0; i<lf->fp.nv; i++)
  { PROC_VERTEX(des,lf,i);
    if (lf_error) return;
    lf->evs.s[i] = 0;
  }
  for (i=0; i<nc; i++)
    triang_grow(des,lf,&ce[i*vc],NULL,NULL);
  lf->evs.nce = nc;
}

double triang_cubicint(v,vv,w,d,nc,xxa)
double *v, *vv, *xxa;
int *w, d, nc;
{ double sa, lb, *vert0, *vert1, *vals0, *vals1, deriv0, deriv1;
  int i, j, k;
  if (nc==1) /* linear interpolate */
  { sa = 0;
    for (i=0; i<=d; i++) sa += xxa[i]*vv[i];
    return(sa);
  }
  sa = 1.0;
  for (j=d; j>0; j--)  /* eliminate v[w[j]] */
  { lb = xxa[j]/sa;
    for (k=0; k<j; k++) /* Interpolate edge v[w[k]],v[w[j]] */
    { vert0 = &v[w[k]*d];
      vert1 = &v[w[j]*d];
      vals0 = &vv[k*nc];
      vals1 = &vv[j*nc];
      deriv0 = deriv1 = 0;
      for (i=0; i<d; i++)
      { deriv0 += (vert1[i]-vert0[i])*vals0[i+1];
        deriv1 += (vert1[i]-vert0[i])*vals1[i+1];
      }
      vals0[0] = cubic_interp(lb,vals0[0],vals1[0],deriv0,deriv1);
      for (i=1; i<=d; i++)
        vals0[i] = (1-lb)*((1-lb)*vals0[i]+lb*vals1[i]);
    }
    sa -= xxa[j];
    if (sa<=0) j = 0;
  }
  return(vals0[0]);
}

double triang_clotoch(xev,vv,ce,p,xxa)
double *xev, *vv, *xxa;
int *ce, p;
{ double cfo[3], cfe[3], cg[9], *va, *vb, *vc,
    l0, nm[3], na, nb, nc, *xl, *xr, *xz, d0, d1, lb, dlt, gam;
  int i, w[3], cfl, cfr;
  if (p==1)
    return(xxa[0]*vv[0]+xxa[1]*vv[1]+xxa[2]*vv[2]);
  if (xxa[2]<=MIN(xxa[0],xxa[1]))
  { va = &xev[2*ce[0]]; vb = &xev[2*ce[1]]; vc = &xev[2*ce[2]];
    w[0] = 0; w[1] = 3; w[2] = 6;
  }
  else
  if (xxa[1]<xxa[0])
  { w[0] = 0; w[1] = 6; w[2] = 3;
    va = &xev[2*ce[0]]; vb = &xev[2*ce[2]]; vc = &xev[2*ce[1]];
    lb = xxa[1]; xxa[1] = xxa[2]; xxa[2] = lb;
  }
  else
  { w[0] = 6; w[1] = 3; w[2] = 0;
    va = &xev[2*ce[2]]; vb = &xev[2*ce[1]]; vc = &xev[2*ce[0]];
    lb = xxa[0]; xxa[0] = xxa[2]; xxa[2] = lb;
  }
  
/* set cg to values and derivatives on standard triangle */
  for (i=0; i<3; i++)
  { cg[3*i] = vv[w[i]];
    cg[3*i+1] = ((vb[0]-va[0])*vv[w[i]+1]
                +(vb[1]-va[1])*vv[w[i]+2])/2;  /* df/dx */
    cg[3*i+2] = ((2*vc[0]-vb[0]-va[0])*vv[w[i]+1]
                +(2*vc[1]-vb[1]-va[1])*vv[w[i]+2])/2.0; /* sqrt{3} df/dy */
  }
  dlt = (vb[0]-va[0])*(vc[1]-va[1])-(vc[0]-va[0])*(vb[1]-va[1]);
  /* Twice area; +ve if abc antic.wise  -ve is abc c.wise */
  cfo[0] = (cg[0]+cg[3]+cg[6])/3;
  cfo[1] = (2*cg[0]-cg[3]-cg[6])/4;
  cfo[2] = (2*cg[3]-cg[0]-cg[6])/4;
  na = -cg[1]+cg[2];  /* perp. deriv, rel. length 2 */
  nb = -cg[4]-cg[5];
  nc = 2*cg[7];
  cfo[1] += (nb-nc)/16;
  cfo[2] += (nc-na)/16;
  na = -cg[1]-cg[2]/3.0;  /* derivatives back to origin */
  nb =  cg[4]-cg[5]/3.0;
  nc =        cg[8]/1.5;
  cfo[0] -= (na+nb+nc)*7/54;
  cfo[1] += 13*(nb+nc-2*na)/144;
  cfo[2] += 13*(na+nc-2*nb)/144;
  for (i=0; i<3; i++)
  { /* Outward normals by linear interpolation on original triangle.
       Convert to outward normals on standard triangle.
       Actually, computed to opposite corner */
    switch(i)
    { case 0: xl = vc; xr = vb; xz = va; cfl = w[2]; cfr = w[1];
              break;
      case 1: xl = va; xr = vc; xz = vb; cfl = w[0]; cfr = w[2];
              break;
      case 2: xl = vb; xr = va; xz = vc; cfl = w[1]; cfr = w[0];
              break;
    }
    na = xr[0]-xl[0]; nb = xr[1]-xl[1];
    lb = na*na+nb*nb;
    d0 = 1.5*(vv[cfr]-vv[cfl]) - 0.25*(na*(vv[cfl+1]+vv[cfr+1])
         +nb*(vv[cfl+2]+vv[cfr+2]));
    d1 = 0.5*( na*(vv[cfl+2]+vv[cfr+2])-nb*(vv[cfl+1]+vv[cfr+1]) );
    l0 = (xz[0]-xl[0])*na+(xz[1]-xl[1])*nb-lb/2;
    nm[i] = (d1*dlt-l0*d0)/lb;
  }
  cfo[0] -= (nm[0]+nm[1]+nm[2])*4/81;
  cfo[1] += (2*nm[0]-nm[1]-nm[2])/27;
  cfo[2] += (2*nm[1]-nm[0]-nm[2])/27;

  gam = xxa[0]+xxa[1]-2*xxa[2];
  if (gam==0) return(cfo[0]);
  lb = (xxa[0]-xxa[2])/gam;
  d0 = -2*cg[4]; d1 = -2*cg[1];
  cfe[0] = cubic_interp(lb,cg[3],cg[0],d0,d1);
  cfe[1] = cubintd(lb,cg[3],cg[0],d0,d1);
  cfe[2] = -(1-lb)*(1-2*lb)*cg[5] + 4*lb*(1-lb)*nm[2] - lb*(2*lb-1)*cg[2];
  d0 = 2*(lb*cfo[1]+(1-lb)*cfo[2]);
  d1 = (lb-0.5)*cfe[1]+cfe[2]/3.0;
  return(cubic_interp(gam,cfo[0],cfe[0],d0,d1));
}

int triang_getvertexvals(fp,evs,vv,i,what)
fitpt *fp;
evstruc *evs;
double *vv;
int i, what;
{ double dx, P, le, vl[1+MXDIM], vh[1+MXDIM];
  int d, il, ih, j, nc;
  d = fp->d;
  if (evs->s[i]==0) return(exvval(fp,vv,i,d,what,0));

  il = evs->lo[i]; nc = triang_getvertexvals(fp,evs,vl,il,what);
  ih = evs->hi[i]; nc = triang_getvertexvals(fp,evs,vh,ih,what);
  vv[0] = (vl[0]+vh[0])/2;
  if (nc==1) return(nc);
  P = 1.5*(vh[0]-vl[0]);
  le = 0.0;
  for (j=0; j<d; j++)
  { dx = evptx(fp,ih,j)-evptx(fp,il,j);
    vv[0] += dx*(vl[j+1]-vh[j+1])/8;
    vv[j+1] = (vl[j+1]+vh[j+1])/2;
    P -= 1.5*dx*vv[j+1];
    le += dx*dx;
  }
  for (j=0; j<d; j++)
    vv[j+1] += P*(evptx(fp,ih,j)-evptx(fp,il,j))/le;
  return(nc);
}

double triang_int(lf,x,what)
lfit *lf;
double *x;
int what;
{
  int d, i, j, k, vc, nc;
  int *ce, nce[1+MXDIM];
  double xa[1+MXDIM], vv[(1+MXDIM)*(1+MXDIM)], lb;
fitpt *fp;
evstruc *evs;
fp = &lf->fp;
evs= &lf->evs;

  d = fp->d; vc = d+1;
  ce = evs->ce;
  i = 0;
  while ((i<evs->nce) && (!intri(x,&ce[i*vc],evp(fp),xa,d))) i++;
  if (i==evs->nce) return(NOSLN);
  i *= vc;
  for (j=0; j<vc; j++) nce[j] = ce[i+j];
  triang_descend(lf,xa,nce);

  /* order the vertices -- needed for asymmetric interptr */
  do
  { k=0;
    for (i=0; i<d; i++)
      if (nce[i]>nce[i+1])
      { j=nce[i]; nce[i]=nce[i+1]; nce[i+1]=j; k=1;
        lb = xa[i]; xa[i] = xa[i+1]; xa[i+1] = lb;
      }
  } while(k);
  nc = 0;
  for (i=0; i<vc; i++)
    nc =  triang_getvertexvals(fp,evs,&vv[i*nc],nce[i],what);
  return((d==2) ? triang_clotoch(evp(fp),vv,nce,nc,xa) :
                 triang_cubicint(evp(fp),vv,nce,d,nc,xa));
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
/*
 * Functions for computing residuals and fitted values from
 * the locfit object.
 *
 * fitted(lf,fit,what,cv,ty) computes fitted values from the
 *   fit structure in lf. 
 * resid(y,c,w,th,mi,ty) converts fitted values to residuals
*/

#include "lfev.h"

#define NRT 8
static char *rtype[NRT] = { "deviance", "d2",    "pearson", "raw",
                          "ldot",     "lddot", "fit",     "mean" };
static int   rvals[NRT] = { RDEV, RDEV2, RPEAR, RRAW, RLDOT, RLDDT, RFIT, RMEAN};
int restyp(z)
char *z;
{ int val;

  val = pmatch(z, rtype, rvals, NRT, -1);
  if (val==-1) LERR(("Unknown type = %s",z));
  return(val);
}

double resid(y,w,th,fam,ty,res)
int fam, ty;
double y, w, th, *res;
{ double raw;

  fam = fam & 63;
  if ((fam==TGAUS) | (fam==TROBT) | (fam==TCAUC))
    raw = y-res[ZMEAN];
  else
    raw = y-w*res[ZMEAN];
  switch(ty)
  { case RDEV:
      if (res[ZDLL]>0) return(sqrt(-2*res[ZLIK]));
            else return(-sqrt(-2*res[ZLIK]));
    case RPEAR:
      if (res[ZDDLL]<=0)
      { if (res[ZDLL]==0) return(0);
        return(NOSLN);
      }
      return(res[ZDLL]/sqrt(res[ZDDLL]));
    case RRAW:  return(raw);
    case RLDOT: return(res[ZDLL]);
    case RDEV2: return(-2*res[ZLIK]);
    case RLDDT: return(res[ZDDLL]);
    case RFIT:  return(th);
    case RMEAN: return(res[ZMEAN]);
    default: LERR(("resid: unknown residual type %d",ty));
  }
  return(0.0);
}

double studentize(res,inl,var,ty,link)
double res, inl, var, *link;
int ty;
{ double den;
  inl *= link[ZDDLL];
  var = var*var*link[ZDDLL];
  if (inl>1) inl = 1;
  if (var>inl) var = inl;
  den = 1-2*inl+var;
  if (den<0) return(0.0);
  switch(ty)
  { case RDEV:
    case RPEAR:
    case RRAW:
    case RLDOT:
      return(res/sqrt(den));
    case RDEV2:
      return(res/den);
    default: return(res);
  }
}

void fitted(lf,fit,what,cv,st,ty)
lfit *lf;
double *fit;
int what, cv, st, ty;
{ int i, j, d, n, evo;
  double xx[MXDIM], th, inl, var, link[LLEN];
  n = lf->lfd.n;
  d = lf->lfd.d;
  evo = ev(&lf->evs);
  cv &= (evo!=ECROS);
  if ((evo==EDATA)|(evo==ECROS)) evo = EFITP;
  setfamily(&lf->sp);

  for (i=0; i<n; i++)
  { for (j=0; j<d; j++) xx[j] = datum(&lf->lfd,j,i);
    th = dointpoint(lf,xx,what,evo,i);
    if ((what==PT0)|(what==PVARI)) th = th*th;
    if (what==PCOEF)
    {
      th += base(&lf->lfd,i);
      stdlinks(link,&lf->lfd,&lf->sp,i,th,rsc(&lf->fp));
      if ((cv)|(st))
      { inl = dointpoint(lf,xx,PT0,evo,i);
        inl = inl*inl;
        if (cv)
        { th -= inl*link[ZDLL];
          stdlinks(link,&lf->lfd,&lf->sp,i,th,rsc(&lf->fp));
        }
        if (st) var = dointpoint(lf,xx,PNLX,evo,i);
      }
      fit[i] = resid(resp(&lf->lfd,i),prwt(&lf->lfd,i),th,fam(&lf->sp),ty,link);
      if (st) fit[i] = studentize(fit[i],inl,var,ty,link);
    } else fit[i] = th;
    if (lf_error) return;
  }
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
#include "lfev.h"

extern double robscale;

/* special version of ressumm to estimate sigma^2, with derivative estimation */
void ressummd(lf)
lfit *lf;
{ int i;
  double s0, s1;
  s0 = s1 = 0.0;
  if ((fam(&lf->sp)&64)==0)
  { rv(&lf->fp) = 1.0;
    return;
  }
  for (i=0; i<lf->fp.nv; i++)
  { s0 += lf->fp.lik[2*lf->fp.nvm+i];
    s1 += lf->fp.lik[i];
  }
  if (s0==0.0)
    rv(&lf->fp) = 0.0;
  else
    rv(&lf->fp) = -2*s1/s0;
}

/*
 * res[0] = log-likelihood.
 * res[1] = df0 = tr(H)
 * res[2] = df0 = tr(H'H)
 * res[3] = s^2.
 * res[5] = robustscale.
 */
void ressumm(lf,des,res)
lfit *lf;
design *des;
double *res;
{ int i, j, evo, tg;
  double *oy, pw, r1, r2, rdf, t0, t1, u[MXDIM], link[LLEN];
  fitpt *fp;

  res[0] = res[1] = res[2] = 0.0;

  evo = ev(&lf->evs);
  if ((evo==EKDCE) | (evo==EPRES))
  { res[3] = 1.0;
    return;
  }
  if (lf->dv.nd>0)
  { ressummd(lf);
    return;
  }

  r1 = r2 = 0.0;
  if ((evo==EDATA) | (evo==ECROS)) evo = EFITP;

  for (i=0; i<lf->lfd.n; i++)
  { for (j=0; j<lf->lfd.d; j++) u[j] = datum(&lf->lfd,j,i);
    fitv(des,i) = base(&lf->lfd,i)+dointpoint(lf,u,PCOEF,evo,i);
    des->wd[i] = resp(&(lf->lfd),i) - fitv(des,i);
    wght(des,i) = 1.0;
    des->ind[i] = i;
  }

  tg = fam(&lf->sp);
  res[5] = 1.0;
  if ((tg==TROBT+64) | (tg==TCAUC+64)) /* global robust scale */
  { oy = lf->lfd.y; lf->lfd.y = des->wd;
    des->xev = lf->pc.xbar;
    locfit(&lf->lfd,des,&lf->sp,1,0,0);
    lf->lfd.y = oy;
    res[5] = robscale;
  }

  for (i=0; i<lf->lfd.n; i++)
  { for (j=0; j<lf->lfd.d; j++) u[j] = datum(&lf->lfd,j,i);
    t0 = dointpoint(lf,u,PT0,evo,i);
    t1 = dointpoint(lf,u,PNLX,evo,i);
    stdlinks(link,&lf->lfd,&lf->sp,i,fitv(des,i),res[5]);
    t1 = t1*t1*link[ZDDLL];
    t0 = t0*t0*link[ZDDLL];
    if (t1>1) t1 = 1;
    if (t0>1) t0 = 1; /* no observation gives >1 deg.free */
    res[0] += link[ZLIK];
    res[1] += t0;
    res[2] += t1;
    pw = prwt(&lf->lfd,i);
    if (pw>0)
    { r1 += link[ZDLL]*link[ZDLL]/pw;
      r2 += link[ZDDLL]/pw;
    }
  }

  res[3] = 1.0;
  if ((fam(&lf->sp)&64)==64) /* quasi family */
  { rdf = lf->lfd.n-2*res[1]+res[2];
    if (rdf<1.0)
    { WARN(("Estimated rdf < 1.0; not estimating variance"));
    }
    else
      res[3] = r1/r2 * lf->lfd.n / rdf;
  }

  /* try to ensure consistency for family="circ"! */
  if (((fam(&lf->sp)&63)==TCIRC) & (lf->lfd.d==1))
  { int *ind, nv;
    double dlt, th0, th1;
    ind = des->ind;
    nv = fp->nv;
    for (i=0; i<nv; i++) ind[i] = i;
    lforder(ind,evp(fp),0,nv-1);
    for (i=1; i<nv; i++)
    { dlt = evptx(fp,ind[i],0)-evptx(fp,ind[i-1],0);
      th0 = fp->coef[ind[i]]-dlt*fp->coef[ind[i]+nv]-fp->coef[ind[i-1]];
      th1 = fp->coef[ind[i]]-dlt*fp->coef[ind[i-1]+nv]-fp->coef[ind[i-1]];
      if ((th0>PI)&(th1>PI))
      { for (j=0; j<i; j++)
          fp->coef[ind[j]] += 2*PI;
        i--;
      }
      if ((th0<(-PI))&(th1<(-PI)))
      { for (j=0; j<i; j++)
          fp->coef[ind[j]] -= 2*PI;
        i--;
      }
    }
  }

}

double rss(lf,des,df)
lfit *lf;
design *des;
double *df;
{ double ss, res[10];
  ss = 0;
  ressumm(lf,des,res);
  *df = lf->lfd.n - 2*res[1] + res[2];
  return(-2*res[0]);
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
/*
 *
 *   Derivative corrections. The local slopes are not the derivatives
 *   of the local likelihood estimate; the function dercor() computes
 *   the adjustment to get the correct derivatives under the assumption
 *   that h is constant.
 *
 *   By differentiating the local likelihood equations, one obtains
 *
 *     d ^      ^       T      -1   T  d    .       ^
 *    -- a   =  a  -  (X W V X)    X  -- W  l( Y, X a)
 *    dx  0      1                    dx
 */

#include "lfev.h"
extern double robscale;

void dercor(lfd,sp,des,coef)
lfdata *lfd;
smpar *sp;
design *des;
double *coef;
{ double s1, dc[MXDIM], wd, link[LLEN];
  int i, ii, j, m, d, p;
  if (fam(sp)<=THAZ) return;
  if (ker(sp)==WPARM) return;

  d = lfd->d;
  p = des->p; m = des->n;

  if (lf_debug>1) mut_printf("  Correcting derivatives\n");
  fitfun(lfd, sp, des->xev,des->xev,des->f1,NULL);
  jacob_solve(&des->xtwx,des->f1);
  setzero(dc,d);

  /* correction term is e1^T (XTWVX)^{-1} XTW' ldot. */
  for (i=0; i<m; i++)
  { s1 = innerprod(des->f1,d_xi(des,i),p);
    ii = des->ind[i];
    stdlinks(link,lfd,sp,ii,fitv(des,ii),robscale);
    for (j=0; j<d; j++)
    { wd = wght(des,ii)*weightd(datum(lfd,j,ii)-des->xev[j],lfd->sca[j],
        d,ker(sp),kt(sp),des->h,lfd->sty[j],dist(des,ii));
      dc[j] += s1*wd*link[ZDLL];
    }

  }
  for (j=0; j<d; j++) coef[j+1] += dc[j];
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
#include "lfev.h"

void allocallcf(lf)
lfit *lf;
{ lf->fp.coef = VVEC(lf,0);
  lf->fp.h    = VVEC(lf,NPAR(lf));
}

int procvallcf(des,lf,v)
design *des;
lfit *lf;
int v;
{ int i, p, lf_status;

  lf_status = procv_nov(des,lf,v);
  if (lf_error) return(lf_status);

  p = NPAR(lf);
  for (i=0; i<p; i++) VVAL(lf,v,i) = des->cf[i];
  lf->fp.h[v] = des->h;

  return(0);
}

void initallcf(lf)
lfit *lf;
{ PROCV(lf) = procvallcf;
  ALLOC(lf) = allocallcf;
  PPROC(lf) = NULL;
  KEEPV(lf) = NPAR(lf)+1;
  KEEPC(lf) = 0;
  NOPC(lf)  = 1;
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
#include "lfev.h"

void pprocgam(lf,des,res)
lfit *lf;
design *des;
double *res;
{ int i, j, n, evo, op;
  double *fv, *vr, df, t0, t1, u[MXDIM], link[LLEN];

  n = lf->lfd.n;
  fv = res;
  vr = &res[n];
  df = 0.0;

  evo = ev(&lf->evs);
  if (evo==EDATA) evo = EFITP;
  
  for (i=0; i<n; i++)
  { for (j=0; j<lf->lfd.d; j++)  u[j] = datum(&lf->lfd,j,i);
    fitv(des,i) = base(&lf->lfd,i)+dointpoint(lf,u,PCOEF,evo,i);
    lf->lfd.y[i] -= fitv(des,i);
    wght(des,i) = 1.0;
    des->ind[i] = i;

    t0 = dointpoint(lf,u,PT0,evo,i);
    t1 = dointpoint(lf,u,PNLX,evo,i);
    stdlinks(link,&lf->lfd,&lf->sp,i,fitv(des,i),1.0);
    t0 = t0*t0*link[ZDDLL];
    t1 = t1*t1*link[ZDDLL];
    if (t0>1) t0 = 1; /* no observation gives >1 deg.free */
    if (t1>1) t1 = 1; /* no observation gives >1 deg.free */
    vr[i] = t1;
    df += t0;
  }

  des->n = n;
  deg(&lf->sp) = 1;
  op = npar(&lf->sp);
  npar(&lf->sp) = des->p = 1+lf->lfd.d;
  des->xev = lf->pc.xbar;
  locfit(&lf->lfd,des,&lf->sp,1,0,0);

  for (i=0; i<n; i++) fv[i] = resp(&lf->lfd,i) - fitv(des,i);
  for (i=0; i<=lf->lfd.d; i++)
    lf->pc.coef[i] += des->cf[i];
  res[2*n] = df-2;
  npar(&lf->sp) = op;
}

void initgam(lf)
lfit *lf;
{ initstd(lf);
  PPROC(lf) = pprocgam;
  KEEPC(lf) = 2*NOBS(lf)+1;
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
#include "lfev.h"

int procvhatm(des,lf,v)
design *des;
lfit *lf;
int v;
{ int k;
  double *l;
  l = &lf->fp.coef[v*lf->lfd.n];
  if ((ker(&lf->sp)!=WPARM) | (!haspc(&lf->pc)))
  { k = procv_nov(des,lf,v);
    wdiag(&lf->lfd,&lf->sp,des,l,&lf->dv,0,1,1);
  }
  else
    wdiagp(&lf->lfd,&lf->sp,des,l,&lf->pc,&lf->dv,0,1,1);
  lf->fp.h[v] = des->h;
  return(k);
}

void allochatm(lf)
lfit *lf;
{ lf->fp.coef = VVEC(lf,0);
  lf->fp.h    = VVEC(lf,NOBS(lf));
}

void pprochatm(lf,des,res)
lfit *lf;
design *des;
double *res;
{ transpose(lf->fp.coef,lf->fp.nvm,lf->lfd.n);
}

void inithatm(lf)
lfit *lf;
{ PROCV(lf) = procvhatm;
  ALLOC(lf) = allochatm;
  PPROC(lf) = pprochatm;
  KEEPV(lf) = NOBS(lf)+1;
  KEEPC(lf) = 1;
  NOPC(lf) = 1; /* could use pc if mi[MKER] == WPARM */
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
#include "lfev.h"

static lfit *lf_scb;
static lfdata *lfd_scb;
static smpar  *scb_sp;
static design *des_scb;

int scbfitter(x,l,reqd)
double *x, *l;
int reqd;
{
  int m;
  des_scb->xev = x;
  if ((ker(scb_sp)!=WPARM) | (!haspc(&lf_scb->pc)))
  { locfit(lfd_scb,des_scb,&lf_scb->sp,1,1,0);
    m = wdiag(lfd_scb, scb_sp, des_scb,l,&lf_scb->dv,reqd,2,0);
  }
  else
    m = wdiagp(lfd_scb, scb_sp, des_scb,l,&lf_scb->pc,&lf_scb->dv,reqd,2,0);
  return(m);
}

int constants(lf,des,res)
lfit *lf;
design *des;
double *res;
{
  int d, m, nt;
  double *wk;
  evstruc *evs;

  lf_scb = lf;
  des_scb = des;
  lfd_scb = &lf->lfd;
  scb_sp  = &lf->sp;

  evs = &lf->evs;
  d = lfd_scb->d;
  m = lfd_scb->n;
  trchck(lf,0,0,0);

  if (lf_error) return(0);
  if ((ker(scb_sp) != WPARM) && (lf->sp.nn>0))
    WARN(("constants are approximate for varying h"));
  npar(scb_sp) = calcp(scb_sp,lf->lfd.d);
  des_init(des,m,npar(scb_sp));
  set_scales(&lf->lfd);
  set_flim(&lf->lfd,&lf->evs);
  compparcomp(des,&lf->lfd,&lf->sp,&lf->pc,ker(scb_sp)!=WPARM);
  
  wk = &res[d+1];
  nt = tube_constants(scbfitter,d,m,ev(evs),mg(evs),evs->fl,
    res,wk,(d>3) ? 4 : d+1,0);
  lf->evs.nce = nt;   /* cheat way to return it to S. */
  return(nt);
}

void initkappa(lf)
lfit *lf;
{ PROCV(lf) = NULL;
  ALLOC(lf) = NULL;
  PPROC(lf) = (void *)constants;
  KEEPV(lf) = 0;
  KEEPC(lf) = NVAR(lf)+1+k0_reqd(NVAR(lf),NOBS(lf),0);
  NOPC(lf) = 0;
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
#include "lfev.h"

/*  fix df computation for link=IDENT. */
void pproclscv(lf,des,res)
lfit *lf;
design *des;
double *res;
{ double df, fh, fh_cv, infl, z0, z1, x[MXDIM];
  int i, n, j, evo;
  z1 = df = 0.0;
  evo = ev(&lf->evs);
  n = lf->lfd.n;
  if ((evo==EDATA) | (evo==ECROS)) evo = EFITP;

  z0 = dens_integrate(lf,des,2);

  for (i=0; i<n; i++)
  { for (j=0; j<lf->lfd.d; j++) x[j] = datum(&lf->lfd,j,i);
    fh = base(&lf->lfd,i)+dointpoint(lf,x,PCOEF,evo,i);
    if (link(&lf->sp)==LLOG) fh = exp(fh);
    infl = dointpoint(lf,x,PT0,evo,i);
    infl = infl * infl;
    if (infl>1) infl = 1;
    fh_cv = (link(&lf->sp) == LIDENT) ?
       (n*fh - infl) / (n-1.0) : fh*(1-infl)*n/(n-1.0);
    z1 += fh_cv;
    df += infl;
  }

  res[0] = z0-2*z1/n;
  res[1] = df;
}

void initlscv(lf)
lfit *lf;
{ initstd(lf);
  KEEPC(lf) = 2;
  PPROC(lf) = pproclscv;
  NOPC(lf) = 1;
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
#include "lfev.h"

static double pen, sig2;

void goldensec(f,des,tr,eps,xm,ym,meth)
double (*f)(), eps, *xm, *ym;
int meth;
design *des;
lfit *tr;
{ double x[4], y[4], xx[11], yy[11];
  int i, im;
  xx[0] = tr->sp.fixh;
  if (xx[0]<=0)
  { LERR(("regband: initialize h>0"));
    return;
  }
  for (i=0; i<=10; i++)
  { if (i>0) xx[i] = (1+gold_rat)*xx[i-1];
    yy[i] = f(xx[i],des,tr,meth);
    if ((i==0) || (yy[i]<yy[im])) im = i;
  }
  if (im==0) im = 1;
  if (im==10)im = 9;
  x[0] = xx[im-1]; y[0] = yy[im-1];
  x[1] = xx[im];   y[1] = yy[im];
  x[3] = xx[im+1]; y[3] = yy[im+1];
  x[2] = gold_rat*x[3]+(1-gold_rat)*x[0];
  y[2] = f(x[2],des,tr,meth);
  while (x[3]-x[0]>eps)
  { if (y[1]<y[2])
    { x[3] = x[2]; y[3] = y[2];
      x[2] = x[1]; y[2] = y[1];
      x[1] = gold_rat*x[0]+(1-gold_rat)*x[3];
      y[1] = f(x[1],des,tr,meth);
    }
    else
    { x[0] = x[1]; y[0] = y[1];
      x[1] = x[2]; y[1] = y[2];
      x[2] = gold_rat*x[3]+(1-gold_rat)*x[0];
      y[2] = f(x[2],des,tr,meth);
    }
  }
  im = 0;
  for (i=1; i<4; i++) if (y[i]<y[im]) im = i;
  *xm = x[im]; *ym = y[im];
}

double dnk(x,k)
double x;
int k;
{ double f;
  switch(k)
  { case 0: f = 1; break;
    case 1: f = -x; break;
    case 2: f = x*x-1; break;
    case 3: f = x*(x*x-3); break;
    case 4: f = 3-x*x*(6-x*x); break;
    case 5: f = -x*(15-x*x*(10-x*x)); break;
    case 6: f = -15+x*x*(45-x*x*(15-x*x)); break;
    default: LERR(("dnk: k=%d too large",k)); return(0.0);
  }
  return(f*exp(-x*x/2)/S2PI);
}

double locai(h,des,lf)
double h;
design *des;
lfit *lf;
{ double cp, res[10];
  nn(&lf->sp) = h;
  lf->mdl.procv = procvstd;
  nstartlf(des,lf);
  ressumm(lf,des,res);
  cp = -2*res[0] + pen*res[1];
  return(cp);
}

static int fc;

double loccp(h,des,lf,m) /* m=1: cp    m=2: gcv */
double h;
design *des;
lfit *lf;
int m;
{ double cp, res[10];
  int dg, n;

  n = lf->lfd.n;
  nn(&lf->sp) = 0;
  fixh(&lf->sp) = h;
  lf->mdl.procv = procvstd;
  nstartlf(des,lf);
  ressumm(lf,des,res);
  if (m==1)
  { if (fc) sig2 = res[3]; /* first call - set sig2 */
    cp = -2*res[0]/sig2 - n + 2*res[1];
  }
  else
    cp = -2*n*res[0]/((n-res[1])*(n-res[1]));
  fc = 0;
  return(cp);
}

double cp(des,lf,meth)
design *des;
lfit *lf;
int meth;
{ double hm, ym;
  fc = 1;
  goldensec(loccp,des,lf,0.001,&hm,&ym,meth);
  return(hm);
}

double gkk(des,lf)
design *des;
lfit *lf;
{ double h, h5, nf, th;
  int i, j, n, dg0, dg1;
  ev(&lf->evs)= EDATA;
  nn(&lf->sp) = 0;
  n = lf->lfd.n;
  dg0 = deg0(&lf->sp);     /* target degree */
  dg1 = dg0+1+(dg0%2==0);  /* pilot degree */
  nf = exp(log(1.0*n)/10); /* bandwidth inflation factor */
  h = lf->sp.fixh;         /* start bandwidth */
  for (i=0; i<=10; i++)
  { deg(&lf->sp) = dg1;
    lf->sp.fixh = h*nf;
    lf->mdl.procv = procvstd;
    nstartlf(des,lf);
    th = 0;
    for (j=10; j<n-10; j++)
      th += lf->fp.coef[dg1*n+j]*lf->fp.coef[dg1*n+j];
th *= n/(n-20.0);
    h5 = sig2 * Wikk(ker(&lf->sp),dg0) / th;
    h = exp(log(h5)/(2*dg1+1));
    if (lf_error) return(0.0);
/* mut_printf("pilot %8.5f  sel %8.5f\n",lf->sp.fixh,h); */
  }
  return(h);
}

double rsw(des,lf)
design *des;
lfit *lf;
{ int i, j, k, nmax, nvm, n, mk, evo, dg0, dg1;
  double rss[6], cp[6], th22, dx, d2, hh;
  nmax = 5;
  evo = ev(&lf->evs); ev(&lf->evs) = EGRID;
  mk = ker(&lf->sp);  ker(&lf->sp) = WRECT;
  dg0 = deg0(&lf->sp);
  dg1 = 1 + dg0 + (dg0%2==0);
  deg(&lf->sp) = 4;
  for (k=nmax; k>0; k--)
  { lf->evs.mg[0] = k;
    lf->evs.fl[0] = 1.0/(2*k);
    lf->evs.fl[1] = 1-1.0/(2*k);
    nn(&lf->sp) = 0;
    fixh(&lf->sp) = 1.0/(2*k);
    lf->mdl.procv = procvstd;
    nstartlf(des,lf);
    nvm = lf->fp.nvm;
    rss[k] = 0;
    for (i=0; i<k; i++) rss[k] += -2*lf->fp.lik[i];
  }
  n = lf->lfd.n; k = 1;
  for (i=1; i<=nmax; i++)
  { /* cp[i] = (n-5*nmax)*rss[i]/rss[nmax]-(n-10*i); */
    cp[i] = rss[i]/sig2-(n-10*i);
    if (cp[i]<cp[k]) k = i;
  }
  lf->evs.mg[0] = k;
  lf->evs.fl[0] = 1.0/(2*k);
  lf->evs.fl[1] = 1-1.0/(2*k);
  nn(&lf->sp) = 0;
  fixh(&lf->sp) = 1.0/(2*k);
  lf->mdl.procv = procvstd;
  nstartlf(des,lf);
  ker(&lf->sp) = mk; ev(&lf->evs) = evo;
  nvm = lf->fp.nvm;
  th22 = 0;
  for (i=10; i<n-10; i++)
  { j = floor(k*datum(&lf->lfd,0,i));
    if (j>=k) j = k-1;
    dx = datum(&lf->lfd,0,i)-evptx(&lf->fp,0,j);
    if (dg1==2)
      d2 = lf->fp.coef[2*nvm+j]+dx*lf->fp.coef[3*nvm+j]+dx*dx*lf->fp.coef[4*nvm+j]/2;
    else d2 = lf->fp.coef[4*nvm+j];
    th22 += d2*d2;
  }
  hh = Wikk(mk,dg0)*sig2/th22*(n-20.0)/n;
  return(exp(log(hh)/(2*dg1+1)));
}

#define MAXMETH 10

void rband(lf,des,hhat)
lfit *lf;
design *des;
double *hhat;
{ int i, dg, nmeth, meth[MAXMETH];
  double h0, res[10];

  nmeth = lf->dv.nd;
  for (i=0; i<nmeth; i++) meth[i] = lf->dv.deriv[i];
  lf->dv.nd = 0;

/* first, estimate sigma^2.
 * this is ridiculously bad. lf->sp.fixh = 0.05????
 */
/*  dg = deg(&lf->sp); deg(&lf->sp) = 2;
  h0 = lf->sp.fixh;  lf->sp.fixh = 0.05;
  lf->mdl.procv = procvstd;
  nstartlf(des,lf);
  ressumm(lf,des,res);
  deg(&lf->sp) = dg; lf->sp.fixh = h0;
  sig2 = res[3];  */

  for (i=0; i<nmeth; i++)
  {
    switch(meth[i])
    { case 0: hhat[i] = cp(des,lf,1);
              break;
      case 1: hhat[i] = cp(des,lf,2);
              break;
      case 2: hhat[i] = gkk(des,lf);
              break;
      case 3: hhat[i] = rsw(des,lf);
              break;
      default: hhat[i] = 0;
              mut_printf("Unknown method %d\n",meth[i]);
    }
    if (lf_error) i = nmeth;
  }
  lf->dv.nd = nmeth;
}

void initrband(lf)
lfit *lf;
{
  initstd(lf);
  KEEPC(lf) = MAXMETH;
  PROCV(lf) = NULL;
  PPROC(lf) = rband;
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
#include "lfev.h"
static double scb_crit, *x, c[10], kap[5], kaq[5], max_p2;
static int side, type;
design *scb_des;

double covar_par(lf,des,x1,x2)
lfit *lf;
design *des;
double x1, x2;
{ double *v1, *v2, *wk;
  paramcomp *pc;
  int i, j, p, ispar;

  v1 = des->f1; v2 = des->ss; wk = des->oc;
  ispar = (ker(&lf->sp)==WPARM) && (haspc(&lf->pc));
  p = npar(&lf->sp);

/*  for parametric models, the covariance is
 *  A(x1)^T (X^T W V X)^{-1} A(x2)
 *  which we can find easily from the parametric component.
 */
  if (ispar)
  { pc = &lf->pc;
    fitfun(&lf->lfd, &lf->sp, &x1,pc->xbar,v1,NULL);
    fitfun(&lf->lfd, &lf->sp, &x2,pc->xbar,v2,NULL);
    jacob_hsolve(&lf->pc.xtwx,v1);
    jacob_hsolve(&lf->pc.xtwx,v2);
  }

/*  for non-parametric models, we must use the cholseky decomposition
 *  of M2 = X^T W^2 V X. Courtesy of lf_vcov caclulations, should have
 *  des->P = M2^{1/2} M1^{-1}.
 */
  if (!ispar)
  { fitfun(&lf->lfd, &lf->sp, &x1,des->xev,wk,NULL);
    for (i=0; i<p; i++)
    { v1[i] = 0;
      for (j=0; j<p; j++) v1[i] += des->P[i*p+j]*wk[j];
    }
    fitfun(&lf->lfd, &lf->sp, &x2,des->xev,wk,NULL);
    for (i=0; i<p; i++)
    { v2[i] = 0;
      for (j=0; j<p; j++) v2[i] += des->P[i*p+j]*wk[j];
    }
  }

  return(innerprod(v1,v2,p));
}

void cumulant(lf,des,sd)
lfit *lf;
design *des;
double sd;
{ double b2i, b3i, b3j, b4i;
  double ss, si, sj, uii, uij, ujj, k1;
  int ii, i, j, jj;
  for (i=1; i<10; i++) c[i] = 0.0;
  k1 = 0;

  /* ss = sd*sd; */
  ss = covar_par(lf,des,des->xev[0],des->xev[0]);

/*
 * this isn't valid for nonparametric models. At a minimum,
 * the sums would have to include weights. Still have to work
 * out the right way.
 */
  for (i=0; i<lf->lfd.n; i++)
  { ii = des->ind[i];
    b2i = b2(fitv(des,ii),fam(&lf->sp),prwt(&lf->lfd,ii));
    b3i = b3(fitv(des,ii),fam(&lf->sp),prwt(&lf->lfd,ii));
    b4i = b4(fitv(des,ii),fam(&lf->sp),prwt(&lf->lfd,ii));
    si = covar_par(lf,des,des->xev[0],datum(&lf->lfd,0,ii));
    uii= covar_par(lf,des,datum(&lf->lfd,0,ii),datum(&lf->lfd,0,ii));
    if (lf_error) return;

    c[2] += b4i*si*si*uii;
    c[6] += b4i*si*si*si*si;
    c[7] += b3i*si*uii;
    c[8] += b3i*si*si*si;
    /* c[9] += b2i*si*si*si*si;
       c[9] += b2i*b2i*si*si*si*si; */
    k1 += b3i*si*(si*si/ss-uii);

    /* i=j components */
    c[1] += b3i*b3i*si*si*uii*uii;
    c[3] += b3i*b3i*si*si*si*si*uii;
    c[4] += b3i*b3i*si*si*uii*uii;

    for (j=i+1; j<lf->lfd.n; j++)
    { jj = des->ind[j];
      b3j = b3(fitv(des,ii),fam(&lf->sp),prwt(&lf->lfd,jj));
      sj = covar_par(lf,des,des->xev[0],datum(&lf->lfd,0,jj));
      uij= covar_par(lf,des,datum(&lf->lfd,0,ii),datum(&lf->lfd,0,jj));
      ujj= covar_par(lf,des,datum(&lf->lfd,0,jj),datum(&lf->lfd,0,jj));

      c[1] += 2*b3i*b3j*si*sj*uij*uij;
      c[3] += 2*b3i*b3j*si*si*sj*sj*uij;
      c[4] += b3i*b3j*uij*(si*si*ujj+sj*sj*uii);
      if (lf_error) return;
    }
  }
  c[5] = c[1];
  c[7] = c[7]*c[8];
  c[8] = c[8]*c[8];

  c[1] /= ss; c[2] /= ss; c[3] /= ss*ss; c[4] /= ss;
  c[5] /= ss; c[6] /= ss*ss; c[7] /= ss*ss;
  c[8] /= ss*ss*ss; c[9] /= ss*ss;

/* constants used in p(x,z) computation */
  kap[1] = k1/(2*sqrt(ss));
  kap[2] = 1 + 0.5*(c[1]-c[2]+c[4]-c[7]) - 3*c[3] + c[6] + 1.75*c[8];
  kap[4] = -9*c[3] + 3*c[6] + 6*c[8] + 3*c[9];

/* constants used in q(x,u) computation */
  kaq[2] = c[3] - 1.5*c[8] - c[5] - c[4] + 0.5*c[7] + c[6] - c[2];
  kaq[4] = -3*c[3] - 6*c[4] - 6*c[5] + 3*c[6] + 3*c[7] - 3*c[8] + 3*c[9];
}

/* q2(u) := u+q2(x,u) in paper */
double q2(u)
double u;
{ return(u-u*(36.0*kaq[2] + 3*kaq[4]*(u*u-3) + c[8]*((u*u-10)*u*u+15))/72.0);
}

/*  p2(u) := p2(x,u) in paper */
double p2(u)
double u;
{ return( -u*( 36*(kap[2]-1+kap[1]*kap[1])
     + 3*(kap[4]+4*kap[1]*sqrt(kap[3]))*(u*u-3)
     + c[8]*((u*u-10)*u*u+15) ) / 72 );
}

extern int likereg();
double gldn_like(a)
double a;
{ int err;

  scb_des->fix[0] = 1;
  scb_des->cf[0] = a;
  max_nr(likereg, scb_des->cf, scb_des->oc, scb_des->res, scb_des->f1,
    &scb_des->xtwx, scb_des->p, lf_maxit, 1.0e-6, &err); 
  scb_des->fix[0] = 0;

  return(scb_des->llk);
}

/* v1/v2 is correct for deg=0 only */
void get_gldn(fp,des,lo,hi,v)
fitpt *fp;
design *des;
double *lo, *hi;
int v;
{ double v1, v2, c, tlk;
  int err;

  v1 = fp->nlx[v];
  v2 = fp->t0[v];
  c = scb_crit * v1 / v2;
  tlk = des->llk - c*c/2;
mut_printf("v %8.5f %8.5f  c %8.5f  tlk %8.5f   llk %8.5f\n",v1,v2,c,tlk,des->llk);

  /* want: { a : l(a) >= l(a-hat) - c*c/2 } */
  lo[v] = fp->coef[v] - scb_crit*v1;
  hi[v] = fp->coef[v] + scb_crit*v1;

  err = 0;

mut_printf("lo %2d\n",v);
  lo[v] = solve_secant(gldn_like,tlk,lo[v],fp->coef[v],1e-8,BDF_EXPLEFT,&err);
  if (err>0) mut_printf("solve_secant error\n");
mut_printf("hi %2d\n",v);
  hi[v] = solve_secant(gldn_like,tlk,fp->coef[v],hi[v],1e-8,BDF_EXPRIGHT,&err);
  if (err>0) mut_printf("solve_secant error\n");
}

int procvscb2(des,lf,v)
design *des;
lfit *lf;
int v;
{ double thhat, sd, *lo, *hi, u;
  int err, st, tmp;
  x = des->xev = evpt(&lf->fp,v);
  tmp = haspc(&lf->pc);
  /* if ((ker(&lf->sp)==WPARM) && (haspc(&lf->pc)))
  { lf->coef[v] = thhat = addparcomp(lf,des->xev,PCOEF);
    lf->nlx[v] = lf->t0[v] = sd = addparcomp(lf,des->xev,PNLX);
  }
  else */
  { haspc(&lf->pc) = 0;
    st = procvstd(des,lf,v);
    thhat = lf->fp.coef[v];
    sd = lf->fp.nlx[v];
  }
  if ((type==GLM2) | (type==GLM3) | (type==GLM4))
  { if (ker(&lf->sp) != WPARM)
      WARN(("nonparametric fit; correction is invalid"));
    cumulant(lf,des,sd);
  }
  haspc(&lf->pc) = tmp;
  lo = lf->fp.t0;
  hi = &lo[lf->fp.nvm];
  switch(type)
  {
    case GLM1:
      return(st);
    case GLM2: /* centered scr */
      lo[v] = kap[1];
      hi[v] = sqrt(kap[2]);
      return(st);
    case GLM3: /* corrected 2 */
      lo[v] = solve_secant(q2,scb_crit,0.0,2*scb_crit,0.000001,BDF_NONE,&err);
      return(st);
    case GLM4: /* corrected 2' */
      u = fabs(p2(scb_crit));
      max_p2 = MAX(max_p2,u);
      return(st);
    case GLDN:
      get_gldn(&lf->fp,des,lo,hi,v);
      return(st);
  }
  LERR(("procvscb2: invalid type"));
  return(st);
}

void scb(lf,des,res)
lfit *lf;
design *des;
double *res;
{ double k1, k2, *lo, *hi, sig, thhat, nlx, rss[10];
  int i, nterms;

  scb_des= des;

  npar(&lf->sp) = calcp(&lf->sp,lf->lfd.d);
  des_init(des,lf->lfd.n,npar(&lf->sp));

  type = geth(&lf->fp);

  if (type >= 80) /* simultaneous */
  {
    nterms = constants(lf,des,res);
    scb_crit = critval(0.05,res,nterms,lf->lfd.d,TWO_SIDED,0.0,GAUSS);
    type -= 10;
  }
  else /* pointwise */
  { res[0] = 1;
    scb_crit = critval(0.05,res,1,lf->lfd.d,TWO_SIDED,0.0,GAUSS);
  }

  max_p2 = 0.0;
  lf->mdl.procv = procvscb2;
  nstartlf(des,lf);

  if ((fam(&lf->sp)&64)==64)
  { i = haspc(&lf->pc); haspc(&lf->pc) = 0;
    ressumm(lf,des,rss);
    haspc(&lf->pc) = i;
    sig = sqrt(rss[3]);
  }
  else sig = 1.0;

  lo = lf->fp.t0;
  hi = &lo[lf->fp.nvm];
  for (i=0; i<lf->fp.nv; i++)
  { thhat = lf->fp.coef[i];
    nlx = lf->fp.nlx[i];
    switch(type)
    {
      case GLM1:  /* basic scb */
        lo[i] = thhat - scb_crit * sig * nlx;
        hi[i] = thhat + scb_crit * sig * nlx;
        break;
      case GLM2:
        k1 = lo[i];
        k2 = hi[i];
        lo[i] = thhat - k1*nlx - scb_crit*nlx*k2;
        hi[i] = thhat - k1*nlx + scb_crit*nlx*k2;
        break;
      case GLM3:
        k1 = lo[i];
        lo[i] = thhat - k1*nlx;
        hi[i] = thhat + k1*nlx;
      case GLM4:  /* corrected 2' */
        lo[i] = thhat - (scb_crit-max_p2)*lf->fp.nlx[i];
        hi[i] = thhat + (scb_crit-max_p2)*lf->fp.nlx[i];
        break;
      case GLDN:
        break;
    }
  }
}

void initscb(lf)
lfit *lf;
{ initstd(lf);
  PROCV(lf) = NULL;
  KEEPC(lf) = NVAR(lf)+1+k0_reqd(NVAR(lf),NOBS(lf),0);
  PPROC(lf) = scb;
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
#include "lfev.h"

int procvsimple(des,lf,v)
design *des;
lfit *lf;
int v;
{ int lf_status;
  lf_status = procv_nov(des,lf,v);
  VVAL(lf,v,0) = des->cf[cfn(des,0)];
  return(lf_status);
}

void allocsimple(lf)
lfit *lf;
{ lf->fp.coef = VVEC(lf,0);
  lf->fp.h = NULL;
}

void initsimple(lf)
lfit *lf;
{
  PROCV(lf) = procvsimple;
  ALLOC(lf) = allocsimple;
  PPROC(lf) = NULL;
  KEEPV(lf) = 1;
  KEEPC(lf) = 1;
  NOPC(lf)  = 1;
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
#include "lfev.h"

/* 3d+8 variables to keep:
 * d+1 coef+derivs.
 * d+1 sd's + derivs.
 * d+1 infl + derivs.
 *   3 likelihood and d.f's.
 *   1 bandwidth h
 *   1 degree.
 */

void allocstd(lf)
lfit *lf;
{ int d;
  d = NVAR(lf);
  lf->fp.coef = VVEC(lf,0);
  lf->fp.nlx  = VVEC(lf,d+1);
  lf->fp.t0   = VVEC(lf,2*d+2);
  lf->fp.lik  = VVEC(lf,3*d+3);
  lf->fp.h    = VVEC(lf,3*d+6);
  lf->fp.deg  = VVEC(lf,3*d+7);
}

int procvstd(des,lf,v)
design *des;
lfit *lf;
int v;
{ int d, p, nvm, i, k;
  double t0[1+MXDIM], vari[1+MXDIM];
  k = procv_var(des,lf,v);
  if (lf_error) return(k);
   
  d = lf->lfd.d;
  p = npar(&lf->sp);
  nvm = lf->fp.nvm;

  if (k != LF_OK) lf_status_msg(k);

  lf->fp.lik[v] = des->llk;
  lf->fp.lik[nvm+v] = des->tr2;
  lf->fp.lik[2*nvm+v] = des->tr0 - des->tr2;

  for (i=0; i<des->ncoef; i++)
    vari[i] = des->V[p*cfn(des,0) + cfn(des,i)];
  vari[0] = sqrt(vari[0]);
  if (vari[0]>0) for (i=1; i<des->ncoef; i++) vari[i] /= vari[0];

  t0[0] = sqrt(des->f1[0]);
  if (t0[0]>0) for (i=1; i<des->ncoef; i++) t0[i] = des->f1[i]/t0[0];

  if (dc(&lf->fp)) dercor(&lf->lfd,&lf->sp,des,des->cf);
  subparcomp(des,lf,des->cf);
  for (i=0; i<des->ncoef; i++)
    lf->fp.coef[i*lf->fp.nvm+v] =  des->cf[cfn(des,i)];

  subparcomp2(des,lf,vari,t0);
  for (i=0; i<des->ncoef; i++)
  { lf->fp.nlx[i*nvm+v] = vari[i];
    lf->fp.t0[i*nvm+v]  = t0[i];
  }

  lf->fp.deg[v] = deg(&lf->sp);
  return(k);
}

void pprocstd(lf,des,res)
lfit *lf;
design *des;
double *res;
{
  ressumm(lf,des,res);
}

void initstd(lf)
lfit *lf;
{ PROCV(lf) = procvstd;
  ALLOC(lf) = allocstd;
  PPROC(lf) = pprocstd;
  KEEPV(lf) = 3*NVAR(lf) + 8;
  KEEPC(lf) = 6;
  NOPC(lf)  = 0;
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
#include "lfev.h"

extern void procstd(), allocstd();
extern double robscale;

double vocri(lk,t0,t2,pen)
double lk, t0, t2, pen;
{ if (pen==0) return(-2*t0*lk/((t0-t2)*(t0-t2)));
  return((-2*lk+pen*t2)/t0);
}

double intvo(des,lf,c0,c1,a,p,t0,t20,t21)
design *des;
lfit *lf;
double *c0, *c1, a, t0, t20, t21;
int p;
{ double th, lk, link[LLEN];
  int i, ii;
  lk = 0;
  for (i=0; i<des->n; i++)
  { ii = des->ind[i];
    th = (1-a)*innerprod(c0,d_xi(des,i),p) + a*innerprod(c1,d_xi(des,i),p);
    stdlinks(link,&lf->lfd,&lf->sp,ii,th,robscale);
    lk += wght(des,ii)*link[ZLIK];
  }
  des->llk = lk;
  return(vocri(des->llk,t0,(1-a)*t20+a*t21,pen(&lf->sp)));
}

int procvvord(des,lf,v)
design *des;
lfit *lf;
int v;
{ double tr[6], gcv, g0, ap, coef[4][10], t2[4], th, md;
  int i, j, k, d1, i0, p1, ip;
  des->xev = evpt(&lf->fp,v);

  ap = pen(&lf->sp);
  if ((ap==0) & ((fam(&lf->sp)&63)!=TGAUS)) ap = 2.0;
  d1 = deg(&lf->sp); p1 = npar(&lf->sp);
  for (i=0; i<p1; i++) coef[0][i] = coef[1][i] = coef[2][i] = coef[3][i] = 0.0;
  i0 = 0; g0 = 0;
  ip = 1;

  for (i=deg0(&lf->sp); i<=d1; i++)
  { deg(&lf->sp) = i;
    des->p = npar(&lf->sp) = calcp(&lf->sp,lf->lfd.d);
    k = locfit(&lf->lfd,des,&lf->sp,0, i==deg0(&lf->sp),0);

    local_df(&lf->lfd,&lf->sp,des,tr);
    gcv = vocri(des->llk,tr[0],tr[2],ap);
    if ((i==deg0(&lf->sp)) || (gcv<g0)) { i0 = i; g0 = gcv; md = i; }

    for (j=0; j<des->p; j++) coef[i][j] = des->cf[j];
    t2[i] = tr[2];

#ifdef RESEARCH
    if ((ip) && (i>deg0(&lf->sp)))
    { for (j=1; j<10; j++)
      { gcv = intvo(des,lf,coef[i-1],coef[i],j/10.0,des->p,tr[0],t2[i-1],t2[i]);
        if (gcv<g0) { g0 = gcv; md = i-1+j/10.0; }
      }
    }
#endif
  }
  lf->fp.h[v] = des->h;
  if (lf->fp.h[v]<=0) WARN(("zero bandwidth in procvvord"));

  if (i0<d1) /* recompute the best fit */
  { deg(&lf->sp) = i0;
    des->p = npar(&lf->sp) = calcp(&lf->sp,lf->lfd.d);
    k = locfit(&lf->lfd,des,&lf->sp,0,0,0);
    for (i=npar(&lf->sp); i<p1; i++) des->cf[i] = 0.0;
    i0 = md; if (i0==d1) i0--;
    th = md-i0;
    for (i=0; i<p1; i++) des->cf[i] = (1-th)*coef[i0][i]+th*coef[i0+1][i];
    deg(&lf->sp) = d1; npar(&lf->sp) = p1;
  }

  for (i=0; i<p1; i++) lf->fp.coef[i*lf->fp.nvm+v] = des->cf[i];
  lf->fp.deg[v] = md;
  return(k);
}

void initvord(lf)
lfit *lf;
{ initstd(lf);
  PROCV(lf) = procvvord;
  ALLOC(lf) = allocstd;
  PPROC(lf) = NULL;
  KEEPC(lf) = 0;
  NOPC(lf)  = 1;
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
/*
 * functions for computing and subtracting, adding the
 * parametric component
 */

#include "lfev.h"

int noparcomp(sp)
smpar *sp;
{ int tg;
  if (ubas(sp)) return(1);
  tg = fam(sp) & 63;
  if (tg<=THAZ) return(1);
  if (tg==TROBT) return(1);
  if (tg==TCAUC) return(1);
  if (tg==TQUANT) return(1);
  return(0);
}

int pc_reqd(d,p)
int d, p;
{ return(d + 2*p + jac_reqd(p));
}

void pcchk(pc,d,p,lc)
paramcomp *pc;
int d, p, lc;
{ int rw;
  double *z;

  rw = pc_reqd(d,p);
  if (pc->lwk < rw)
  { pc->wk = (double *)calloc(rw,sizeof(double));
    if ( pc->wk == NULL ) {
      printf("Problem allocating memory for pc->wk\n");fflush(stdout);
    }
    pc->lwk= rw;
  }
  z = pc->wk;

  pc->xbar = z; z += d;
  pc->coef = z; z += p;
  pc->f    = z; z += p;

  z = jac_alloc(&pc->xtwx,p,z);
  pc->xtwx.p = p;
}

void compparcomp(des,lfd,sp,pc,nopc)
design *des;
lfdata *lfd;
smpar *sp;
paramcomp *pc;
int nopc;
{ int i, j, k, p;
  double wt, sw;

  if (lf_debug>1) mut_printf(" compparcomp:\n");
  p = des->p;
  pcchk(pc,lfd->d,p,1);

  for (i=0; i<lfd->d; i++) pc->xbar[i] = 0.0;
  sw = 0.0;
  for (i=0; i<lfd->n; i++)
  { 
    wt = prwt(lfd,i);
    sw += wt;
    for (j=0; j<lfd->d; j++)
      pc->xbar[j] += datum(lfd,j,i)*wt;
    des->ind[i] = i;
    wght(des,i) = 1.0;
  }
  for (i=0; i<lfd->d; i++) pc->xbar[i] /= sw;
  if ((nopc) || noparcomp(sp))
  { haspc(pc) = 0;
    return;
  }
  haspc(pc) = 1;
  des->xev = pc->xbar;
  k = locfit(lfd,des,sp,0,0,0);
  if (k != LF_OK) lf_status_msg(k);
  if (lf_error) return;
  switch(k)
  { case LF_NOPT: return;
    case LF_INFA: return;
    case LF_NCON: return;
    case LF_OOB: return;
    case LF_NSLN: return;
    case LF_PF:
      WARN(("compparcomp: perfect fit"));
    case LF_OK:
    case LF_DONE:
      for (i=0; i<p; i++)
      { pc->coef[i] = des->cf[i];
        pc->xtwx.dg[i] = des->xtwx.dg[i];
        pc->xtwx.wk[i] = des->xtwx.wk[i];
      }
      for (i=0; i<p*p; i++)
      { pc->xtwx.Z[i] = des->xtwx.Z[i];
        pc->xtwx.Q[i] = des->xtwx.Q[i];
      }
      pc->xtwx.sm = des->xtwx.sm;
      pc->xtwx.st = des->xtwx.st;
      return;
    default:
      LERR(("compparcomp: locfit unknown return status %d",k));
      return;
  }
}

void subparcomp(des,lf,coef)
design *des;
lfit *lf;
double *coef;
{ int i, nd;
  deriv *dv;
  paramcomp *pc;

  pc = &lf->pc;
  if (!haspc(pc)) return;

  dv = &lf->dv; nd = dv->nd;
  fitfun(&lf->lfd, &lf->sp, des->xev,pc->xbar,des->f1,dv);
  coef[0] -= innerprod(pc->coef,des->f1,pc->xtwx.p);
  if (des->ncoef == 1) return;

  dv->nd = nd+1;
  for (i=0; i<lf->lfd.d; i++)
  { dv->deriv[nd] = i;
    fitfun(&lf->lfd, &lf->sp, des->xev,pc->xbar,des->f1,dv);
    coef[i+1] -= innerprod(pc->coef,des->f1,pc->xtwx.p);
  }
  dv->nd = nd;
}

void subparcomp2(des,lf,vr,il)
design *des;
lfit *lf;
double *vr, *il;
{ double t0, t1;
  int i, nd;
  deriv *dv;
  paramcomp *pc;

  pc = &lf->pc;
  if (!haspc(pc)) return;

  dv = &lf->dv; nd = dv->nd;

  fitfun(&lf->lfd, &lf->sp, des->xev,pc->xbar,des->f1,dv);
  for (i=0; i<npar(&lf->sp); i++) pc->f[i] = des->f1[i];
  jacob_solve(&pc->xtwx,des->f1);
  t0 = sqrt(innerprod(pc->f,des->f1,pc->xtwx.p));
  vr[0] -= t0;
  il[0] -= t0;
  if ((t0==0) | (des->ncoef==1)) return;

  dv->nd = nd+1;
  for (i=0; i<lf->lfd.d; i++)
  { dv->deriv[nd] = i;
    fitfun(&lf->lfd, &lf->sp, des->xev,pc->xbar,pc->f,dv);
    t1 = innerprod(pc->f,des->f1,pc->xtwx.p)/t0;
    vr[i+1] -= t1;
    il[i+1] -= t1;
  }
  dv->nd = nd;
}

double addparcomp(lf,x,c)
lfit *lf;
double *x;
int c;
{ double y;
  paramcomp *pc;

  pc = &lf->pc;
  if (!haspc(pc)) return(0.0);
  fitfun(&lf->lfd, &lf->sp, x,pc->xbar,pc->f,&lf->dv);
  if (c==PCOEF) return(innerprod(pc->coef,pc->f,pc->xtwx.p));
  if ((c==PNLX)|(c==PT0)|(c==PVARI))
  { y = sqrt(jacob_qf(&pc->xtwx,pc->f));
    return(y);
  }
  return(0.0);
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
#include "lfev.h"

/*
  preplot():  interpolates the fit to a new set of points.
  lf  -- the fit structure.
  x   -- the points to predict at.
  f   -- vector to return the predictions.
  se  -- vector to return std errors (NULL if not req'd)
  band-- char for conf band type. ('n'=none, 'g'=global etc.)
  n   -- no of predictions (or vector of margin lengths for grid)
  where -- where to predict:
           1 = points in the array x.
           2 = grid defined by margins in x.
           3 = data points from lf (ignore x).
           4 = fit points from lf (ignore x).
  what -- what to predict.
           (PCOEF etc; see lfcons.h file)

*/

#define NWH 8
static char *whtyp[NWH] = { "coef", "nlx", "infl", "band",
                          "degr", "like", "rdf", "vari" };
static int   whval[NWH] = { PCOEF, PNLX, PT0, PBAND, PDEGR, PLIK, PRDF, PVARI };
int ppwhat(z)
char *z;
{ int val;

  val = pmatch(z, whtyp, whval, NWH, -1);
  if (val==-1) LERR(("Unknown what = %s",z));
  return(val);
}

static char cb;
double *sef, *fit, sigmahat;

void predptall(lf,x,what,ev,i)
lfit *lf;
double *x;
int what, ev, i;
{ double lik, rdf;
  fit[i] = dointpoint(lf,x,what,ev,i);
  if (cb=='n') return;
  sef[i] = dointpoint(lf,x,PNLX,ev,i);
  if (cb=='g')
  { sef[i] *= sigmahat;
    return;
  }
  if (cb=='l')
  { lik = dointpoint(lf,x,PLIK,ev,i);
    rdf = dointpoint(lf,x,PRDF,ev,i);
    sef[i] *= sqrt(-2*lik/rdf);
    return;
  }
  if (cb=='p')
  { sef[i] = sigmahat*sqrt(1+sef[i]*sef[i]);
    return;
  }
}

void predptdir(lf,des,x,what,i)
lfit *lf;
design *des;
double *x;
int what, i;
{ int needcv;
  des->xev = x;
  needcv = (what==PVARI) | (what==PNLX) | (what==PT0) | (what==PRDF);
  locfit(&lf->lfd,des,&lf->sp,0,1,needcv);
  switch(what)
  { case PCOEF: fit[i] = des->cf[0]; break;
    case PVARI: fit[i] = des->V[0]; break;
    case PNLX:  fit[i] = sqrt(des->V[0]); break;
    case PT0:   fit[i] = des->f1[0]; break;
    case PBAND: fit[i] = des->h; break;
    case PDEGR: fit[i] = deg(&lf->sp); break;
    case PLIK:  fit[i] = des->llk; break;
    case PRDF:  fit[i] = des->tr0 - des->tr2; break;
    default: LERR(("unknown what in predptdir"));
  }
}

void prepvector(lf,des,x,n,what,dir) /* interpolate a vector */
lfit *lf;
design *des;
double **x;
int n, what, dir;
{ int i, j;
  double xx[MXDIM];
  for (i=0; i<n; i++)
  { for (j=0; j<lf->fp.d; j++) xx[j] = x[j][i];
    if (dir)
      predptdir(lf,des,xx,what,i);
    else
      predptall(lf,xx,what,ev(&lf->evs),i);
    if (lf_error) return;
  }
}

void prepfitp(lf,what)
lfit *lf;
int what;
{ int  i;
  for (i=0; i<lf->fp.nv; i++)
  { predptall(lf,evpt(&lf->fp,i),what,EFITP,i);
    if (lf_error) return;
  }
}

void prepgrid(lf,des,x,mg,n,what,dir) /* interpolate a grid given margins */
lfit *lf;
design *des;
double **x;
int *mg, dir, n, what;
{ int i, ii, j, d;
  double xv[MXDIM];
  d = lf->fp.d;
  for (i=0; i<n; i++)
  { ii = i;
    for (j=0; j<d; j++)
    { xv[j] = x[j][ii%mg[j]];
      ii /= mg[j];
    }
    if (dir)
      predptdir(lf,des,xv,what,i);
    else
      predptall(lf,xv,what,ev(&lf->evs),i);
    if (lf_error) return;
  }
}

void preplot(lf,x,f,se,band,mg,where,what,dir)
lfit *lf;
double **x, *f, *se;
int *mg, where, what, dir;
char band;
{ int d, i, n;
  double *xx[MXDIM];
  design ppdes;
  d = lf->fp.d;
  fit = f;
  sef = se;
  cb = band;
  if (cb!='n') sigmahat = sqrt(rv(&lf->fp));
  if (dir) des_init(&ppdes,lf->lfd.n,npar(&lf->sp));

  switch(where)
  { case 1: /* vector */
      n = mg[0];
      prepvector(lf,&ppdes,x,n,what,dir);
      break;
    case 2: /* grid */
      n = 1;
      for (i=0; i<d; i++) n *= mg[i];
      prepgrid(lf,&ppdes,x,mg,n,what,dir);
      break;
    case 3: /* data */
      n = lf->lfd.n;
      if ((ev(&lf->evs)==EDATA) | (ev(&lf->evs)==ECROS))
      { prepfitp(lf,what);
        dir = 0;
      }
      else
      { for (i=0; i<d; i++) xx[i] = dvari(&lf->lfd,i);
        prepvector(lf,&ppdes,xx,n,what,dir);
      }
      break;
    case 4: /* fit points */
      n = lf->fp.nv; dir = 0;
      prepfitp(lf,what);
      break;
    default:
      LERR(("unknown where in preplot"));
  }

  if ((!dir) && ((what==PT0)|(what==PVARI)))
    for (i=0; i<n; i++) f[i] = f[i]*f[i];
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
#include "lfev.h"

int procv_nov(des,lf,v)
design *des;
lfit *lf;
int v;
{ int lf_status;

  if (lf_debug>1) mut_printf(" procveraw: %d\n",v);
  des->xev = evpt(&lf->fp,v);

  if (acri(&lf->sp)==ANONE)
    lf_status = locfit(&lf->lfd,des,&lf->sp,0,1,0);
  else
    lf_status = alocfit(&lf->lfd,&lf->sp,&lf->dv,des,0);
  if (lf->fp.h != NULL) lf->fp.h[v] = des->h;

  return(lf_status);
}

int procv_var(des,lf,v)
design *des;
lfit *lf;
int v;
{ int i, lf_status;

  if (lf_debug>1) mut_printf(" procvraw: %d\n",v);
  des->xev = evpt(&lf->fp,v);

  if (acri(&lf->sp)==ANONE)
    lf_status = locfit(&lf->lfd,des,&lf->sp,0,1,1);
  else
    lf_status = alocfit(&lf->lfd,&lf->sp,&lf->dv,des,1);
  if (lf->fp.h != NULL) lf->fp.h[v] = des->h;

  return(lf_status);
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
/*
 * startmodule(lf,des,mod,dir) -- the standard entry point.
 *   des and lf are pointers to the design and fit structures.
 *   mod - module name. Set to NULL if the module is already
 *         initialized.
 *   dir - for dynamic modules, the directory.
 *
 * initmodule(mdl,mod,dir,lf)
 *   direct call for module initialization.
 *
 */

#include "lfev.h"

#ifdef WINDOWS

#define DIRSEP '\\'
#define PATHSEP ';'

#else

#define DIRSEP '/'
#define PATHSEP ':'

#endif


#ifdef ALLOW_MODULES

#ifdef WINDOWS
#include <windows.h>
#define DLEXT "dll"
#define DLOPEN(x) LoadLibrary(x)
#define DLSYM GetProcAddress

#else

#include <dlfcn.h>
#define DLEXT "so"
#define DLOPEN(x) dlopen(x,RTLD_LAZY)
#define DLSYM dlsym
#endif

#endif

static double fpkap[6];
void fitpt_init(fp)
fitpt *fp;
{ 
  dc(fp) = 0;
  geth(fp) = GSTD;
  fp->nv = fp->nvm = 0;
  if (fp->kap==NULL) fp->kap = fpkap;
}

void lfit_init(lf)
lfit *lf;
{
  lfdata_init(&lf->lfd);
  evstruc_init(&lf->evs);
  smpar_init(&lf->sp,&lf->lfd);
  deriv_init(&lf->dv);
  fitpt_init(&lf->fp);
  lf->mdl.np = 0;
}

void fitdefault(lf)
lfit *lf;
{ WARN(("fitdefault deprecated -- use lfit_init()"));
  lfit_init(lf);
}

void set_flim(lfd,evs)
lfdata *lfd;
evstruc *evs;
{ int i, j, d, n;
  double z, mx, mn, *bx;

  if (ev(evs)==ESPHR) return;
  d = lfd->d; n = lfd->n;
  bx = evs->fl;
  for (i=0; i<d; i++)
    if (bx[i]==bx[i+d])
    { if (lfd->sty[i]==STANGL)
      { bx[i] = 0.0; bx[i+d] = 2*PI*lfd->sca[i];
      }
      else
      { mx = mn = datum(lfd,i,0);
        for (j=1; j<n; j++)
        { mx = MAX(mx,datum(lfd,i,j));
          mn = MIN(mn,datum(lfd,i,j));
        }
        if (lfd->xl[i]<lfd->xl[i+d]) /* user set xlim; maybe use them. */
        { z = mx-mn;
          if (mn-0.2*z < lfd->xl[i]) mn = lfd->xl[i];
          if (mx+0.2*z > lfd->xl[i+d]) mx = lfd->xl[i+d];
        }
        bx[i] = mn;
        bx[i+d] = mx;
      }
    }
}

double vvari(v,n)
double *v;
int n;
{ int i;
  double xb, s2;
  xb = s2 = 0.0;
  for (i=0; i<n; i++) xb += v[i];
  xb /= n;
  for (i=0; i<n; i++) s2 += SQR(v[i]-xb);
  return(s2/(n-1));
}

void set_scales(lfd)
lfdata *lfd;
{ int i;
  for (i=0; i<lfd->d; i++)
    if (lfd->sca[i]<=0) /* set automatic scales */
    { if (lfd->sty[i]==STANGL)
        lfd->sca[i] = 1.0;
      else lfd->sca[i] = sqrt(vvari(lfd->x[i],lfd->n));
    }
}

void nstartlf(des,lf)
design *des;
lfit *lf;
{ int i, d, n;

  if (lf_debug>0) mut_printf("nstartlf\n");
  n = lf->lfd.n;
  d = lf->lfd.d;
  npar(&lf->sp) = calcp(&lf->sp,d);

  des_init(des,n,npar(&lf->sp));
  set_scales(&lf->lfd);
  set_flim(&lf->lfd,&lf->evs);
  compparcomp(des,&lf->lfd,&lf->sp,&lf->pc,lf->mdl.nopc);
  if (lf_error) return;
  makecfn(&lf->sp,des,&lf->dv,lf->lfd.d);

  lf->lfd.ord = 0;
  if ((d==1) && (lf->lfd.sty[0]!=STANGL))
  { i = 1;
    while ((i<n) && (datum(&lf->lfd,0,i)>=datum(&lf->lfd,0,i-1))) i++;
    lf->lfd.ord = (i==n);
  }
  for (i=0; i<npar(&lf->sp); i++) des->fix[i] = 0;

  lf->fp.d = lf->lfd.d;
  lf->fp.hasd = (des->ncoef==(1+lf->fp.d));
  lf->fp.nv = lf->evs.nce = 0;

  if (lf_debug>1) mut_printf("call eval structure %d\n",ev(&lf->evs));
  switch(ev(&lf->evs))
  { case EPHULL: triang_start(des,lf); break;
    case EDATA:  dataf(des,lf); break;
    case ECROS:  crossf(des,lf); break;
    case EGRID:  gridf(des,lf); break;
    case ETREE:  atree_start(des,lf); break;
    case EKDCE:  kt(&lf->sp) = KCE;
    case EKDTR:  kdtre_start(des,lf); break;
    case EPRES:  preset(des,lf); break;
    case EXBAR:  xbarf(des,lf); break;
    case ENONE:  return;
    case ESPHR:  sphere_start(des,lf); break;
    case ESPEC:  lf->evs.espec(des,lf); break;
    default: LERR(("startlf: Invalid evaluation structure %d",ev(&lf->evs)));
  }

  /* renormalize for family=density */
  if ((de_renorm) && (fam(&lf->sp)==TDEN)) dens_renorm(lf,des);
}

/*
 * getnextdir() gets the next dir from a string dirpath="dir1:dir2:dir3:..."
 *   (;-separated on windows).
 *   The directory is returned through dirnext, and the function returns
 *   a pointer to the next string.
 *   typical usage is recursive, dirpath = getnextdir(dirpath,dirnext).
 *   with the above example, this sets dirnext="dir1" and dirpath="dir2:dir3:...".
 * if the input dirpath has no :, then it is copied to dirnext, and return is "".
 * if input dirpath is "", dirnext is set to "", and null pointer returned.
 */
char *getnextdir(dirpath,dirnext)
char *dirpath, *dirnext;
{ char *z;
  if (strlen(dirpath)==0)
  { sprintf(dirnext,"");
    return(NULL);
  }

  z = strchr(dirpath,PATHSEP);
  if (z==NULL)
  { sprintf(dirnext,"%s%c",dirpath,DIRSEP);
    return(&dirpath[strlen(dirnext)]);
  }

  *z = '\0';
  sprintf(dirnext,"%s%c",dirpath,DIRSEP);
  return(++z);
}

int initmodule(mdl, mod, dir, lf)
module *mdl;
lfit *lf;
char *mod, *dir;
{ int n, d, p;
#ifdef ALLOW_MODULES
#ifdef WINDOWS
HINSTANCE res;
typedef void (CALLBACK* DLLFN)();
DLLFN init;
#else
void *res;
void (*init)();
#endif
  char distname[500];
#endif

  n = lf->lfd.n;
  d = lf->lfd.d;
  p = npar(&lf->sp) = calcp(&lf->sp,d);

  mdl->isset = 1;
  PPROC(lf) = NULL;
  if (strcmp(mod,"std")==0)    { initstd(lf); return(1); }
  if (strcmp(mod,"simple")==0) { initsimple(lf); return(1); }
  if (strcmp(mod,"allcf")==0)  { initallcf(lf); return(1); }
  if (strcmp(mod,"hatm")==0)   { inithatm(lf); return(1); }
  if (strcmp(mod,"kappa")==0)  { initkappa(lf); return(1); }
  if (strcmp(mod,"lscv")==0)   { initlscv(lf); return(1); }
  if (strcmp(mod,"gamf")==0)   { initgam(lf); return(1); }
  if (strcmp(mod,"gamp")==0)   { initgam(lf); return(1); }
  if (strcmp(mod,"rband")==0)  { initrband(lf); return(1); }
  if (strcmp(mod,"scb")==0)    { initscb(lf); return(1); }
  if (strcmp(mod,"vord")==0)   { initvord(lf); return(1); }

#ifdef ALLOW_MODULES
  while (dir != NULL)
  {
    dir = getnextdir(dir,distname);
    sprintf(&distname[strlen(distname)],"mod%s.%s",mod,DLEXT);
    res = DLOPEN(distname);
    if (res==NULL)
    {
#ifdef WINDOWS
      mut_printf("LoadLibrary failed: %s, %d\n",distname,GetLastError());
#else
      mut_printf("dlopen failed: %s\n",dlerror());
#endif
    }
    else
    {
#ifdef WINDOWS
      mut_printf("LoadLibrary success: %s\n",distname);
#else
      mut_printf("dlopen success: %s\n",distname);
#endif
      sprintf(distname,"init%s",mod);
      init = (void *)DLSYM(res,distname);
      if (init==NULL)
      { mut_printf("I can't find %s() function.\n",distname);
        mdl->isset = 0;
        return(0);
      }
      init(lf);
      return(1);
    }
  }
#endif

  mdl->isset = 0;
  return(0);
}

/*
 * startmodule is the entry point to launch the fit.
 * if mod is provided, will first initialize the module.
 * if mod==NULL, assumes the module has been initialized separately.
 */
void startmodule(lf,des,mod,dir)
lfit *lf;
design *des;
char *mod, *dir;
{ int z;

  if (mod != NULL)
  { z = initmodule(&lf->mdl,mod,dir,lf);
    if (!z) return;
  }

  lf->fp.nv = lf->evs.nce = 0;
  if (lf_error) return;
  if (PROCV(lf) != NULL) nstartlf(des,lf);
  if (lf_error) return;
  if (PPROC(lf) != NULL) PPROC(lf)(lf,des,lf->fp.kap);
}

/* for compatability, more or less. */
void startlf(des,lf,vfun,nopc)
design *des;
lfit *lf;
int (*vfun)(), nopc;
{ int z;
  z = initmodule(&lf->mdl,"std",NULL,lf);
  if (!z) return;
  lf->mdl.procv = vfun;
  lf->mdl.nopc = nopc;
  nstartlf(des,lf);
}
