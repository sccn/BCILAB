/*
 * Copyright 1996-2006 Catherine Loader.
 */

#include "mex.h"
/*
 * Copyright 1996-2006 Catherine Loader.
 */
/*
 *   Integration for hazard rate estimation. The functions in this
 *   file are used to evaluate
 *      sum int_0^{Ti} W_i(t,x) A()A()' exp( P() ) dt
 *   for hazard rate models.
 *
 *   These routines assume the weight function is supported on [-1,1].
 *   hasint_sph multiplies by exp(base(lf,i)), which allows estimating
 *   the baseline in a proportional hazards model, when the covariate
 *   effect base(lf,i) is known.
 *
 *   TODO:
 *     hazint_sph, should be able to reduce mint in some cases with
 *       small integration range. onedint could be used for beta-family
 *       (RECT,EPAN,BISQ,TRWT) kernels.
 *     hazint_prod, restrict terms from the sum based on x values.
 *       I should count obs >= max, and only do that integration once.
 */

#include "locf.h"

static double ilim[2*MXDIM], *ff, tmax;
static lfdata *haz_lfd;
static smpar  *haz_sp;

/*
 *  hrao returns 0 if integration region is empty.
 *               1 otherwise.
 */
int haz_sph_int(dfx,cf,h,r1)
double *dfx, *cf, h, *r1;
{ double s, t0, t1, wt, th;
  int j, dim, p;
  s = 0; p = npar(haz_sp);
  dim = haz_lfd->d;
  for (j=1; j<dim; j++) s += SQR(dfx[j]/(h*haz_lfd->sca[j]));
  if (s>1) return(0);

  setzero(r1,p*p);
  t1 = sqrt(1-s)*h*haz_lfd->sca[0];
  t0 = -t1;
  if (t0<ilim[0])   t0 = ilim[0];
  if (t1>ilim[dim]) t1 = ilim[dim];
  if (t1>dfx[0]) t1 = dfx[0];
  if (t1<t0) return(0);

/*  Numerical integration by Simpson's rule.
 */
  for (j=0; j<=de_mint; j++)
  { dfx[0] = t0+(t1-t0)*j/de_mint;
    wt = weight(haz_lfd, haz_sp, dfx, NULL, h, 0, 0.0);
    fitfun(haz_lfd, haz_sp, dfx,NULL,ff,NULL);
    th = innerprod(cf,ff,p);
    if (link(haz_sp)==LLOG) th = exp(th);
    wt *= 2+2*(j&1)-(j==0)-(j==de_mint);
    addouter(r1,ff,ff,p,wt*th);
  }
  multmatscal(r1,(t1-t0)/(3*de_mint),p*p);

  return(1);
}

int hazint_sph(t,resp,r1,cf,h)
double *t, *resp, *r1, *cf, h;
{ int i, j, n, p, st;
  double dfx[MXDIM], eb, sb;
  p = npar(haz_sp);
  setzero(resp,p*p);
  sb = 0.0;

  n = haz_lfd->n;
  for (i=0; i<=n; i++)
  {
    if (i==n)
    { dfx[0] = tmax-t[0];
      for (j=1; j<haz_lfd->d; j++) dfx[j] = 0.0;
      eb = exp(sb/n);
    }
    else
    { eb = exp(base(haz_lfd,i)); sb += base(haz_lfd,i);
      for (j=0; j<haz_lfd->d; j++) dfx[j] = datum(haz_lfd,j,i)-t[j];
    }

    st = haz_sph_int(dfx,cf,h,r1);
    if (st)
      for (j=0; j<p*p; j++) resp[j] += eb*r1[j];
  }
  return(LF_OK);
}

int hazint_prod(t,resp,x,cf,h)
double *t, *resp, *x, *cf, h;
{ int d, p, i, j, k, st;
  double dfx[MXDIM], t_prev,
         hj, hs, ncf[MXDEG], ef, il1;
  double prod_wk[MXDIM][2*MXDEG+1], eb, sb;

  p = npar(haz_sp);
  d = haz_lfd->d;
  setzero(resp,p*p);
  hj = hs = h*haz_lfd->sca[0];

  ncf[0] = cf[0];
  for (i=1; i<=deg(haz_sp); i++)
  { ncf[i] = hj*cf[(i-1)*d+1]; hj *= hs;
  }

/*   for i=0..n....
 *     First we compute prod_wk[j], j=0..d.
 *     For j=0, this is int_0^T_i (u-t)^k W((u-t)/h) exp(b0*(u-t)) du
 *     For remaining j,   (x(i,j)-x(j))^k Wj exp(bj*(x..-x.))
 *
 *     Second, we add to the integration (exp(a) incl. in integral)
 *     with the right factorial denominators.
 */
  t_prev = ilim[0]; sb = 0.0;
  for (i=0; i<=haz_lfd->n; i++)
  { if (i==haz_lfd->n)
    { dfx[0] = tmax-t[0];
      for (j=1; j<d; j++) dfx[j] = 0.0;
      eb = exp(sb/haz_lfd->n);
    }
    else
    { eb = exp(base(haz_lfd,i)); sb += base(haz_lfd,i);
      for (j=0; j<d; j++) dfx[j] = datum(haz_lfd,j,i)-t[j];
    }

    if (dfx[0]>ilim[0]) /* else it doesn't contribute */
    {
/* time integral */
      il1 = (dfx[0]>ilim[d]) ? ilim[d] : dfx[0];
      if (il1 != t_prev) /* don't repeat! */
      { st = onedint(haz_sp,ncf,ilim[0]/hs,il1/hs,prod_wk[0]);
        if (st>0) return(st);
        hj = eb;
        for (j=0; j<=2*deg(haz_sp); j++)
        { hj *= hs;
          prod_wk[0][j] *= hj;
        }
        t_prev = il1;
      }

/* covariate terms */
      for (j=1; j<d; j++)
      {
        ef = 0.0;
        for (k=deg(haz_sp); k>0; k--) ef = (ef+dfx[j])*cf[1+(k-1)*d+j];
        ef = exp(ef);
        prod_wk[j][0] = ef * W(dfx[j]/(h*haz_lfd->sca[j]),ker(haz_sp));
        for (k=1; k<=2*deg(haz_sp); k++)
          prod_wk[j][k] = prod_wk[j][k-1] * dfx[j];
      }

/*  add to the integration.  */
      prodintresp(resp,prod_wk,d,deg(haz_sp),p);
    } /* if dfx0 > ilim0 */
  } /* n loop */

/* symmetrize */
  for (k=0; k<p; k++)
    for (j=k; j<p; j++)
      resp[j*p+k] = resp[k*p+j];
  return(LF_OK);
}

int hazint(t,resp,resp1,cf,h)
double *t, *resp, *resp1, *cf, h;
{ if (haz_lfd->d==1) return(hazint_prod(t,resp,resp1,cf,h));
  if (kt(haz_sp)==KPROD) return(hazint_prod(t,resp,resp1,cf,h));

  return(hazint_sph(t,resp,resp1,cf,h));
}

void haz_init(lfd,des,sp,il)
lfdata *lfd;
design *des;
smpar *sp;
double *il;
{ int i;
  
  haz_lfd = lfd;
  haz_sp  = sp;

  tmax = datum(lfd,0,0);
  for (i=1; i<lfd->n; i++) tmax = MAX(tmax,datum(lfd,0,i));
  ff = des->xtwx.wk;
  for (i=0; i<2*lfd->d; i++) ilim[i] = il[i];
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
/*
 *
 *  Routines for one-dimensional numerical integration
 *  in density estimation. The entry point is
 *
 *  onedint(cf,mi,l0,l1,resp)
 *
 *  which evaluates int W(u)u^j exp( P(u) ), j=0..2*deg.
 *  P(u) = cf[0] + cf[1]u + cf[2]u^2/2 + ... + cf[deg]u^deg/deg!
 *  l0 and l1 are the integration limits.
 *  The results are returned through the vector resp.
 *
 */

#include "locf.h"

static int debug;

int exbctay(b,c,n,z) /* n-term taylor series of e^(bx+cx^2) */
double b, c, *z;
int n;
{ double ec[20];
  int i, j;
  z[0] = 1;
  for (i=1; i<=n; i++) z[i] = z[i-1]*b/i;
  if (c==0.0) return(n);
  if (n>=40)
  { WARN(("exbctay limit to n<40"));
    n = 39;
  }
  ec[0] = 1;
  for (i=1; 2*i<=n; i++) ec[i] = ec[i-1]*c/i;
  for (i=n; i>1; i--)
    for (j=1; 2*j<=i; j++)
      z[i] += ec[j]*z[i-2*j];
  return(n);
}

double explinjtay(l0,l1,j,cf)
/* int_l0^l1 x^j e^(a+bx+cx^2); exbctay aroud l1 */
double l0, l1, *cf;
int j;
{ double tc[40], f, s;
  int k, n;
  if ((l0!=0.0) | (l1!=1.0)) WARN(("explinjtay: invalid l0, l1"));
  n = exbctay(cf[1]+2*cf[2]*l1,cf[2],20,tc);
  s = tc[0]/(j+1);
  f = 1/(j+1);
  for (k=1; k<=n; k++)
  { f *= -k/(j+k+1.0);
    s += tc[k]*f;
  }
  return(f);
}

void explint1(l0,l1,cf,I,p) /* int x^j exp(a+bx); j=0..p-1 */
double l0, l1, *cf, *I;
int p;
{ double y0, y1, f;
  int j, k, k1;
  y0 = mut_exp(cf[0]+l0*cf[1]);
  y1 = mut_exp(cf[0]+l1*cf[1]);
  if (p<2*fabs(cf[1])) k = p; else k = (int)fabs(cf[1]);

  if (k>0)
  { I[0] = (y1-y0)/cf[1];
    for (j=1; j<k; j++) /* forward steps for small j */
    { y1 *= l1; y0 *= l0;
      I[j] = (y1-y0-j*I[j-1])/cf[1];
    }
    if (k==p) return;
    y1 *= l1; y0 *= l0;
  }

  f = 1; k1 = k;
  while ((k<50) && (f>1.0e-8)) /* initially Ik = diff(x^{k+1}e^{a+bx}) */
  { y1 *= l1; y0 *= l0;
    I[k] = y1-y0;
    if (k>=p) f *= fabs(cf[1])/(k+1);
    k++;
  }
  if (k==50) WARN(("explint1: want k>50"));
  I[k] = 0.0;
  for (j=k-1; j>=k1; j--) /* now do back step recursion */
    I[j] = (I[j]-cf[1]*I[j+1])/(j+1);
}

void explintyl(l0,l1,cf,I,p) /* small c, use taylor series and explint1 */
double l0, l1, *cf, *I;
int p;
{ int i;
  double c;
  explint1(l0,l1,cf,I,p+8);
  c = cf[2];
  for (i=0; i<p; i++)
    I[i] = (((I[i+8]*c/4+I[i+6])*c/3+I[i+4])*c/2+I[i+2])*c+I[i];
}

void solvetrid(X,y,m)
double *X, *y;
int m;
{ int i;
  double s;
  for (i=1; i<m; i++)
  { s = X[3*i]/X[3*i-2];
    X[3*i] = 0; X[3*i+1] -= s*X[3*i-1];
    y[i] -= s*y[i-1];
  }
  for (i=m-2; i>=0; i--)
  { s = X[3*i+2]/X[3*i+4];
    X[3*i+2] = 0;
    y[i] -= s*y[i+1];
  }
  for (i=0; i<m; i++) y[i] /= X[3*i+1];
}

void initi0i1(I,cf,y0,y1,l0,l1)
double *I, *cf, y0, y1, l0, l1;
{ double a0, a1, c, d, bi;
  d = -cf[1]/(2*cf[2]); c = sqrt(2*fabs(cf[2]));
  a0 = c*(l0-d); a1 = c*(l1-d);
  if (cf[2]<0)
  { bi = mut_exp(cf[0]+cf[1]*d+cf[2]*d*d)/c;
    if (a0>0)
    { if (a0>6) I[0] = (y0*ptail(-a0)-y1*ptail(-a1))/c;
      else I[0] = S2PI*(mut_pnorm(-a0)-mut_pnorm(-a1))*bi;
    }
    else
    { if (a1< -6) I[0] = (y1*ptail(a1)-y0*ptail(a0))/c;
      else I[0] = S2PI*(mut_pnorm(a1)-mut_pnorm(a0))*bi;
    }
  }
  else
    I[0] = (y1*mut_daws(a1)-y0*mut_daws(a0))/c;
  I[1] = (y1-y0)/(2*cf[2])+d*I[0];
}

void explinsid(l0,l1,cf,I,p) /* large b; don't use fwd recursion */
double l0, l1, *cf, *I;
int p;
{ int k, k0, k1, k2;
  double y0, y1, Z[150];
if (debug) mut_printf("side: %8.5f %8.5f %8.5f    limt %8.5f %8.5f  p %2d\n",cf[0],cf[1],cf[2],l0,l1,p);
 
  k0 = 2;
  k1 = (int)(fabs(cf[1])+fabs(2*cf[2]));
  if (k1<2) k1 = 2;
  if (k1>p+20) k1 = p+20;
  k2 = p+20;

if (k2>50) { mut_printf("onedint: k2 warning\n"); k2 = 50; }
  if (debug) mut_printf("k0 %2d  k1 %2d  k2 %2d  p %2d\n",k0,k1,k2,p);

  y0 = mut_exp(cf[0]+l0*(cf[1]+l0*cf[2]));
  y1 = mut_exp(cf[0]+l1*(cf[1]+l1*cf[2]));
  initi0i1(I,cf,y0,y1,l0,l1);
if (debug) mut_printf("i0 %8.5f  i1 %8.5f\n",I[0],I[1]);

  y1 *= l1; y0 *= l0; /* should be x^(k1)*exp(..) */
  if (k0<k1) /* center steps; initially x^k*exp(...) */
    for (k=k0; k<k1; k++)
    { y1 *= l1; y0 *= l0;
      I[k] = y1-y0;
      Z[3*k] = k; Z[3*k+1] = cf[1]; Z[3*k+2] = 2*cf[2];
    }
   
  y1 *= l1; y0 *= l0; /* should be x^(k1)*exp(..) */
if (debug) mut_printf("k1 %2d  y0 %8.5f  y1 %8.5f\n",k1,y0,y1);
  for (k=k1; k<k2; k++)
  { y1 *= l1; y0 *= l0;
    I[k] = y1-y0;
  }
  I[k2] = I[k2+1] = 0.0;
  for (k=k2-1; k>=k1; k--)
    I[k] = (I[k]-cf[1]*I[k+1]-2*cf[2]*I[k+2])/(k+1);

  if (k0<k1)
  { I[k0] -= k0*I[k0-1];
    I[k1-1] -= 2*cf[2]*I[k1];
    Z[3*k0] = Z[3*k1-1] = 0;
    solvetrid(&Z[3*k0],&I[k0],k1-k0);
  }
if (debug)
{ mut_printf("explinsid:\n");
  for (k=0; k<p; k++) mut_printf("  %8.5f\n",I[k]);
}
}

void explinbkr(l0,l1,cf,I,p) /* small b,c; use back recursion */
double l0, l1, *cf, *I;
int p;
{ int k, km;
  double y0, y1;
  y0 = mut_exp(cf[0]+l0*(cf[1]+cf[2]*l0));
  y1 = mut_exp(cf[0]+l1*(cf[1]+cf[2]*l1));
  km = p+10;
  for (k=0; k<=km; k++)
  { y1 *= l1; y0 *= l0;
    I[k] = y1-y0;
  }
  I[km+1] = I[km+2] = 0;
  for (k=km; k>=0; k--)
    I[k] = (I[k]-cf[1]*I[k+1]-2*cf[2]*I[k+2])/(k+1);
}

void explinfbk0(l0,l1,cf,I,p) /* fwd and bac recur; b=0; c<0 */
double l0, l1, *cf, *I;
int p;
{ double y0, y1, f1, f2, f, ml2;
  int k, ks;

  y0 = mut_exp(cf[0]+l0*l0*cf[2]);
  y1 = mut_exp(cf[0]+l1*l1*cf[2]);
  initi0i1(I,cf,y0,y1,l0,l1);

  ml2 = MAX(l0*l0,l1*l1);
  ks = 1+(int)(2*fabs(cf[2])*ml2);
  if (ks<2) ks = 2;
  if (ks>p-3) ks = p;

  /* forward recursion for k < ks */
  for (k=2; k<ks; k++)
  { y1 *= l1; y0 *= l0;
    I[k] = (y1-y0-(k-1)*I[k-2])/(2*cf[2]);
  }
  if (ks==p) return;

  y1 *= l1*l1; y0 *= l0*l0;
  for (k=ks; k<p; k++) /* set I[k] = x^{k+1}e^(a+cx^2) | {l0,l1} */
  { y1 *= l1; y0 *= l0;
    I[k] = y1-y0;
  }

  /* initialize I[p-2] and I[p-1] */
  f1 = 1.0/p; f2 = 1.0/(p-1);
  I[p-1] *= f1; I[p-2] *= f2;
  k = p; f = 1.0;
  while (f>1.0e-8)
  { y1 *= l1; y0 *= l0;
    if ((k-p)%2==0) /* add to I[p-2] */
    { f2 *= -2*cf[2]/(k+1);
      I[p-2] += (y1-y0)*f2;
    }
    else /* add to I[p-1] */
    { f1 *= -2*cf[2]/(k+1);
      I[p-1] += (y1-y0)*f1;
      f *= 2*fabs(cf[2])*ml2/(k+1);
    }
    k++;
  }
  
  /* use back recursion for I[ks..(p-3)] */
  for (k=p-3; k>=ks; k--)
    I[k] = (I[k]-2*cf[2]*I[k+2])/(k+1);
}

void explinfbk(l0,l1,cf,I,p) /* fwd and bac recur; b not too large */
double l0, l1, *cf, *I;
int p;
{ double y0, y1;
  int k, ks, km;

  y0 = mut_exp(cf[0]+l0*(cf[1]+l0*cf[2]));
  y1 = mut_exp(cf[0]+l1*(cf[1]+l1*cf[2]));
  initi0i1(I,cf,y0,y1,l0,l1);

  ks = (int)(3*fabs(cf[2]));
  if (ks<3) ks = 3;
  if (ks>0.75*p) ks = p; /* stretch the forward recurs as far as poss. */
  /* forward recursion for k < ks */
  for (k=2; k<ks; k++)
  { y1 *= l1; y0 *= l0;
    I[k] = (y1-y0-cf[1]*I[k-1]-(k-1)*I[k-2])/(2*cf[2]);
  }
  if (ks==p) return;

  km = p+15;
  y1 *= l1*l1; y0 *= l0*l0;
  for (k=ks; k<=km; k++)
  { y1 *= l1; y0 *= l0;
    I[k] = y1-y0;
  }
  I[km+1] = I[km+2] = 0.0;
  for (k=km; k>=ks; k--)
    I[k] = (I[k]-cf[1]*I[k+1]-2*cf[2]*I[k+2])/(k+1);
}

void recent(I,resp,wt,p,s,x)
double *I, *resp, *wt, x;
int p, s;
{ int i, j;

  /* first, use W taylor series I -> resp */
  for (i=0; i<=p; i++)
  { resp[i] = 0.0;
    for (j=0; j<s; j++) resp[i] += wt[j]*I[i+j];
  }

  /* now, recenter x -> 0 */
  if (x==0) return;
  for (j=0; j<=p; j++) for (i=p; i>j; i--) resp[i] += x*resp[i-1];
}

void recurint(l0,l2,cf,resp,p,ker)
double l0, l2, *cf, *resp;
int p, ker;
{ int i, s;
  double l1, d0, d1, d2, dl, z0, z1, z2, wt[20], ncf[3], I[50], r1[5], r2[5];
if (debug) mut_printf("\nrecurint: %8.5f %8.5f %8.5f   %8.5f %8.5f\n",cf[0],cf[1],cf[2],l0,l2);

  if (cf[2]==0) /* go straight to explint1 */
  { s = wtaylor(wt,0.0,ker);
if (debug) mut_printf("case 1\n");
    explint1(l0,l2,cf,I,p+s);
    recent(I,resp,wt,p,s,0.0);
    return;
  }

  dl = l2-l0;
  d0 = cf[1]+2*l0*cf[2];
  d2 = cf[1]+2*l2*cf[2];
  z0 = cf[0]+l0*(cf[1]+l0*cf[2]);
  z2 = cf[0]+l2*(cf[1]+l2*cf[2]);

  if ((fabs(cf[1]*dl)<1) && (fabs(cf[2]*dl*dl)<1))
  { ncf[0] = z0; ncf[1] = d0; ncf[2] = cf[2];
if (debug) mut_printf("case 2\n");
    s = wtaylor(wt,l0,ker);
    explinbkr(0.0,dl,ncf,I,p+s);
    recent(I,resp,wt,p,s,l0);
    return;
  }

  if (fabs(cf[2]*dl*dl)<0.001) /* small c, use explint1+tay.ser */
  { ncf[0] = z0; ncf[1] = d0; ncf[2] = cf[2];
if (debug) mut_printf("case small c\n");
    s = wtaylor(wt,l0,ker);
    explintyl(0.0,l2-l0,ncf,I,p+s);
    recent(I,resp,wt,p,s,l0);
    return;
  }

  if (d0*d2<=0) /* max/min in [l0,l2] */
  { l1 = -cf[1]/(2*cf[2]);
    z1 = cf[0]+l1*(cf[1]+l1*cf[2]);
    d1 = 0.0;
    if (cf[2]<0) /* peak, integrate around l1 */
    { s = wtaylor(wt,l1,ker);
      ncf[0] = z1; ncf[1] = 0.0; ncf[2] = cf[2];
if (debug) mut_printf("case peak  p %2d  s %2d\n",p,s);
      explinfbk0(l0-l1,l2-l1,ncf,I,p+s);
      recent(I,resp,wt,p,s,l1);
      return;
    }
  }

  if ((d0-2*cf[2]*dl)*(d2+2*cf[2]*dl)<0) /* max/min is close to [l0,l2] */
  { l1 = -cf[1]/(2*cf[2]);
    z1 = cf[0]+l1*(cf[1]+l1*cf[2]);
    if (l1<l0) { l1 = l0; z1 = z0; }
    if (l1>l2) { l1 = l2; z1 = z2; }

    if ((z1>=z0) & (z1>=z2)) /* peak; integrate around l1 */
    { s = wtaylor(wt,l1,ker);
if (debug) mut_printf("case 4\n");
      d1 = cf[1]+2*l1*cf[2];
      ncf[0] = z1; ncf[1] = d1; ncf[2] = cf[2];
      explinfbk(l0-l1,l2-l1,ncf,I,p+s);
      recent(I,resp,wt,p,s,l1);
      return;
    }

    /* trough; integrate [l0,l1] and [l1,l2] */
    for (i=0; i<=p; i++) r1[i] = r2[i] = 0.0;
    if (l0<l1)
    { s = wtaylor(wt,l0,ker);
if (debug) mut_printf("case 5\n");
      ncf[0] = z0; ncf[1] = d0; ncf[2] = cf[2];
      explinfbk(0.0,l1-l0,ncf,I,p+s);
      recent(I,r1,wt,p,s,l0);
    }
    if (l1<l2)
    { s = wtaylor(wt,l2,ker);
if (debug) mut_printf("case 6\n");
      ncf[0] = z2; ncf[1] = d2; ncf[2] = cf[2];
      explinfbk(l1-l2,0.0,ncf,I,p+s);
      recent(I,r2,wt,p,s,l2);
    }
    for (i=0; i<=p; i++) resp[i] = r1[i]+r2[i];
    return;
  }

  /* Now, quadratic is monotone on [l0,l2]; big b; moderate c */
  if (z2>z0+3) /* steep increase, expand around l2 */
  { s = wtaylor(wt,l2,ker);
if (debug) mut_printf("case 7\n");


    ncf[0] = z2; ncf[1] = d2; ncf[2] = cf[2];
    explinsid(l0-l2,0.0,ncf,I,p+s);
    recent(I,resp,wt,p,s,l2);
if (debug) mut_printf("7 resp: %8.5f %8.5f %8.5f %8.5f\n",resp[0],resp[1],resp[2],resp[3]);
    return;
  }

  /* bias towards expansion around l0, because it's often 0 */
if (debug) mut_printf("case 8\n");
  s = wtaylor(wt,l0,ker);
  ncf[0] = z0; ncf[1] = d0; ncf[2] = cf[2];
  explinsid(0.0,l2-l0,ncf,I,p+s);
  recent(I,resp,wt,p,s,l0);
  return;
}

int onedexpl(cf,deg,resp)
double *cf, *resp;
int deg;
{ int i;
  double f0, fr, fl;
  if (deg>=2) LERR(("onedexpl only valid for deg=0,1"));
  if (fabs(cf[1])>=EFACT) return(LF_BADP);

  f0 = exp(cf[0]); fl = fr = 1.0;
  for (i=0; i<=2*deg; i++)
  { f0 *= i+1;
    fl /=-(EFACT+cf[1]);
    fr /=  EFACT-cf[1];
    resp[i] = f0*(fr-fl);
  }
  return(LF_OK);
}

int onedgaus(cf,deg,resp)
double *cf, *resp;
int deg;
{ int i;
  double f0, mu, s2;
  if (deg==3)
  { LERR(("onedgaus only valid for deg=0,1,2"));
    return(LF_ERR);
  }
  if (2*cf[2]>=GFACT*GFACT) return(LF_BADP);

  s2 = 1/(GFACT*GFACT-2*cf[2]);
  mu = cf[1]*s2;
  resp[0] = 1.0;
  if (deg>=1)
  { resp[1] = mu;
    resp[2] = s2+mu*mu;
    if (deg==2)
    { resp[3] = mu*(3*s2+mu*mu);
      resp[4] = 3*s2*s2 + mu*mu*(6*s2+mu*mu);
    }
  }
  f0 = S2PI * exp(cf[0]+mu*mu/(2*s2))*sqrt(s2);
  for (i=0; i<=2*deg; i++) resp[i] *= f0;
  return(LF_OK);
}

int onedint(sp,cf,l0,l1,resp) /* int W(u)u^j exp(..), j=0..2*deg */
smpar *sp;
double *cf, l0, l1, *resp;
{ double u, uj, y, ncf[4], rr[5];
  int i, j;

if (debug) mut_printf("onedint: %f %f %f   %f %f\n",cf[0],cf[1],cf[2],l0,l1);

  if (deg(sp)<=2)
  { for (i=0; i<3; i++) ncf[i] = (i>deg(sp)) ? 0.0 : cf[i];
    ncf[2] /= 2;

    if (ker(sp)==WEXPL) return(onedexpl(ncf,deg(sp),resp));
    if (ker(sp)==WGAUS) return(onedgaus(ncf,deg(sp),resp));

    if (l1>0)
      recurint(MAX(l0,0.0),l1,ncf,resp,2*deg(sp),ker(sp));
    else for (i=0; i<=2*deg(sp); i++) resp[i] = 0;

    if (l0<0)
    { ncf[1] = -ncf[1];
      l0 = -l0; l1 = -l1;
      recurint(MAX(l1,0.0),l0,ncf,rr,2*deg(sp),ker(sp));
    }
    else for (i=0; i<=2*deg(sp); i++) rr[i] = 0.0;

    for (i=0; i<=2*deg(sp); i++)
      resp[i] += (i%2==0) ? rr[i] : -rr[i];

    return(LF_OK);
  }

  /* For degree >= 3, we use Simpson's rule. */
  for (j=0; j<=2*deg(sp); j++) resp[j] = 0.0;
  for (i=0; i<=de_mint; i++)
  { u = l0+(l1-l0)*i/de_mint;
    y = cf[0]; uj = 1;
    for (j=1; j<=deg(sp); j++)
    { uj *= u;
      y += cf[j]*uj/fact[j];
    }
    y = (4-2*(i%2==0)-(i==0)-(i==de_mint)) *
          W(fabs(u),ker(sp))*exp(MIN(y,300.0));
    for (j=0; j<=2*deg(sp); j++)
    { resp[j] += y;
      y *= u;
    }
  }
  for (j=0; j<=2*deg(sp); j++) resp[j] = resp[j]*(l1-l0)/(3*de_mint);
  return(LF_OK);
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
#include "locf.h"

extern int lf_status;
static double u[MXDIM], ilim[2*MXDIM], *ff, hh, *cff;
static lfdata *den_lfd;
static design *den_des;
static smpar *den_sp;
int fact[] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800};
int de_mint  = 20;
int de_itype = IDEFA;
int de_renorm= 0;

int multint(), prodint(), gausint(), mlinint();

#define NITYPE 7
static char *itype[NITYPE] = { "default", "multi", "product", "mlinear",
                          "hazard",  "sphere", "monte" };
static int   ivals[NITYPE] =
   { IDEFA, IMULT, IPROD, IMLIN, IHAZD, ISPHR, IMONT };
int deitype(char *z)
{ return(pmatch(z, itype, ivals, NITYPE, IDEFA));
}

void prresp(coef,resp,p)
double *coef, *resp;
int p;
{ int i, j;
  mut_printf("Coefficients:\n");
  for (i=0; i<p; i++) mut_printf("%8.5f ",coef[i]);
  mut_printf("\n");
  mut_printf("Response matrix:\n");
  for (i=0; i<p; i++)
  { for (j=0; j<p; j++) mut_printf("%9.6f, ",resp[i+j*p]);
    mut_printf("\n");
  }
}

int mif(u,d,resp,M)
double *u, *resp, *M;
int d;
{ double wt;
  int i, j, p;

  p = den_des->p;
  wt = weight(den_lfd, den_sp, u, NULL, hh, 0, 0.0);
  if (wt==0)
  { setzero(resp,p*p);
    return(p*p);
  }

  fitfun(den_lfd, den_sp, u,NULL,ff,NULL);
  if (link(den_sp)==LLOG)
    wt *= mut_exp(innerprod(ff,cff,p));
  for (i=0; i<p; i++)
    for (j=0; j<p; j++)
      resp[i*p+j] = wt*ff[i]*ff[j];
  return(p*p);
}

int multint(t,resp1,resp2,cf,h)
double *t, *resp1, *resp2, *cf, h;
{ int d, i, mg[MXDIM];

  if (ker(den_sp)==WGAUS) return(gausint(t,resp1,resp2,cf,h,den_lfd->sca));

  d = den_lfd->d;
  for (i=0; i<d; i++) mg[i] = de_mint;

  hh = h;
  cff= cf;
  simpsonm(mif,ilim,&ilim[d],d,resp1,mg,resp2);
  return(LF_OK);
}

int mlinint(t,resp1,resp2,cf,h)
double *t, *resp1, *resp2, *cf, h;
{
  double hd, nb, wt, wu, g[4], w0, w1, v, *sca;
  int d, p, i, j, jmax, k, l, z, jj[2];

  d = den_lfd->d; p = den_des->p; sca = den_lfd->sca;
  hd = 1;
  for (i=0; i<d; i++) hd *= h*sca[i];

  if (link(den_sp)==LIDENT)
  { setzero(resp1,p*p);
    resp1[0] = wint(d,NULL,0,ker(den_sp))*hd;
    if (deg(den_sp)==0) return(LF_OK);
    jj[0] = 2; w0 = wint(d,jj,1,ker(den_sp))*hd*h*h;
    for (i=0; i<d; i++) resp1[(i+1)*p+i+1] = w0*sca[i]*sca[i];
    if (deg(den_sp)==1) return(LF_OK);
    for (i=0; i<d; i++)
    { j = p-(d-i)*(d-i+1)/2;
      resp1[j] = resp1[p*j] = w0*sca[i]*sca[i]/2;
    }
    if (d>1)
    { jj[1] = 2;
      w0 = wint(d,jj,2,ker(den_sp)) * hd*h*h*h*h;
    }
    jj[0] = 4;
    w1 = wint(d,jj,1,ker(den_sp)) * hd*h*h*h*h/4;
    z = d+1;
    for (i=0; i<d; i++)
    { k = p-(d-i)*(d-i+1)/2;
      for (j=i; j<d; j++)
      { l = p-(d-j)*(d-j+1)/2;
        if (i==j) resp1[z*p+z] = w1*SQR(sca[i])*SQR(sca[i]);
        else
        { resp1[z*p+z] = w0*SQR(sca[i])*SQR(sca[j]);
          resp1[k*p+l] = resp1[k+p*l] = w0/4*SQR(sca[i])*SQR(sca[j]);
        }
        z++;
    } }
    return(LF_OK);
  }
  switch(deg(den_sp))
  { case 0:
      resp1[0] = mut_exp(cf[0])*wint(d,NULL,0,ker(den_sp))*hd;
      return(LF_OK);
    case 1:
      nb = 0.0;
      for (i=1; i<=d; i++)
      { v = h*cf[i]*sca[i-1];
        nb += v*v;
      }
      if (ker(den_sp)==WGAUS)
      { w0 = 1/(GFACT*GFACT);
        g[0] = mut_exp(cf[0]+w0*nb/2+d*log(S2PI/2.5));
        g[1] = g[3] = g[0]*w0;
        g[2] = g[0]*w0*w0;
      }
      else
      { wt = wu = mut_exp(cf[0]);
        w0 = wint(d,NULL,0,ker(den_sp)); g[0] = wt*w0;
        g[1] = g[2] = g[3] = 0.0;
        j = 0; jmax = (d+2)*de_mint;
        while ((j<jmax) && (wt*w0/g[0]>1.0e-8))
        { j++;
          jj[0] = 2*j; w0 = wint(d,jj,1,ker(den_sp));
          if (d==1) g[3] += wt * w0;
          else
          { jj[0] = 2; jj[1] = 2*j-2; w1 = wint(d,jj,2,ker(den_sp));
            g[3] += wt*w1;
            g[2] += wu*(w0-w1);
          }
          wt /= (2*j-1.0); g[1] += wt*w0;
          wt *= nb/(2*j); g[0] += wt*w0;
          wu /= (2*j-1.0)*(2*j);
          if (j>1) wu *= nb;
        }
        if (j==jmax) WARN(("mlinint: series not converged"));
      }
      g[0] *= hd; g[1] *= hd;
      g[2] *= hd; g[3] *= hd;
      resp1[0] = g[0];
      for (i=1; i<=d; i++)
      { resp1[i] = resp1[(d+1)*i] = cf[i]*SQR(h*sca[i-1])*g[1];
        for (j=1; j<=d; j++)
        { resp1[(d+1)*i+j] = (i==j) ? g[3]*SQR(h*sca[i-1]) : 0;
          resp1[(d+1)*i+j] += g[2]*SQR(h*h*sca[i-1]*sca[j-1])*cf[i]*cf[j];
        }
      }
      return(LF_OK);
  }
  LERR(("mlinint: deg=0,1 only"));
  return(LF_ERR);
}

void prodintresp(resp,prod_wk,dim,deg,p)
double *resp, prod_wk[MXDIM][2*MXDEG+1];
int dim, deg, p;
{ double prod;
  int i, j, k, j1, k1;

  prod = 1.0;
  for (i=0; i<dim; i++) prod *= prod_wk[i][0];
  resp[0] += prod;
  if (deg==0) return;

  for (j1=1; j1<=deg; j1++)
  { for (j=0; j<dim; j++)
    { prod = 1.0;
      for (i=0; i<dim; i++) prod *= prod_wk[i][j1*(j==i)];
      prod /= fact[j1];
      resp[1 + (j1-1)*dim +j] += prod;
    }
  }

  for (k1=1; k1<=deg; k1++)
    for (j1=k1; j1<=deg; j1++)
    { for (k=0; k<dim; k++)
        for (j=0; j<dim; j++)
        { prod = 1.0;
          for (i=0; i<dim; i++) prod *= prod_wk[i][k1*(k==i) + j1*(j==i)];
          prod /= fact[k1]*fact[j1];
          resp[ (1+(k1-1)*dim+k)*p + 1+(j1-1)*dim+j] += prod;
        }
    }
}

int prodint(t,resp,resp2,coef,h)
double *t, *resp, *resp2, *coef, h;
{ int dim, p, i, j, k, st;
  double cf[MXDEG+1], hj, hs, prod_wk[MXDIM][2*MXDEG+1];

  dim = den_lfd->d;
  p = den_des->p;
  for (i=0; i<p*p; i++) resp[i] = 0.0;
  cf[0] = coef[0];

/*  compute the one dimensional terms
 */
  for (i=0; i<dim; i++)
  { hj = 1; hs = h*den_lfd->sca[i];
    for (j=0; j<deg(den_sp); j++)
    { hj *= hs;
      cf[j+1] = hj*coef[ j*dim+i+1 ];
    }
    st = onedint(den_sp,cf,ilim[i]/hs,ilim[i+dim]/hs,prod_wk[i]);
    if (st==LF_BADP) return(st);
    hj = 1;
    for (j=0; j<=2*deg(den_sp); j++)
    { hj *= hs;
      prod_wk[i][j] *= hj;
    }
    cf[0] = 0.0; /* so we only include it once, when d>=2 */
  }

/*  transfer to the resp array
 */
  prodintresp(resp,prod_wk,dim,deg(den_sp),p);

/* Symmetrize.
*/
  for (k=0; k<p; k++)
    for (j=k; j<p; j++)
      resp[j*p+k] = resp[k*p+j];

  return(st);
}

int gausint(t,resp,C,cf,h,sca)
double *t, *resp, *C, *cf, h, *sca;
{ double nb, det, z, *P;
  int d, p, i, j, k, l, m1, m2, f;
  d = den_lfd->d; p = den_des->p;
  m1 = d+1; nb = 0;
  P = &C[d*d];
  resp[0] = 1;
  for (i=0; i<d; i++)
  { C[i*d+i] = SQR(GFACT/(h*sca[i]))-cf[m1++];
    for (j=i+1; j<d; j++) C[i*d+j] = C[j*d+i] = -cf[m1++];
  }
  eig_dec(C,P,d);
  det = 1;
  for (i=1; i<=d; i++)
  { det *= C[(i-1)*(d+1)];
    if (det <= 0) return(LF_BADP);
    resp[i] = cf[i];
    for (j=1; j<=d; j++) resp[j+i*p] = 0;
    resp[i+i*p] = 1;
    svdsolve(&resp[i*p+1],u,P,C,P,d,0.0);
  }
  svdsolve(&resp[1],u,P,C,P,d,0.0);
  det = sqrt(det);
  for (i=1; i<=d; i++)
  { nb += cf[i]*resp[i];
    resp[i*p] = resp[i];
    for (j=1; j<=d; j++)
      resp[i+p*j] += resp[i]*resp[j];
  }
  m1 = d;
  for (i=1; i<=d; i++)
    for (j=i; j<=d; j++)
    { m1++; f = 1+(i==j);
      resp[m1] = resp[m1*p] = resp[i*p+j]/f;
      m2 = d;
      for (k=1; k<=d; k++)
      { resp[m1+k*p] = resp[k+m1*p] =
        ( resp[i]*resp[j*p+k] + resp[j]*resp[i*p+k]
        + resp[k]*resp[i*p+j] - 2*resp[i]*resp[j]*resp[k] )/f;
        for (l=k; l<=d; l++)
        { m2++; f = (1+(i==j))*(1+(k==l));
          resp[m1+m2*p] = resp[m2+m1*p] = ( resp[i+j*p]*resp[k+l*p]
            + resp[i+k*p]*resp[j+l*p] + resp[i+l*p]*resp[j+k*p]
            - 2*resp[i]*resp[j]*resp[k]*resp[l] )/f;
    } } }
  z = mut_exp(d*0.918938533+cf[0]+nb/2)/det;
  multmatscal(resp,z,p*p);
  return(LF_OK);
}

int likeden(coef, lk0, f1, A)
double *coef, *lk0, *f1, *A;
{ double lk, r;
  int i, j, p, rstat;

  lf_status = LF_OK;
  p = den_des->p;
  if ((link(den_sp)==LIDENT) && (coef[0] != 0.0)) return(NR_BREAK);
  lf_status = (den_des->itype)(den_des->xev,A,den_des->xtwx.Q,coef,den_des->h);
  if (lf_error) lf_status = LF_ERR;
  if (lf_status==LF_BADP)
  { *lk0 = -1.0e300;
    return(NR_REDUCE);
  }
  if (lf_status!=LF_OK) return(NR_BREAK);
  if (lf_debug>2) prresp(coef,A,p);

  den_des->xtwx.p = p;
  rstat = NR_OK;
  switch(link(den_sp))
  { case LLOG:
      r = den_des->ss[0]/A[0];
      coef[0] += log(r);
      multmatscal(A,r,p*p);
      A[0] = den_des->ss[0];
      lk = -A[0];
      if (fabs(coef[0]) > 700)
      { lf_status = LF_OOB;
        rstat = NR_REDUCE;
      }
      for (i=0; i<p; i++)
      { lk += coef[i]*den_des->ss[i];
        f1[i] = den_des->ss[i]-A[i];
      }
      break;
    case LIDENT:
      lk = 0.0;
      for (i=0; i<p; i++)
      { f1[i] = den_des->ss[i];
        for (j=0; j<p; j++)
          den_des->res[i] -= A[i*p+j]*coef[j];
      }
      break;
  }
  *lk0 = den_des->llk = lk;

  return(rstat);
}

int inre(x,bound,d)
double *x, *bound;
int d;
{ int i, z;
  z = 1;
  for (i=0; i<d; i++)
    if (bound[i]<bound[i+d])
      z &= (x[i]>=bound[i]) & (x[i]<=bound[i+d]);
  return(z);
}

int setintlimits(lfd, x, h, ang, lset)
lfdata *lfd;
int *ang, *lset;
double *x, h;
{ int d, i;
  d = lfd->d;
  *ang = *lset = 0;
  for (i=0; i<d; i++)
  { if (lfd->sty[i]==STANGL)
    { ilim[i+d] = ((h<2) ? 2*asin(h/2) : PI)*lfd->sca[i];
      ilim[i] = -ilim[i+d];
      *ang = 1;
    }
    else
    { ilim[i+d] = h*lfd->sca[i];
      ilim[i] = -ilim[i+d];

      if (lfd->sty[i]==STLEFT) { ilim[i+d] = 0; *lset = 1; }
      if (lfd->sty[i]==STRIGH) { ilim[i] = 0;   *lset = 1; }

      if (lfd->xl[i]<lfd->xl[i+d]) /* user limits for this variable */
      { if (lfd->xl[i]-x[i]> ilim[i])
        { ilim[i] = lfd->xl[i]-x[i]; *lset=1; }
        if (lfd->xl[i+d]-x[i]< ilim[i+d])
        { ilim[i+d] = lfd->xl[i+d]-x[i]; *lset=1; }
      }
    }
    if (ilim[i]==ilim[i+d]) return(LF_DEMP); /* empty integration */
  }
  return(LF_OK);
}

int selectintmeth(itype,lset,ang)
int itype, lset, ang;
{
  if (itype==IDEFA) /* select the default method */
  { if (fam(den_sp)==THAZ)
    { if (ang) return(IDEFA);
      return( IHAZD );
    }

    if (ubas(den_sp)) return(IMULT);

    if (ang) return(IMULT);

    if (iscompact(ker(den_sp)))
    { if (kt(den_sp)==KPROD) return(IPROD);
      if (lset)
        return( (den_lfd->d==1) ? IPROD : IMULT );
      if (deg(den_sp)<=1) return(IMLIN);
      if (den_lfd->d==1) return(IPROD);
      return(IMULT);
    }

    if (ker(den_sp)==WGAUS)
    { if (lset) WARN(("Integration for Gaussian weights ignores limits"));
      if ((den_lfd->d==1)|(kt(den_sp)==KPROD)) return(IPROD);
      if (deg(den_sp)<=1) return(IMLIN);
      if (deg(den_sp)==2) return(IMULT);
    }

    return(IDEFA);
  }

  /* user provided an integration method, check it is valid */

  if (fam(den_sp)==THAZ)
  { if (ang) return(INVLD);
    if (!iscompact(ker(den_sp))) return(INVLD);
    return( ((kt(den_sp)==KPROD) | (kt(den_sp)==KSPH)) ? IHAZD : INVLD );
  }

  if ((ang) && (itype != IMULT)) return(INVLD);

  switch(itype)
  { case IMULT:
      if (ker(den_sp)==WGAUS) return(deg(den_sp)==2);
      return( iscompact(ker(den_sp)) ? IMULT : INVLD );
    case IPROD: return( ((den_lfd->d==1) | (kt(den_sp)==KPROD)) ? IPROD : INVLD );
    case IMLIN: return( ((kt(den_sp)==KSPH) && (!lset) &&
      (deg(den_sp)<=1)) ? IMLIN : INVLD );
  }

  return(INVLD);
}

extern double lf_tol;

int densinit(lfd,des,sp)
lfdata *lfd;
design *des;
smpar *sp;
{ int p, i, ii, j, nnz, rnz, ang, lset, status;
  double w, *cf;

  den_lfd = lfd;
  den_des = des;
  den_sp  = sp;
  cf = des->cf;

  lf_tol = (link(sp)==LLOG) ? 1.0e-6 : 0.0;

  p = des->p;
  ff = des->xtwx.wk;
  cf[0] = NOSLN;
  for (i=1; i<p; i++) cf[i] = 0.0;

  if (!inre(des->xev,lfd->xl,lfd->d)) return(LF_XOOR);

  status = setintlimits(lfd,des->xev,des->h,&ang,&lset);
  if (status != LF_OK) return(status);

  switch(selectintmeth(de_itype,lset,ang))
  { case IMULT: des->itype = multint; break;
    case IPROD: des->itype = prodint; break;
    case IMLIN: des->itype = mlinint; break;
    case IHAZD: des->itype = hazint; break;
    case INVLD: LERR(("Invalid integration method %d",de_itype));
                break;
    case IDEFA: LERR(("No integration type available for this model"));
                break;
    default: LERR(("densinit: unknown integral type"));
  }

  switch(deg(den_sp))
  { case 0: rnz = 1; break;
    case 1: rnz = 1; break;
    case 2: rnz = lfd->d+1; break;
    case 3: rnz = lfd->d+2; break;
    default: LERR(("densinit: invalid degree %d",deg(den_sp)));
  }
  if (lf_error) return(LF_ERR);

  setzero(des->ss,p);
  nnz = 0;
  for (i=0; i<des->n; i++)
  { ii = des->ind[i];
    if (!cens(lfd,ii))
    { w = wght(des,ii)*prwt(lfd,ii);
      for (j=0; j<p; j++) des->ss[j] += d_xij(des,ii,j)*w;
      if (wght(des,ii)>0.00001) nnz++;
  } }

  if (fam(den_sp)==THAZ) haz_init(lfd,des,sp,ilim);
/* this should really only be done once. Not sure how to enforce that,
 * esp. when locfit() has been called directly.
 */
  if (fam(den_sp)==TDEN)
    des->smwt = (lfd->w==NULL) ? lfd->n : vecsum(lfd->w,lfd->n);

  if (lf_debug>2)
  { mut_printf("    LHS: ");
    for (i=0; i<p; i++) mut_printf(" %8.5f",des->ss[i]);
    mut_printf("\n");
  }

  switch(link(den_sp))
  { case LIDENT:
      cf[0] = 0.0;
      return(LF_OK);
    case LLOG:
      if (nnz<rnz) { cf[0] = -1000; return(LF_DNOP); }
      cf[0] = 0.0;
      return(LF_OK);
    default:
      LERR(("unknown link in densinit"));
      return(LF_ERR);
  }
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
#include "locf.h"

int bino_vallink(link)
int link;
{ return((link==LLOGIT) | (link==LIDENT) | (link==LASIN));
}

int bino_fam(y,p,th,link,res,cens,w)
double y, p, th, *res, w;
int link, cens;
{ double wp;
  if (link==LINIT)
  { if (y<0) y = 0;
    if (y>w) y = w;
    res[ZDLL] = y;
    return(LF_OK);
  }
  wp = w*p;
  if (link==LIDENT)
  { if ((p<=0) && (y>0)) return(LF_BADP);
    if ((p>=1) && (y<w)) return(LF_BADP);
    res[ZLIK] = res[ZDLL] = res[ZDDLL] = 0.0;
    if (y>0)
    { res[ZLIK] += y*log(wp/y);
      res[ZDLL] += y/p;
      res[ZDDLL]+= y/(p*p);
    }
    if (y<w)
    { res[ZLIK] += (w-y)*log((w-wp)/(w-y));
      res[ZDLL] -= (w-y)/(1-p);
      res[ZDDLL]+= (w-y)/SQR(1-p);
    }
    return(LF_OK);
  }
  if (link==LLOGIT)
  { if ((y<0) | (y>w)) /* goon observation; delete it */
    { res[ZLIK] = res[ZDLL] = res[ZDDLL] = 0.0;
      return(LF_OK);
    }
    res[ZLIK] = (th<0) ? th*y-w*log(1+exp(th)) : th*(y-w)-w*log(1+exp(-th));
    if (y>0) res[ZLIK] -= y*log(y/w);
    if (y<w) res[ZLIK] -= (w-y)*log(1-y/w);
    res[ZDLL] = (y-wp);
    res[ZDDLL]= wp*(1-p);
    return(LF_OK);
  }
  if (link==LASIN)
  { if ((p<=0) && (y>0)) return(LF_BADP);
    if ((p>=1) && (y<w)) return(LF_BADP);
    if ((th<0) | (th>PI/2)) return(LF_BADP);
    res[ZDLL] = res[ZDDLL] = res[ZLIK] = 0;
    if (y>0)
    { res[ZDLL] += 2*y*sqrt((1-p)/p);
      res[ZLIK] += y*log(wp/y);
    }
    if (y<w)
    { res[ZDLL] -= 2*(w-y)*sqrt(p/(1-p));
      res[ZLIK] += (w-y)*log((w-wp)/(w-y));
    }
    res[ZDDLL] = 4*w;
    return(LF_OK);
  }
  LERR(("link %d invalid for binomial family",link));
  return(LF_LNK);
}

int bino_check(sp,des,lfd)
smpar *sp;
design *des;
lfdata *lfd;
{ int i, ii;
  double t0, t1;

  if (fabs(des->cf[0])>700) return(LF_OOB);

  /* check for separation.
   * this won't detect separation if there's boundary points with
   *   both 0 and 1 responses.
   */
  t0 = -1e100; t1 = 1e100;
  for (i=0; i<des->n; i++)
  { ii = des->ind[i];
    if ((resp(lfd,ii)<prwt(lfd,ii)) && (fitv(des,ii) > t0)) t0 = fitv(des,ii);
    if ((resp(lfd,ii)>0) && (fitv(des,ii) < t1)) t1 = fitv(des,ii);
    if (t1 <= t0) return(LF_OK);
  }
  mut_printf("separated %8.5f %8.5f\n",t0,t1);
  return(LF_NSLN);
}

void setfbino(fam)
family *fam;
{ fam->deflink = LLOGIT;
  fam->canlink = LLOGIT;
  fam->vallink = bino_vallink;
  fam->family  = bino_fam;
  fam->pcheck  = bino_check;
}

int rbin_vallink(link)
int link;
{ return(link==LLOGIT);
}

int rbin_fam(y,p,th,link,res,cens,w)
double y, p, th, *res, w;
int link, cens;
{ double s2y;
  if (link==LINIT)
  { res[ZDLL] = y;
    return(LF_OK);
  }
  if ((y<0) | (y>w)) /* goon observation; delete it */
  { res[ZLIK] = res[ZDLL] = res[ZDDLL] = 0.0;
    return(LF_OK);
  }
  res[ZLIK] = (th<0) ? th*y-w*log(1+exp(th)) : th*(y-w)-w*log(1+exp(-th));
  if (y>0) res[ZLIK] -= y*log(y/w);
  if (y<w) res[ZLIK] -= (w-y)*log(1-y/w);
  res[ZDLL] = (y-w*p);
  res[ZDDLL]= w*p*(1-p);
  if (-res[ZLIK]>HUBERC*HUBERC/2.0)
  { s2y = sqrt(-2*res[ZLIK]);
    res[ZLIK] = HUBERC*(HUBERC/2.0-s2y);
    res[ZDLL] *= HUBERC/s2y;
    res[ZDDLL] = HUBERC/s2y*(res[ZDDLL]-1/(s2y*s2y)*w*p*(1-p));
  }
  return(LF_OK);
}

void setfrbino(fam)
family *fam;
{ fam->deflink = LLOGIT;
  fam->canlink = LLOGIT;
  fam->vallink = rbin_vallink;
  fam->family  = rbin_fam;
  fam->pcheck  = bino_check;
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
#include "locf.h"

int circ_vallink(link)
int link;
{ return(link==LIDENT);
}

int circ_fam(y,mean,th,link,res,cens,w)
double y, mean, th, *res, w;
int link, cens;
{ if (link==LINIT)
  { res[ZDLL] = w*sin(y);
    res[ZLIK] = w*cos(y);
    return(LF_OK);
  }
  res[ZDLL] = w*sin(y-mean);
  res[ZDDLL]= w*cos(y-mean);
  res[ZLIK] = res[ZDDLL]-w;
  return(LF_OK);
}

extern double lf_tol;
int circ_init(lfd,des,sp)
lfdata *lfd;
design *des;
smpar *sp;
{ int i, ii;
  double s0, s1;
  s0 = s1 = 0.0;
  for (i=0; i<des->n; i++)
  { ii = des->ind[i];
    s0 += wght(des,ii)*prwt(lfd,ii)*sin(resp(lfd,ii)-base(lfd,ii));
    s1 += wght(des,ii)*prwt(lfd,ii)*cos(resp(lfd,ii)-base(lfd,ii));
  }
  des->cf[0] = atan2(s0,s1);
  for (i=1; i<des->p; i++) des->cf[i] = 0.0;
  lf_tol = 1.0e-6;
  return(LF_OK);
}


void setfcirc(fam)
family *fam;
{ fam->deflink = LIDENT;
  fam->canlink = LIDENT;
  fam->vallink = circ_vallink;
  fam->family  = circ_fam;
  fam->initial = circ_init;
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
#include "locf.h"

int dens_vallink(link)
int link;
{ return((link==LIDENT) | (link==LLOG));
}

int dens_fam(y,mean,th,link,res,cens,w)
double y, mean, th, *res, w;
int link, cens;
{ if (cens)
    res[ZLIK] = res[ZDLL] = res[ZDDLL] = 0.0;
  else
  { res[ZLIK] = w*th;
    res[ZDLL] = res[ZDDLL] = w;
  }
  return(LF_OK);
}

void setfdensity(fam)
family *fam;
{ fam->deflink = LLOG;
  fam->canlink = LLOG;
  fam->vallink = dens_vallink;
  fam->family  = dens_fam;
  fam->initial = densinit;
  fam->like = likeden;
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
#include "locf.h"

int gamma_vallink(link)
int link;
{ return((link==LIDENT) | (link==LLOG) | (link==LINVER));
}

int gamma_fam(y,mean,th,link,res,cens,w)
double y, mean, th, *res, w;
int link, cens;
{ double lb, pt, dg;
  if (link==LINIT)
  { res[ZDLL] = MAX(y,0.0);
    return(LF_OK);
  }
  res[ZLIK] = res[ZDLL] = res[ZDDLL] = 0.0;
  if (w==0.0) return(LF_OK);
  if ((mean<=0) & (y>0)) return(LF_BADP);
  if (link==LIDENT) lb = 1/th;
  if (link==LINVER) lb = th;
  if (link==LLOG)   lb = mut_exp(-th);
  if (cens)
  { if (y<=0) return(LF_OK);
    pt = 1-igamma(lb*y,w);
    dg = dgamma(lb*y,w,1.0,0);
    res[ZLIK] = log(pt);
    res[ZDLL] = -y*dg/pt;
/*
 * res[ZDLL]  = -y*dg/pt * dlb/dth.
 * res[ZDDLL] =  y*dg/pt * (d2lb/dth2 + ((w-1)/lb-y)*(dlb/dth)^2)
 *              + res[ZDLL]^2.
 */
    if (link==LLOG)       /* lambda = exp(-theta) */
    { res[ZDLL] *= -lb;
      res[ZDDLL] = dg*y*lb*(w-lb*y)/pt + SQR(res[ZDLL]);
      return(LF_OK);
    }
    if (link==LINVER)     /* lambda = theta */
    { res[ZDLL] *= 1.0;
      res[ZDDLL] = dg*y*((w-1)*mean-y)/pt + SQR(res[ZDLL]);
      return(LF_OK);
    }
    if (link==LIDENT)     /* lambda = 1/theta */
    { res[ZDLL] *= -lb*lb;
      res[ZDDLL] = dg*y*lb*lb*lb*(1+w-lb*y)/pt + SQR(res[ZDLL]);
      return(LF_OK);
    }
  }
  else
  { if (y<0) WARN(("Negative Gamma observation"));
    if (link==LLOG)
    { res[ZLIK] = -lb*y+w*(1-th);
      if (y>0) res[ZLIK] += w*log(y/w);
      res[ZDLL] = lb*y-w;
      res[ZDDLL]= lb*y;
      return(LF_OK);
    }
    if (link==LINVER)
    { res[ZLIK] = -lb*y+w-w*log(mean);
      if (y>0) res[ZLIK] += w*log(y/w);
      res[ZDLL] = -y+w*mean;
      res[ZDDLL]= w*mean*mean;
      return(LF_OK);
    }
    if (link==LIDENT)
    { res[ZLIK] = -lb*y+w-w*log(mean);
      if (y>0) res[ZLIK] += w*log(y/w);
      res[ZDLL] = lb*lb*(y-w*mean);
      res[ZDDLL]= lb*lb*lb*(2*y-w*mean);
      return(LF_OK);
    }
  }
  LERR(("link %d invalid for Gamma family",link));
  return(LF_LNK);
}

void setfgamma(fam)
family *fam;
{ fam->deflink = LLOG;
  fam->canlink = LINVER;
  fam->vallink = gamma_vallink;
  fam->family  = gamma_fam;
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
#include "locf.h"

int gaus_vallink(link)
int link;
{ return((link==LIDENT) | (link==LLOG) | (link==LLOGIT));
}

int gaus_fam(y,mean,th,link,res,cens,w)
double y, mean, th, *res, w;
int link, cens;
{ double z, pz, dp;
  if (link==LINIT)
  { res[ZDLL] = w*y;
    return(LF_OK);
  }
  z = y-mean;
  if (cens)
  { if (link!=LIDENT)
    { LERR(("Link invalid for censored Gaussian family"));
      return(LF_LNK);
    }
    pz = mut_pnorm(-z);
    dp = ((z>6) ? ptail(-z) : exp(-z*z/2)/pz)/2.5066283;
    res[ZLIK] = w*log(pz);
    res[ZDLL] = w*dp;
    res[ZDDLL]= w*dp*(dp-z);
    return(LF_OK);
  }
  res[ZLIK] = -w*z*z/2; 
  switch(link)
  { case LIDENT:
      res[ZDLL] = w*z;
      res[ZDDLL]= w;
      break;
    case LLOG:
      res[ZDLL] = w*z*mean;
      res[ZDDLL]= w*mean*mean;
      break;
    case LLOGIT:
      res[ZDLL] = w*z*mean*(1-mean);
      res[ZDDLL]= w*mean*mean*(1-mean)*(1-mean);
      break;
    default:
      LERR(("Invalid link for Gaussian family"));
      return(LF_LNK);
  }
  return(LF_OK);
}

int gaus_check(sp,des,lfd)
smpar *sp;
design *des;
lfdata *lfd;
{ int i, ii;
  if (fami(sp)->robust) return(LF_OK);
  if (link(sp)==LIDENT)
  { for (i=0; i<des->n; i++)
    { ii = des->ind[i];
      if (cens(lfd,ii)) return(LF_OK);
    }
    return(LF_DONE);
  }
  return(LF_OK);
}

void setfgauss(fam)
family *fam;
{ fam->deflink = LIDENT;
  fam->canlink = LIDENT;
  fam->vallink = gaus_vallink;
  fam->family  = gaus_fam;
  fam->pcheck  = gaus_check;
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
#include "locf.h"

int geom_vallink(link)
int link;
{ return((link==LIDENT) | (link==LLOG));
}

int geom_fam(y,mean,th,link,res,cens,w)
double y, mean, th, *res, w;
int link, cens;
{ double p, pt, dp, p1;
  if (link==LINIT)
  { res[ZDLL] = MAX(y,0.0);
    return(LF_OK);
  }
  p = 1/(1+mean);
  if (cens) /* censored observation */
  { if (y<=0)
    { res[ZLIK] = res[ZDLL] = res[ZDDLL] = 0;
      return(LF_OK);
    }
    p1 = (link==LIDENT) ? -p*p : -p*(1-p);
    pt = 1-ibeta(p,w,y);
    dp = dbeta(p,w,y,0)/pt;
    res[ZLIK] = log(pt);
    res[ZDLL] = -dp*p1;
    res[ZDDLL] = dp*dp*p1*p1;
    if (link==LIDENT)
      res[ZDDLL] += dp*p*p*p*(1+w*(1-p)-p*y)/(1-p);
    else
      res[ZDDLL] += dp*p*(1-p)*(w*(1-p)-p*y);
    return(LF_OK);
  }
  else
  { res[ZLIK] = (y+w)*log((y/w+1)/(mean+1));
    if (y>0) res[ZLIK] += y*log(w*mean/y);
    if (link==LLOG)
    { res[ZDLL] = (y-w*mean)*p;
      res[ZDDLL]= (y+w)*p*(1-p);
      return(LF_OK);
    }
    if (link==LIDENT)
    { res[ZDLL] = (y-w*mean)/(mean*(1+mean));
      res[ZDDLL]= w/(mean*(1+mean));
      return(LF_OK);
    }
  }
  LERR(("link %d invalid for geometric family",link));
  return(LF_LNK);
}

void setfgeom(fam)
family *fam;
{ fam->deflink = LLOG;
  fam->canlink = LIDENT; /* this isn't correct. I haven't prog. canon */
  fam->vallink = geom_vallink;
  fam->family  = geom_fam;
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
#include "locf.h"

#define HUBERC 2.0

double links_rs;
int inllmix=0;

/*
 * lffamily("name") converts family names into a numeric value.
 * typical usage is  fam(&lf->sp) = lffamily("gaussian");
 * Note that family can be preceded by q and/or r for quasi, robust.
 *
 * link(&lf->sp) = lflink("log") does the same for the link function.
 */
#define NFAMILY 18
static char *famil[NFAMILY] =
  { "density", "ate",   "hazard",    "gaussian", "binomial",
    "poisson", "gamma", "geometric", "circular", "obust", "huber",
    "weibull", "cauchy","probab",    "logistic", "nbinomial",
    "vonmises", "quant" };
static int   fvals[NFAMILY] =
  { TDEN,  TRAT,  THAZ,  TGAUS, TLOGT,
    TPOIS, TGAMM, TGEOM, TCIRC, TROBT, TROBT,
    TWEIB, TCAUC, TPROB, TLOGT, TGEOM, TCIRC, TQUANT };
int lffamily(z)
char *z;
{ int quasi, robu, f;
  quasi = robu = 0;
  while ((z[0]=='q') | (z[0]=='r'))
  { quasi |= (z[0]=='q');
    robu  |= (z[0]=='r');
    z++;
  }
  z[0] = tolower(z[0]);
  f = pmatch(z,famil,fvals,NFAMILY,-1);
  if ((z[0]=='o') | (z[0]=='a')) robu = 0;
  if (f==-1)
  { WARN(("unknown family %s",z));
    f = TGAUS;
  }
  if (quasi) f += 64;
  if (robu)  f += 128;
  return(f);
}

#define NLINKS 8
static char *ltype[NLINKS] = { "default", "canonical", "identity", "log",
                          "logi",    "inverse",   "sqrt",     "arcsin" };
static int   lvals[NLINKS] = { LDEFAU, LCANON, LIDENT, LLOG,
                          LLOGIT, LINVER, LSQRT,  LASIN };
int lflink(char *z)
{ int f;
  if (z==NULL) return(LDEFAU);
  z[0] = tolower(z[0]);
  f = pmatch(z, ltype, lvals, NLINKS, -1);
  if (f==-1)
  { WARN(("unknown link %s",z));
    f = LDEFAU;
  }
  return(f);
}

int defaultlink(link,fam)
int link;
family *fam;
{ if (link==LDEFAU) return(fam->deflink);
  if (link==LCANON) return(fam->canlink);
  return(link);
}

/*
void robustify(res,rs)
double *res, rs;
{ double sc, z;
  sc = rs*HUBERC;
  if (res[ZLIK] > -sc*sc/2) return;
  z = sqrt(-2*res[ZLIK]);
  res[ZDDLL]= -sc*res[ZDLL]*res[ZDLL]/(z*z*z)+sc*res[ZDDLL]/z;
  res[ZDLL]*= sc/z;
  res[ZLIK] = sc*sc/2-sc*z;
}
*/
void robustify(res,rs)
double *res, rs;
{ double sc, z;
  sc = rs*HUBERC;
  if (res[ZLIK] > -sc*sc/2)
  { res[ZLIK] /= sc*sc;
    res[ZDLL] /= sc*sc;
    res[ZDDLL] /= sc*sc;
    return;
  }
  z = sqrt(-2*res[ZLIK]);
  res[ZDDLL]= (-sc*res[ZDLL]*res[ZDLL]/(z*z*z)+sc*res[ZDDLL]/z)/(sc*sc);
  res[ZDLL]*= 1.0/(z*sc);
  res[ZLIK] = 0.5-z/sc;
}

double lf_link(y,lin)
double y;
int lin;
{ switch(lin)
  { case LIDENT: return(y);
    case LLOG:   return(log(y));
    case LLOGIT: return(logit(y));
    case LINVER: return(1/y);
    case LSQRT:  return(sqrt(fabs(y)));
    case LASIN:  return(asin(sqrt(y)));
  }
  LERR(("link: unknown link %d",lin));
  return(0.0);
}

double invlink(th,lin)
double th;
int lin;
{ switch(lin)
  { case LIDENT: return(th);
    case LLOG:   return(mut_exp(th));
    case LLOGIT: return(expit(th));
    case LINVER: return(1/th);
    case LSQRT:  return(th*fabs(th));
    case LASIN:  return(sin(th)*sin(th));
    case LINIT:  return(0.0);
  }
  LERR(("invlink: unknown link %d",lin));
  return(0.0);
}

/* the link and various related functions */
int links(th,y,fam,link,res,c,w,rs)
double th, y, *res, w, rs;
int link, c;
family *fam;
{ double mean;
  int st;

  mean = res[ZMEAN] = invlink(th,link);
  if (lf_error) return(LF_LNK);
  links_rs = rs;
/*  mut_printf("links: rs %8.5f\n",rs); */

  st = fam->family(y,mean,th,link,res,c,w);

  if (st!=LF_OK) return(st);
  if (link==LINIT) return(st);
  if (isrobust(fam)) robustify(res,rs);
  return(st);
}

/*
  stdlinks is a version of links when family, link, response e.t.c
  all come from the standard places.
*/
int stdlinks(res,lfd,sp,i,th,rs)
lfdata *lfd;
smpar *sp;
double th, rs, *res;
int i;
{
  return(links(th,resp(lfd,i),fami(sp),link(sp),res,cens(lfd,i),prwt(lfd,i),rs));
}

/*
 *  functions used in variance, skewness, kurtosis calculations
 *  in scb corrections.
 */

double b2(th,tg,w)
double th, w;
int tg;
{ double y;
  switch(tg&63)
  { case TGAUS: return(w);
    case TPOIS: return(w*mut_exp(th));
    case TLOGT:
      y = expit(th);
      return(w*y*(1-y));
  }
  LERR(("b2: invalid family %d",tg));
  return(0.0);
}

double b3(th,tg,w)
double th, w;
int tg;
{ double y;
  switch(tg&63)
  { case TGAUS: return(0.0);
    case TPOIS: return(w*mut_exp(th));
    case TLOGT:
      y = expit(th);
      return(w*y*(1-y)*(1-2*y));
  }
  LERR(("b3: invalid family %d",tg));
  return(0.0);
}

double b4(th,tg,w)
double th, w;
int tg;
{ double y;
  switch(tg&63)
  { case TGAUS: return(0.0);
    case TPOIS: return(w*mut_exp(th));
    case TLOGT:
      y = expit(th); y = y*(1-y);
      return(w*y*(1-6*y));
  }
  LERR(("b4: invalid family %d",tg));
  return(0.0);
}

int def_check(sp,des,lfd)
smpar *sp;
design *des;
lfdata *lfd;
{ switch(link(sp))
  { case LLOG: if (des->cf[0]>700) return(LF_OOB);
               break;
  }
  return(LF_OK);
}
extern void setfdensity(), setfgauss(), setfbino(), setfpoisson();
extern void setfgamma(), setfgeom(), setfcirc(), setfweibull();
extern void setfrbino(), setfrobust(), setfcauchy(), setfquant();

void setfamily(sp)
smpar *sp;
{ int tg, lnk;
  family *f;

  tg = fam(sp);
  f = fami(sp);
  f->quasi = tg&64;
  f->robust = tg&128;
  f->initial = reginit;
  f->like = likereg;
  f->pcheck = def_check;

  switch(tg&63)
  { case TDEN:
    case THAZ:
    case TRAT:	setfdensity(f); break;
    case TGAUS: setfgauss(f); break;
    case TLOGT: setfbino(f); break;
    case TRBIN: setfrbino(f); break;
    case TPROB:
    case TPOIS: setfpoisson(f); break;
    case TGAMM: setfgamma(f); break;
    case TGEOM: setfgeom(f); break;
    case TWEIB: setfweibull(f);
    case TCIRC: setfcirc(f); break;
    case TROBT: setfrobust(f); break;
    case TCAUC: setfcauchy(f); break;
    case TQUANT: setfquant(f); break;
    default: LERR(("setfamily: unknown family %d",tg&63));
             return;
  }
  
  lnk = defaultlink(link(sp),f);
  if (!f->vallink(lnk))
  { WARN(("setfamily: invalid link %d - revert to default",link(sp)));
    link(sp) = f->deflink;
  }
  else
    link(sp) = lnk;
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
#include "locf.h"

int pois_vallink(link)
int link;
{ return((link==LLOG) | (link==LIDENT) | (link==LSQRT));
}

int pois_fam(y,mean,th,link,res,cens,w)
double y, mean, th, *res, w;
int link, cens;
{ double wmu, pt, dp;
  if (link==LINIT)
  { res[ZDLL] = MAX(y,0.0);
    return(LF_OK);
  }
  wmu = w*mean;
  if (inllmix) y = w*y;
  if (cens)
  { if (y<=0)
    { res[ZLIK] = res[ZDLL] = res[ZDDLL] = 0.0;
      return(LF_OK);
    }
    pt = igamma(wmu,y);
    dp = dgamma(wmu,y,1.0,0)/pt;
    res[ZLIK] = log(pt);
/*
 * res[ZDLL] = dp * w*dmu/dth
 * res[ZDDLL]= -dp*(w*d2mu/dth2 + (y-1)/mu*(dmu/dth)^2) + res[ZDLL]^2
 */
    if (link==LLOG)
    { res[ZDLL] = dp*wmu;
      res[ZDDLL]= -dp*wmu*(y-wmu) + SQR(res[ZDLL]);
      return(LF_OK);
    }
    if (link==LIDENT)
    { res[ZDLL] = dp*w;
      res[ZDDLL]= -dp*(y-1-wmu)*w/mean + SQR(res[ZDLL]);
      return(LF_OK);
    }
    if (link==LSQRT)
    { res[ZDLL] = dp*2*w*th;
      res[ZDDLL]= -dp*w*(4*y-2-4*wmu) + SQR(res[ZDLL]);
      return(LF_OK);
  } }
  if (link==LLOG)
  { if (y<0) /* goon observation - delete it */
    { res[ZLIK] = res[ZDLL] = res[ZDDLL] = 0;
      return(LF_OK);
    }
    res[ZLIK] = res[ZDLL] = y-wmu;
    if (y>0) res[ZLIK] += y*(th-log(y/w));
    res[ZDDLL] = wmu;
    return(LF_OK);
  }
  if (link==LIDENT)
  { if ((mean<=0) && (y>0)) return(LF_BADP);
    res[ZLIK] = y-wmu;
    res[ZDLL] = -w;
    res[ZDDLL] = 0;
    if (y>0)
    { res[ZLIK] += y*log(wmu/y);
      res[ZDLL] += y/mean;
      res[ZDDLL]= y/(mean*mean);
    }
    return(LF_OK);
  }
  if (link==LSQRT)
  { if ((mean<=0) && (y>0)) return(LF_BADP);
    res[ZLIK] = y-wmu;
    res[ZDLL] = -2*w*th;
    res[ZDDLL]= 2*w;
    if (y>0)
    { res[ZLIK] += y*log(wmu/y);
      res[ZDLL] += 2*y/th;
      res[ZDDLL]+= 2*y/mean;
    }
    return(LF_OK);
  }
  LERR(("link %d invalid for Poisson family",link));
  return(LF_LNK);
}

void setfpoisson(fam)
family *fam;
{ fam->deflink = LLOG;
  fam->canlink = LLOG;
  fam->vallink = pois_vallink;
  fam->family  = pois_fam;
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
#include "locf.h"

#define QTOL 1.0e-10
extern int lf_status;
static double q0;

int quant_vallink(int link) { return(1); }

int quant_fam(y,mean,th,link,res,cens,w)
double y, mean, th, *res, w;
int link, cens;
{ double z, p;
  if (link==LINIT)
  { res[ZDLL] = w*y;
    return(LF_OK);
  }
p = 0.5; /* should be pen(sp) */
  z = y-mean;
  res[ZLIK] = (z<0) ? (w*z/p) : (-w*z/(1-p));
  res[ZDLL] = (z<0) ? -w/p : w/(1-p);
  res[ZDDLL]= w/(p*(1-p));
  return(LF_OK);
}

int quant_check(sp,des,lfd)
smpar *sp;
design *des;
lfdata *lfd;
{ return(LF_DONE);
}

void setfquant(fam)
family *fam;
{ fam->deflink = LIDENT;
  fam->canlink = LIDENT;
  fam->vallink = quant_vallink;
  fam->family  = quant_fam;
  fam->pcheck  = quant_check;
}

/*
 * cycling rule for choosing among ties.
 */
int tiecycle(ind,i0,i1,oi)
int *ind, i0, i1, oi;
{ int i, ii, im;
  im = ind[i0];
  for (i=i0+1; i<=i1; i++)
  { ii = ind[i];
    if (im<=oi)
    { if ((ii<im) | (ii>oi)) im = ii;
    }
    else
    { if ((ii<im) & (ii>oi)) im = ii;
    }
  }
  return(im);
}

/*
 * move coefficient vector cf, as far as possible, in direction dc.
 */
int movecoef(lfd,des,p,cf,dc,oi)
lfdata *lfd;
design *des;
double p, *cf, *dc;
int oi;
{ int i, ii, im, i0, i1, j;
  double *lb, *el, e, sp, sn, sw, sum1, sum2, tol1;

  lb = des->th;
  el = des->res;
  sum1 = sum2 = 0.0;

  sp = sn = sw = 0.0;
  for (i=0; i<des->n; i++)
  { ii = des->ind[i];
    lb[ii] = innerprod(dc,d_xi(des,ii),des->p);
    e = resp(lfd,ii) - innerprod(cf,d_xi(des,ii),des->p);
    el[ii] = (fabs(lb[ii])<QTOL) ? 1e100 : e/lb[ii];
    if (lb[ii]>0)
      sp += prwt(lfd,ii)*wght(des,ii)*lb[ii];
    else
      sn -= prwt(lfd,ii)*wght(des,ii)*lb[ii];
    sw += prwt(lfd,ii)*wght(des,ii);
  }
printf("sp %8.5f  sn %8.5f\n",sn,sp);
/* if sn, sp are both zero, should return an LF_PF.
 * but within numerical tolerance? what does it mean?
 */
  if (sn+sp <= QTOL*q0) { lf_status = LF_PF; return(0); }

  sum1 = sp/(1-p) + sn/p;
  tol1 = QTOL*(sp+sn);
  mut_order(el,des->ind,0,des->n-1);

  for (i=0; i<des->n; i++)
  { ii = des->ind[i];
    sum2 += prwt(lfd,ii)*wght(des,ii)*((lb[ii]>0) ? lb[ii]/p : -lb[ii]/(1-p) );
    sum1 -= prwt(lfd,ii)*wght(des,ii)*((lb[ii]>0) ? lb[ii]/(1-p) : -lb[ii]/p );
    if (sum1<=sum2+tol1)
    {
/* determine the range of ties [i0,i1]
 *   el[ind[i0..i1]] = el[ind[i]].
 *   if sum1==sum2, el[ind[i+1]]..el[ind[i1]]] = el[ind[i1]], else i1 = i.
 */
      i0 = i1 = i;
      while ((i0>0) && (el[des->ind[i0-1]]==el[ii])) i0--;
      while ((i1<des->n-1) && (el[des->ind[i1+1]]==el[ii])) i1++;
      if (sum1>=sum2-tol1)
        while ((i1<des->n-1) && (el[des->ind[i1+1]]==el[des->ind[i+1]])) i1++;

      if (i0<i1) ii = tiecycle(des->ind,i0,i1,oi);
      for (j=0; j<des->p; j++) cf[j] += el[ii]*dc[j];
      return(ii);
    }
  }
mut_printf("Big finddlt problem.\n");
ii = des->ind[des->n-1];
for (j=0; j<des->p; j++) cf[j] += el[ii]*dc[j];
return(ii);
}

/*
 * special version of movecoef for min/max.
 */
int movemin(lfd,des,f,cf,dc,oi)
design *des;
lfdata *lfd;
double *cf, *dc, f;
int oi;
{ int i, ii, im, p, s, ssum;
  double *lb, sum, lb0, lb1, z0, z1;

  lb = des->th;
  s = (f<=0.0) ? 1 : -1;

/* first, determine whether move should be in positive or negative direction */
  p = des->p;
  sum = 0;
  for (i=0; i<des->n; i++)
  { ii = des->ind[i];
    lb[ii] = innerprod(dc,d_xi(des,ii),des->p);
    sum += prwt(lfd,ii)*wght(des,ii)*lb[ii];
  }
  if (fabs(sum) <= QTOL*q0)
  { lf_status = LF_PF;
    return(0);
  }
  ssum = (sum<=0.0) ? -1 : 1;
  if (ssum != s)
    for (i=0; i<p; i++) dc[i] = -dc[i];

/* now, move positively. How far can we move? */
  lb0 = 1.0e100; im = oi;
  for (i=0; i<des->n; i++)
  { ii = des->ind[i];
    lb[ii] = innerprod(dc,d_xi(des,ii),des->p); /* must recompute - signs! */
    if (s*lb[ii]>QTOL) /* should have scale-free tolerance here */
    { z0 = innerprod(cf,d_xi(des,ii),p);
      lb1 = (resp(lfd,ii) - z0)/lb[ii];
      if (lb1<lb0)
      { if (fabs(lb1-lb0)<QTOL) /* cycle */
        { if (im<=oi)
          { if ((ii>oi) | (ii<im)) im = ii; }
          else
          { if ((ii>oi) & (ii<im)) im = ii; }
        }
        else
        { im = ii; lb0 = lb1; }
      }
    }
  }

  for (i=0; i<p; i++) cf[i] = cf[i]+lb0*dc[i];
  if (im==-1) lf_status = LF_PF;
  return(im);
}

double qll(lfd,spr,des,cf)
lfdata *lfd;
smpar *spr;
design *des;
double *cf;
{ int i, ii;
  double th, sp, sn, p, e;

  p = pen(spr);
  sp = sn = 0.0;
  for (i=0; i<des->n; i++)
  { ii = des->ind[i];
    th = innerprod(d_xi(des,ii),cf,des->p);
    e = resp(lfd,ii)-th;
    if (e<0) sn -= prwt(lfd,ii)*wght(des,ii)*e;
    if (e>0) sp += prwt(lfd,ii)*wght(des,ii)*e;
  }
  if (p<=0.0) return((sn<QTOL) ? -sp : -1e300);
  if (p>=1.0) return((sp<QTOL) ? -sn : -1e300);
  return(-sp/(1-p)-sn/p);
}

/*
 * running quantile smoother.
 */
void lfquantile(lfd,sp,des,maxit)
lfdata *lfd;
smpar *sp;
design *des;
int maxit;
{ int i, ii, im, j, k, p, *ci, (*mover)();
  double *cf, *db, *dc, *cm, f, q1, q2, l0;

printf("in lfquantile\n");
  f = pen(sp);
  p = des->p;
  cf = des->cf;
  dc = des->oc;
  db = des->ss;
  setzero(cf,p);
  setzero(dc,p);
  cm = des->V;
  setzero(cm,p*p);
  ci = (int *)des->fix;

  q1 = -qll(lfd,sp,des,cf);
  if (q1==0.0) { lf_status = LF_PF; return; }
  for (i=0; i<p; i++) cm[i*(p+1)] = 1;
  mover = movecoef;
  if ((f<=0.0) | (f>=1.0)) mover = movemin;

  dc[0] = 1.0;
  im = mover(lfd,des,f,cf,dc,-1);
  if (lf_status != LF_OK) return;
  ci[0] = im;
printf("init const %2d\n",ci[0]);
  q0 = -qll(lfd,sp,des,cf);
  if (q0<QTOL*q1) { lf_status = LF_PF; return; }

printf("loop 0\n"); fflush(stdout);
  for (i=1; i<p; i++)
  {
printf("i %2d\n",i);
    memcpy(&cm[(i-1)*p],d_xi(des,im),p*sizeof(double));
    setzero(db,p);
    db[i] = 1.0;
    resproj(db,cm,dc,p,i);
printf("call mover\n"); fflush(stdout);
    im = mover(lfd,des,f,cf,dc,-1);
    if (lf_status != LF_OK) return;
printf("mover %2d\n",im); fflush(stdout);
    ci[i] = im;
  }
printf("call qll\n"); fflush(stdout);
  q1 = qll(lfd,sp,des,cf);

printf("loop 1    %d %d %d %d\n",ci[0],ci[1],ci[2],ci[3]); fflush(stdout);
  for (k=0; k<maxit; k++)
  { for (i=0; i<p; i++)
    { for (j=0; j<p; j++)
        if (j!=i) memcpy(&cm[(j-(j>i))*p],d_xi(des,ci[j]),p*sizeof(double));
      memcpy(db,d_xi(des,ci[i]),p*sizeof(double));
      resproj(db,cm,dc,p,p-1);
printf("call mover\n"); fflush(stdout);
      im = mover(lfd,des,f,cf,dc,ci[i]);
      if (lf_status != LF_OK) return;
printf("mover %2d\n",im); fflush(stdout);
      ci[i] = im;
    }
    q2 = qll(lfd,sp,des,cf);
/*
 * convergence: require no change -- reasonable, since discrete?
 * remember we're maximizing, and q's are negative.
 */
     if (q2 <= q1) return;
     q1 = q2;
  }
printf("loop 2\n");
  mut_printf("Warning: lfquantile not converged.\n");
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
#include "locf.h"

extern double links_rs;

int robust_vallink(link)
int link;
{ return(link==LIDENT);
}

int robust_fam(y,mean,th,link,res,cens,w)
double y, mean, th, *res, w;
int link, cens;
{ double z, sw;
  if (link==LINIT)
  { res[ZDLL] = w*y;
    return(LF_OK);
  }
  sw = (w==1.0) ? 1.0 : sqrt(w); /* don't want unnecess. sqrt! */
  z = sw*(y-mean)/links_rs;
  res[ZLIK] = (fabs(z)<HUBERC) ? -z*z/2 : HUBERC*(HUBERC/2.0-fabs(z));
  if (z< -HUBERC)
  { res[ZDLL] = -sw*HUBERC/links_rs;
    res[ZDDLL]= 0.0;
    return(LF_OK);
  }
  if (z> HUBERC)
  { res[ZDLL] = sw*HUBERC/links_rs;
    res[ZDDLL]= 0.0;
    return(LF_OK);
  }
  res[ZDLL] =  sw*z/links_rs;
  res[ZDDLL] = w/(links_rs*links_rs);
  return(LF_OK);
}

int cauchy_fam(y,p,th,link,res,cens,w)
double y, p, th, *res, w;
int link, cens;
{ double z;
  if (link!=LIDENT)
  { LERR(("Invalid link in famcauc"));
    return(LF_LNK);
  }
  z = w*(y-th)/links_rs;
  res[ZLIK] = -log(1+z*z);
  res[ZDLL] = 2*w*z/(links_rs*(1+z*z));
  res[ZDDLL] = 2*w*w*(1-z*z)/(links_rs*links_rs*(1+z*z)*(1+z*z));
  return(LF_OK);
}

extern double lf_tol;
int robust_init(lfd,des,sp)
lfdata *lfd;
design *des;
smpar *sp;
{ int i;
  for (i=0; i<des->n; i++)
  des->res[i] = resp(lfd,(int)des->ind[i]) - base(lfd,(int)des->ind[i]);
  des->cf[0] = median(des->res,des->n);
  for (i=1; i<des->p; i++) des->cf[i] = 0.0;
  lf_tol = 1.0e-6;
  return(LF_OK);
}

void setfrobust(fam)
family *fam;
{ fam->deflink = LIDENT;
  fam->canlink = LIDENT;
  fam->vallink = robust_vallink;
  fam->family  = robust_fam;
  fam->initial = robust_init;
  fam->robust = 0;
}

void setfcauchy(fam)
family *fam;
{ fam->deflink = LIDENT;
  fam->canlink = LIDENT;
  fam->vallink = robust_vallink;
  fam->family  = cauchy_fam;
  fam->initial = robust_init;
  fam->robust = 0;
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
#include "locf.h"

int weibull_vallink(link)
int link;
{ return((link==LIDENT) | (link==LLOG) | (link==LLOGIT));
}

int weibull_fam(y,mean,th,link,res,cens,w)
double y, mean, th, *res, w;
int link, cens;
{ double yy;
  yy = pow(y,w);
  if (link==LINIT)
  { res[ZDLL] = MAX(yy,0.0);
    return(LF_OK);
  }
  if (cens)
  { res[ZLIK] = -yy/mean;
    res[ZDLL] = res[ZDDLL] = yy/mean;
    return(LF_OK);
  }
  res[ZLIK] = 1-yy/mean-th;
  if (yy>0) res[ZLIK] += log(w*yy);
  res[ZDLL] = -1+yy/mean;
  res[ZDDLL]= yy/mean;
  return(LF_OK);
}

void setfweibull(fam)
family *fam;
{ fam->deflink = LLOG;
  fam->canlink = LLOG;
  fam->vallink = weibull_vallink;
  fam->family  = weibull_fam;
  fam->robust = 0;
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
/*
  Functions implementing the adaptive bandwidth selection.
  Will make the final call to nbhd() to set smoothing weights
  for selected bandwidth, But will **not** make the
  final call to locfit().
*/

#include "locf.h"

static double hmin;

#define NACRI 5
static char *atype[NACRI] = { "none", "cp", "ici", "mindex", "ok" };
static int   avals[NACRI] = { ANONE, ACP, AKAT, AMDI, AOK };
int lfacri(char *z)
{ return(pmatch(z, atype, avals, NACRI, ANONE));
}

double adcri(lk,t0,t2,pen)
double lk, t0, t2, pen;
{ double y;
/* return(-2*lk/(t0*exp(pen*log(1-t2/t0)))); */
  /* return((-2*lk+pen*t2)/t0); */
  y = (MAX(-2*lk,t0-t2)+pen*t2)/t0;
  return(y);
}

double mmse(lfd,sp,dv,des)
lfdata *lfd;
smpar *sp;
deriv *dv;
design *des;
{ int i, ii, j, p, p1;
  double sv, sb, *l, dp;

  l = des->wd;
  wdiag(lfd, sp, des,l,dv,0,1,0);
  sv = sb = 0;
  p = npar(sp);
  for (i=0; i<des->n; i++)
  { sv += l[i]*l[i];
    ii = des->ind[i];
    dp = dist(des,ii);
    for (j=0; j<deg(sp); j++) dp *= dist(des,ii);
    sb += fabs(l[i])*dp;
  }
  p1 = factorial(deg(sp)+1);
printf("%8.5f sv %8.5f  sb %8.5f  %8.5f\n",des->h,sv,sb,sv+sb*sb*pen(sp)*pen(sp)/(p1*p1));
  return(sv+sb*sb*pen(sp)*pen(sp)/(p1*p1));
}

static double mcp, clo, cup;

/*
  Initial bandwidth will be (by default)
  k-nearest neighbors for k small, just large enough to
  get defined estimate (unless user provided nonzero nn or fix-h components)
*/

int ainitband(lfd,sp,dv,des)
lfdata *lfd;
smpar *sp;
deriv *dv;
design *des;
{ int lf_status, p, z, cri, noit, redo;
  double ho, t[6];

  if (lf_debug >= 2) mut_printf("ainitband:\n");
  p = des->p;
  cri = acri(sp);
  noit = (cri!=AOK);
  z = (int)(lfd->n*nn(sp));
  if ((noit) && (z<p+2)) z = p+2;
  redo = 0; ho = -1;
  do
  { 
    nbhd(lfd,des,z,redo,sp);
    if (z<des->n) z = des->n;
    if (des->h>ho) lf_status = locfit(lfd,des,sp,noit,0,0);
    z++;
    redo = 1;
  } while ((z<=lfd->n) && ((des->h==0)||(lf_status!=LF_OK)));
  hmin = des->h;

  switch(cri)
  { case ACP:
      local_df(lfd,sp,des,t);
      mcp = adcri(des->llk,t[0],t[2],pen(sp));
      return(lf_status);
    case AKAT:
      local_df(lfd,sp,des,t);
      clo = des->cf[0]-pen(sp)*t[5];
      cup = des->cf[0]+pen(sp)*t[5];
      return(lf_status);
    case AMDI:
      mcp = mmse(lfd,sp,dv,des);
      return(lf_status);
    case AOK: return(lf_status);
  }
  LERR(("aband1: unknown criterion"));
  return(LF_ERR);
}

/*
  aband2 increases the initial bandwidth until lack of fit results,
  or the fit is close to a global fit. Increase h by 1+0.3/d at
  each iteration.
*/

double aband2(lfd,sp,dv,des,h0)
lfdata *lfd;
smpar *sp;
deriv *dv;
design *des;
double h0;
{ double t[6], h1, nu1, cp, ncp, tlo, tup;
  int d, inc, n, p, done;

  if (lf_debug >= 2) mut_printf("aband2:\n");
  d = lfd->d; n = lfd->n; p = npar(sp);
  h1 = des->h = h0;
  done = 0; nu1 = 0.0;
  inc = 0; ncp = 0.0;
  while ((!done) & (nu1<(n-p)*0.95))
  { fixh(sp) = (1+0.3/d)*des->h;
    nbhd(lfd,des,0,1,sp);
    if (locfit(lfd,des,sp,1,0,0) > 0) WARN(("aband2: failed fit"));
    local_df(lfd,sp,des,t);
    nu1 = t[0]-t[2]; /* tr(A) */
    switch(acri(sp))
    { case AKAT:
        tlo = des->cf[0]-pen(sp)*t[5];
        tup = des->cf[0]+pen(sp)*t[5];
/* mut_printf("h %8.5f  tlo %8.5f  tup %8.5f\n",des->h,tlo,tup); */
        done = ((tlo>cup) | (tup<clo));
        if (!done)
        { clo = MAX(clo,tlo);
          cup = MIN(cup,tup);
          h1 = des->h;
        }
        break;
      case ACP:
        cp = adcri(des->llk,t[0],t[2],pen(sp));
/* mut_printf("h %8.5f  lk %8.5f  t0 %8.5f  t2 %8.5f  cp %8.5f\n",des->h,des->llk,t[0],t[2],cp); */
        if (cp<mcp) { mcp = cp; h1 = des->h; }
        if (cp>=ncp) inc++; else inc = 0;
        ncp = cp;
        done = (inc>=10) | ((inc>=3) & ((t[0]-t[2])>=10) & (cp>1.5*mcp));
        break;
      case AMDI:
        cp = mmse(lfd,sp,dv,des);
        if (cp<mcp) { mcp = cp; h1 = des->h; }
        if (cp>ncp) inc++; else inc = 0;
        ncp = cp;
        done = (inc>=3);
        break;
    }
  }
  return(h1);
}

/*
  aband3 does a finer search around best h so far. Try
  h*(1-0.2/d), h/(1-0.1/d), h*(1+0.1/d), h*(1+0.2/d)
*/
double aband3(lfd,sp,dv,des,h0)
lfdata *lfd;
smpar *sp;
deriv *dv;
design *des;
double h0;
{ double t[6], h1, cp, tlo, tup;
  int i, i0, d, n;

  if (lf_debug >= 2) mut_printf("aband3:\n");
  d = lfd->d; n = lfd->n;
  h1 = h0;
  i0 = (acri(sp)==AKAT) ? 1 : -2;
  if (h0==hmin) i0 = 1;

  for (i=i0; i<=2; i++)
  { if (i==0) i++;
    fixh(sp) = h0*(1+0.1*i/d);
    nbhd(lfd,des,0,1,sp);
    if (locfit(lfd,des,sp,1,0,0) > 0) WARN(("aband3: failed fit"));
    local_df(lfd,sp,des,t);
    switch (acri(sp))
    { case AKAT:
        tlo = des->cf[0]-pen(sp)*t[5];
        tup = des->cf[0]+pen(sp)*t[5];
        if ((tlo>cup) | (tup<clo)) /* done */
          i = 2;
        else
        { h1 = des->h;
          clo = MAX(clo,tlo);
          cup = MIN(cup,tup);
        }
        break;
      case ACP:
        cp = adcri(des->llk,t[0],t[2],pen(sp));
        if (cp<mcp) { mcp = cp; h1 = des->h; }
        else
        { if (i>0) i = 2; }
        break;
      case AMDI:
        cp = mmse(lfd,sp,dv,des);
        if (cp<mcp) { mcp = cp; h1 = des->h; }
        else
        { if (i>0) i = 2; }
    }
  }
  return(h1);
}

int alocfit(lfd,sp,dv,des,cv)
lfdata *lfd;
smpar *sp;
deriv *dv;
design *des;
int cv;
{ int lf_status;
  double h0;

  lf_status = ainitband(lfd,sp,dv,des);
  if (lf_error) return(lf_status);
  if (acri(sp) == AOK) return(lf_status);

  h0 = fixh(sp);
  fixh(sp) = aband2(lfd,sp,dv,des,des->h);
  fixh(sp) = aband3(lfd,sp,dv,des,fixh(sp));
  nbhd(lfd,des,0,1,sp);
  lf_status = locfit(lfd,des,sp,0,0,cv);
  fixh(sp) = h0;

  return(lf_status);
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
/*
 *
 *   Evaluate the locfit fitting functions.
 *     calcp(sp,d)
 *       calculates the number of fitting functions.
 *     makecfn(sp,des,dv,d)
 *       makes the coef.number vector.
 *     fitfun(lfd, sp, x,t,f,dv)
 *       lfd is the local fit structure.
 *       sp  smoothing parameter structure.
 *       x is the data point.
 *       t is the fitting point.
 *       f is a vector to return the results.
 *       dv derivative structure.
 *     designmatrix(lfd, sp, des)
 *       is a wrapper for fitfun to build the design matrix.
 *
 */

#include "locf.h"

int calcp(sp,d)
smpar *sp;
int d;
{ int i, k;

  if (ubas(sp)) return(npar(sp));

  switch (kt(sp))
  { case KSPH:
    case KCE:
      k = 1;
      for (i=1; i<=deg(sp); i++) k = k*(d+i)/i;
      return(k);
    case KPROD: return(d*deg(sp)+1);
    case KLM: return(d);
    case KZEON: return(1);
  }
  LERR(("calcp: invalid kt %d",kt(sp)));
  return(0);
}

int coefnumber(dv,kt,d,deg)
int kt, d, deg;
deriv *dv;
{ int d0, d1, t;

  if (d==1)
  { if (dv->nd<=deg) return(dv->nd);
    return(-1);
  }

  if (dv->nd==0) return(0);
  if (deg==0) return(-1);
  if (dv->nd==1) return(1+dv->deriv[0]);
  if (deg==1) return(-1);
  if (kt==KPROD) return(-1);

  if (dv->nd==2)
  { d0 = dv->deriv[0]; d1 = dv->deriv[1];
    if (d0<d1) { t = d0; d0 = d1; d1 = t; }
    return((d+1)*(d0+1)-d0*(d0+3)/2+d1);
  }
  if (deg==2) return(-1);

  LERR(("coefnumber not programmed for nd>=3"));
  return(-1);
}

void makecfn(sp,des,dv,d)
smpar *sp;
design *des;
deriv *dv;
int d;
{ int i, nd;
  
  nd = dv->nd;

  des->cfn[0] = coefnumber(dv,kt(sp),d,deg(sp));
  des->ncoef = 1;
  if (nd >= deg(sp)) return;
  if (kt(sp)==KZEON) return;

  if (d>1)
  { if (nd>=2) return;
    if ((nd>=1) && (kt(sp)==KPROD)) return;
  }

  dv->nd = nd+1;
  for (i=0; i<d; i++)
  { dv->deriv[nd] = i;
    des->cfn[i+1] = coefnumber(dv,kt(sp),d,deg(sp));
  }
  dv->nd = nd;

  des->ncoef = 1+d;
}

void fitfunangl(dx,ff,sca,cd,deg)
double dx, *ff, sca;
int deg, cd;
{
  if (deg>=3) WARN(("Can't handle angular model with deg>=3"));

  switch(cd)
  { case 0:
      ff[0] = 1;
      ff[1] = sin(dx/sca)*sca;
      ff[2] = (1-cos(dx/sca))*sca*sca;
      return;
    case 1:
      ff[0] = 0;
      ff[1] = cos(dx/sca);
      ff[2] = sin(dx/sca)*sca;
      return;
    case 2:
      ff[0] = 0;
      ff[1] = -sin(dx/sca)/sca;
      ff[2] = cos(dx/sca);
      return;
    default: WARN(("Can't handle angular model with >2 derivs"));
  }
}

void fitfun(lfd,sp,x,t,f,dv)
lfdata *lfd;
smpar *sp;
double *x, *t, *f;
deriv *dv;
{ int d, deg, nd, m, i, j, k, ct_deriv[MXDIM];
  double ff[MXDIM][1+MXDEG], dx[MXDIM], *xx[MXDIM];

  if (ubas(sp))
  { for (i=0; i<lfd->d; i++) xx[i] = &x[i];
    i = 0;
    sp->vbasis(xx,t,1,lfd->d,1,npar(sp),f);
    return;
  }

  d = lfd->d;
  deg = deg(sp);
  m = 0;
  nd = (dv==NULL) ? 0 : dv->nd;

  if (kt(sp)==KZEON)
  { f[0] = 1.0;
    return;
  }

  if (kt(sp)==KLM)
  { for (i=0; i<d; i++) f[m++] = x[i];
    return;
  }

  f[m++] = (nd==0);
  if (deg==0) return;

  for (i=0; i<d; i++)
  { ct_deriv[i] = 0;
    dx[i] = (t==NULL) ? x[i] : x[i]-t[i];
  }
  for (i=0; i<nd; i++) ct_deriv[dv->deriv[i]]++;

  for (i=0; i<d; i++)
  { switch(lfd->sty[i])
    {
      case STANGL:
        fitfunangl(dx[i],ff[i],lfd->sca[i],ct_deriv[i],deg(sp));
        break;
      default:
        for (j=0; j<ct_deriv[i]; j++) ff[i][j] = 0.0;
        ff[i][ct_deriv[i]] = 1.0;
        for (j=ct_deriv[i]+1; j<=deg; j++)
          ff[i][j] = ff[i][j-1]*dx[i]/(j-ct_deriv[i]);
    }
  }

/*
 *  Product kernels. Note that if ct_deriv[i] != nd, that implies
 *  there is differentiation wrt another variable, and all components
 *  involving x[i] are 0.
 */
  if ((d==1) || (kt(sp)==KPROD))
  { for (j=1; j<=deg; j++)
      for (i=0; i<d; i++)
        f[m++] = (ct_deriv[i]==nd) ? ff[i][j] : 0.0;
    return;
  }

/*
 *  Spherical kernels with the full polynomial basis.
 *  Presently implemented up to deg=3.
 */
  for (i=0; i<d; i++)
    f[m++] = (ct_deriv[i]==nd) ? ff[i][1] : 0.0;
  if (deg==1) return;

  for (i=0; i<d; i++)
  {
    /* xi^2/2 terms. */
    f[m++] = (ct_deriv[i]==nd) ? ff[i][2] : 0.0;

    /* xi xj terms */
    for (j=i+1; j<d; j++)
      f[m++] = (ct_deriv[i]+ct_deriv[j]==nd) ? ff[i][1]*ff[j][1] : 0.0;
  }
  if (deg==2) return;

  for (i=0; i<d; i++)
  { 
    /* xi^3/6 terms */
    f[m++] = (ct_deriv[i]==nd) ? ff[i][3] : 0.0;

    /* xi^2/2 xk terms */
    for (k=i+1; k<d; k++)
      f[m++] = (ct_deriv[i]+ct_deriv[k]==nd) ? ff[i][2]*ff[k][1] : 0.0;

    /* xi xj xk terms */
    for (j=i+1; j<d; j++)
    { f[m++] = (ct_deriv[i]+ct_deriv[j]==nd) ? ff[i][1]*ff[j][2] : 0.0;
      for (k=j+1; k<d; k++)
        f[m++] = (ct_deriv[i]+ct_deriv[j]+ct_deriv[k]==nd) ?
                    ff[i][1]*ff[j][1]*ff[k][1] : 0.0;
    }
  }
  if (deg==3) return;

  LERR(("fitfun: can't handle deg=%d for spherical kernels",deg));
}

/*
 *  Build the design matrix. Assumes des->ind contains the indices of
 *  the required data points; des->n the number of points; des->xev
 *  the fitting point.
 */
void designmatrix(lfd,sp,des)
lfdata *lfd;
smpar *sp;
design *des;
{ int i, ii, j, p;
  double *X, u[MXDIM];

  X = d_x(des);
  p = des->p;

  if (ubas(sp))
  {
    sp->vbasis(lfd->x,des->xev,lfd->n,lfd->d,des->n,p,X);
    return;
  }

  for (i=0; i<des->n; i++)
  { ii = des->ind[i];
    for (j=0; j<lfd->d; j++) u[j] = datum(lfd,j,ii);
    fitfun(lfd,sp,u,des->xev,&X[ii*p],NULL);
  }
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
/*
 *
 *
 *  Functions for determining bandwidth; smoothing neighborhood
 *  and smoothing weights.
 */

#include "locf.h"

double rho(x,sc,d,kt,sty) /* ||x|| for appropriate distance metric */
double *x, *sc;
int d, kt, *sty;
{ double rhoi[MXDIM], s;
  int i;
  for (i=0; i<d; i++)
  { if (sty!=NULL)
    { switch(sty[i])
      { case STANGL:  rhoi[i] = 2*sin(x[i]/(2*sc[i])); break;
        case STCPAR: rhoi[i] = 0; break;
        default: rhoi[i] = x[i]/sc[i];
    } }
    else rhoi[i] = x[i]/sc[i];
  }

  if (d==1) return(fabs(rhoi[0]));

  s = 0;
  if (kt==KPROD)
  { for (i=0; i<d; i++)
    { rhoi[i] = fabs(rhoi[i]);
      if (rhoi[i]>s) s = rhoi[i];
    }
    return(s);
  }

  if (kt==KSPH)
  { for (i=0; i<d; i++)
      s += rhoi[i]*rhoi[i];
    return(sqrt(s));
  }

  LERR(("rho: invalid kt"));
  return(0.0);
}

double kordstat(x,k,n,ind)
double *x;
int k, n, *ind;
{ int i, i0, i1, l, r;
  double piv;
  if (k<1) return(0.0);
  i0 = 0; i1 = n-1;
  while (1)
  { piv = x[ind[(i0+i1)/2]];
    l = i0; r = i1;
    while (l<=r)
    { while ((l<=i1) && (x[ind[l]]<=piv)) l++;
      while ((r>=i0) && (x[ind[r]]>piv)) r--;
      if (l<=r) ISWAP(ind[l],ind[r]);
    } /* now, x[ind[i0..r]] <= piv < x[ind[l..i1]] */
    if (r<k-1) i0 = l;  /* go right */
    else /* put pivots in middle */
    { for (i=i0; i<=r; )
        if (x[ind[i]]==piv) { ISWAP(ind[i],ind[r]); r--; }
        else i++;
      if (r<k-1) return(piv);
      i1 = r;
    }
  }
}

/* check if i'th data point is in limits */
int inlim(lfd,i)
lfdata *lfd;
int i;
{ int d, j, k;
  double *xlim;

  xlim = lfd->xl;
  d = lfd->d;
  k = 1;
  for (j=0; j<d; j++)
  { if (xlim[j]<xlim[j+d])
      k &= ((datum(lfd,j,i)>=xlim[j]) & (datum(lfd,j,i)<=xlim[j+d]));
  }
  return(k);
}

double compbandwid(di,ind,x,n,d,nn,fxh)
double *di, *x, fxh;
int n, d, nn, *ind;
{ int i;
  double nnh;

  if (nn==0) return(fxh);

  if (nn<n)
    nnh = kordstat(di,nn,n,ind);
  else
  { nnh = 0;
    for (i=0; i<n; i++) nnh = MAX(nnh,di[i]);
    nnh = nnh*exp(log(1.0*nn/n)/d);
  }
  return(MAX(fxh,nnh));
}

/*
  fast version of nbhd for ordered 1-d data
*/
void nbhd1(lfd,sp,des,k)
lfdata *lfd;
smpar *sp;
design *des;
int k;
{ double x, h, *xd, sc;
  int i, l, r, m, n, z;

  n = lfd->n;
  x = des->xev[0];
  xd = dvari(lfd,0);
  sc = lfd->sca[0];

  /* find closest data point to x */
  if (x<=xd[0]) z = 0;
  else
  if (x>=xd[n-1]) z = n-1;
  else
  { l = 0; r = n-1;
    while (r-l>1)
    { z = (r+l)/2;
      if (xd[z]>x) r = z;
              else l = z;
    }
    /* now, xd[0..l] <= x < x[r..n-1] */
    if ((x-xd[l])>(xd[r]-x)) z = r; else z = l;
  }
  /* closest point to x is xd[z] */

  if (nn(sp)<0)  /* user bandwidth */
    h = sp->vb(des->xev);
  else
  { if (k>0) /* set h to nearest neighbor bandwidth */
    { l = r = z;
      if (l==0) r = k-1;
      if (r==n-1) l = n-k;
      while (r-l<k-1)
      { if ((x-xd[l-1])<(xd[r+1]-x)) l--; else r++;
        if (l==0) r = k-1;
        if (r==n-1) l = n-k;
      }
      h = x-xd[l];
      if (h<xd[r]-x) h = xd[r]-x;
    }
    else h = 0;
    h /= sc;
    if (h<fixh(sp)) h = fixh(sp);
  }

  m = 0;
  if (xd[z]>x) z--; /* so xd[z]<=x */
  /* look left */
  for (i=z; i>=0; i--) if (inlim(lfd,i))
  { dist(des,i) = (x-xd[i])/sc;
    wght(des,i) = weight(lfd, sp, &xd[i], &x, h, 1, dist(des,i));
    if (wght(des,i)>0)
    { des->ind[m] = i;
      m++; 
    } else i = 0;
  }
  /* look right */
  for (i=z+1; i<n; i++) if (inlim(lfd,i))
  { dist(des,i) = (xd[i]-x)/sc;
    wght(des,i) = weight(lfd, sp, &xd[i], &x, h, 1, dist(des,i));
    if (wght(des,i)>0)
    { des->ind[m] = i;
      m++; 
    } else i = n;
  }

  des->n = m;
  des->h = h;
}

void nbhd_zeon(lfd,des)
lfdata *lfd;
design *des;
{ int i, j, m, eq;

  m = 0;
  for (i=0; i<lfd->n; i++)
  { eq = 1;
    for (j=0; j<lfd->d; j++) eq = eq && (des->xev[j] == datum(lfd,j,i));
    if (eq)
    { wght(des,i) = 1;
      des->ind[m] = i;
      m++;
    }
  }
  des->n = m;
  des->h = 1.0;
}

void nbhd(lfd,des,nn,redo,sp)
lfdata *lfd;
design *des;
int redo, nn;
smpar *sp;
{ int d, i, j, m, n;
  double h, u[MXDIM];

  if (lf_debug>1) mut_printf("nbhd: nn %d  fixh %8.5f\n",nn,fixh(sp));
  
  d = lfd->d; n = lfd->n;

  if (ker(sp)==WPARM)
  { for (i=0; i<n; i++)
    { wght(des,i) = 1.0;
      des->ind[i] = i;
    }
    des->n = n;
    return;
  }

  if (kt(sp)==KZEON)
  { nbhd_zeon(lfd,des);
    return;
  }

  if (kt(sp)==KCE)
  { des->h = 0.0;
    return;
  }

  /* ordered 1-dim; use fast searches */
  if ((nn<=n) & (lfd->ord) & (ker(sp)!=WMINM) & (lfd->sty[0]!=STANGL))
  { nbhd1(lfd,sp,des,nn);
    return;
  }

  if (!redo)
  { for (i=0; i<n; i++)
    { for (j=0; j<d; j++) u[j] = datum(lfd,j,i)-des->xev[j];
      dist(des,i) = rho(u,lfd->sca,d,kt(sp),lfd->sty);
      des->ind[i] = i;
    }
  }
  else
    for (i=0; i<n; i++) des->ind[i] = i;

  if (ker(sp)==WMINM)
  { des->h = minmax(lfd,des,sp);
    return;
  }

  if (nn<0)
    h = sp->vb(des->xev);
  else
    h = compbandwid(des->di,des->ind,des->xev,n,lfd->d,nn,fixh(sp));
  m = 0;
  for (i=0; i<n; i++) if (inlim(lfd,i))
  { for (j=0; j<d; j++) u[j] = datum(lfd,j,i);
    wght(des,i) = weight(lfd, sp, u, des->xev, h, 1, dist(des,i));
    if (wght(des,i)>0)
    { des->ind[m] = i;
      m++;
    }
  }
  des->n = m;
  des->h = h;
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
/*
 *
 *   This file includes functions to solve for the scale estimate in
 *   local robust regression and likelihood. The main entry point is
 *   lf_robust(lfd,sp,des,mxit),
 *   called from the locfit() function.
 *
 *   The update_rs(x) accepts a residual scale x as the argument (actually,
 *   it works on the log-scale). The function computes the local fit
 *   assuming this residual scale, and re-estimates the scale from this
 *   new fit. The final solution satisfies the fixed point equation
 *   update_rs(x)=x. The function lf_robust() automatically calls
 *   update_rs() through the fixed point iterations.
 *
 *   The estimation of the scale from the fit is based on the sqrt of
 *   the median deviance of observations with non-zero weights (in the
 *   gaussian case, this is the median absolute residual).
 *
 *   TODO:
 *     Should use smoothing weights in the median.
 */

#include "locf.h"

extern int lf_status;
double robscale;

static lfdata *rob_lfd;
static smpar *rob_sp;
static design *rob_des;
static int rob_mxit;

double median(x,n)
double *x;
int n;
{ int i, j, lt, eq, gt;
  double lo, hi, s;
  lo = hi = x[0];
  for (i=0; i<n; i++)
  { lo = MIN(lo,x[i]);
    hi = MAX(hi,x[i]);
  }
  if (lo==hi) return(lo);
  lo -= (hi-lo);
  hi += (hi-lo);
  for (i=0; i<n; i++)
  { if ((x[i]>lo) & (x[i]<hi))
    { s = x[i]; lt = eq = gt = 0;
      for (j=0; j<n; j++)
      { lt += (x[j]<s);
        eq += (x[j]==s);
        gt += (x[j]>s);
      }
      if ((2*(lt+eq)>n) && (2*(gt+eq)>n)) return(s);
      if (2*(lt+eq)<=n) lo = s;
      if (2*(gt+eq)<=n) hi = s;
    }
  }
  return((hi+lo)/2);
}

double nrobustscale(lfd,sp,des,rs)
lfdata *lfd;
smpar *sp;
design *des;
double rs;
{ int i, ii, p;
  double link[LLEN], sc, sd, sw, e;
  p = des->p; sc = sd = sw = 0.0;
  for (i=0; i<des->n; i++)
  { ii = des->ind[i];
    fitv(des,ii) = base(lfd,ii)+innerprod(des->cf,d_xi(des,ii),p);
    e = resp(lfd,ii)-fitv(des,ii);
    stdlinks(link,lfd,sp,ii,fitv(des,ii),rs);
    sc += wght(des,ii)*e*link[ZDLL];
    sd += wght(des,ii)*e*e*link[ZDDLL];
    sw += wght(des,ii);
  }

  /* newton-raphson iteration for log(s)
     -psi(ei/s) - log(s); s = e^{-th}
  */
  rs *= exp((sc-sw)/(sd+sc));
  return(rs);
}

double robustscale(lfd,sp,des)
lfdata *lfd;
smpar *sp;
design *des;
{ int i, ii, p, fam, lin, or;
  double rs, link[LLEN];
  p = des->p;
  fam = fam(sp);
  lin = link(sp);
  or = fami(sp)->robust;
  fami(sp)->robust = 0;

  for (i=0; i<des->n; i++)
  { ii = des->ind[i];
    fitv(des,ii) = base(lfd,ii) + innerprod(des->cf,d_xi(des,ii),p);
    links(fitv(des,ii),resp(lfd,ii),fami(sp),lin,link,cens(lfd,ii),prwt(lfd,ii),1.0);
    des->res[i] = -2*link[ZLIK];
  }
  fami(sp)->robust = or;
  rs = sqrt(median(des->res,des->n));

  if (rs==0.0) rs = 1.0;
  return(rs);
}

double update_rs(x)
double x;
{ double nx;
  if (lf_status != LF_OK) return(x);
  robscale = exp(x);
  lfiter(rob_lfd,rob_sp,rob_des,rob_mxit);
  if (lf_status != LF_OK) return(x);

  nx = log(robustscale(rob_lfd,rob_sp,rob_des));
  if (nx<x-0.2) nx = x-0.2;
  return(nx);
}

void lf_robust(lfd,sp,des,mxit)
lfdata *lfd;
design *des;
smpar *sp;
int mxit;
{ double x;
  rob_lfd = lfd;
  rob_des = des;
  rob_sp = sp;
  rob_mxit = mxit;
  lf_status = LF_OK;

  x = log(robustscale(lfd,sp,des));

  solve_fp(update_rs, x, 1.0e-6, mxit);
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
/*
 *   Post-fitting functions to compute the local variance and
 *   influence functions. Also the local degrees of freedom
 *   calculations for adaptive smoothing.
 */

#include "locf.h"

extern double robscale;

/*
  vmat() computes (after the local fit..) the matrix 
  M2  = X^T W^2 V X.
  M12 = (X^T W V X)^{-1} M2
  Also, for convenience, tr[0] = sum(wi) tr[1] = sum(wi^2).
*/
void vmat(lfd, sp, des, M12, M2)
lfdata *lfd;
smpar *sp;
design *des;
double *M12, *M2;
{ int i, ii, p, nk, ok;
  double link[LLEN], h, ww, tr0, tr1;
  p = des->p;
  setzero(M2,p*p);

  nk = -1;

  /* for density estimation, use integral rather than
     sum form, if W^2 is programmed...
  */
  if ((fam(sp)<=THAZ) && (link(sp)==LLOG))
  { switch(ker(sp))
    { case WGAUS: nk = WGAUS; h = des->h/SQRT2; break;
      case WRECT: nk = WRECT; h = des->h; break;
      case WEPAN: nk = WBISQ; h = des->h; break;
      case WBISQ: nk = WQUQU; h = des->h; break;
      case WTCUB: nk = W6CUB; h = des->h; break;
      case WEXPL: nk = WEXPL; h = des->h/2; break;
    }
  }

  tr0 = tr1 = 0.0;
  if (nk != -1)
  { ok = ker(sp); ker(sp) = nk;
/* compute M2 using integration. Use M12 as work matrix. */
    (des->itype)(des->xev, M2, M12, des->cf, h);
    ker(sp) = ok;
    if (fam(sp)==TDEN) multmatscal(M2,des->smwt,p*p);
    tr0 = des->ss[0];
    tr1 = M2[0]; /* n int W e^<a,A> */
  }
  else
  { for (i=0; i<des->n; i++)
    { ii = des->ind[i];
      stdlinks(link,lfd,sp,ii,fitv(des,ii),robscale);
      ww = SQR(wght(des,ii))*link[ZDDLL];
      tr0 += wght(des,ii);
      tr1 += SQR(wght(des,ii));
      addouter(M2,d_xi(des,ii),d_xi(des,ii),p,ww);
    }
  }
  des->tr0 = tr0;
  des->tr1 = tr1;

  memcpy(M12,M2,p*p*sizeof(double));
  for (i=0; i<p; i++)
    jacob_solve(&des->xtwx,&M12[i*p]);
}

void lf_vcov(lfd,sp,des)
lfdata *lfd;
smpar *sp;
design *des;
{ int i, j, k, p;
  double *M12, *M2;
  M12 = des->V; M2 = des->P; p = des->p;
  vmat(lfd,sp,des,M12,M2); /* M2 = X^T W^2 V X  tr0=sum(W) tr1=sum(W*W) */
  des->tr2 = m_trace(M12,p);   /* tr (XTWVX)^{-1}(XTW^2VX) */

/*
 * Covariance matrix is M1^{-1} * M2 * M1^{-1}
 * We compute this using the cholesky decomposition of
 * M2; premultiplying by M1^{-1} and squaring. This
 * is more stable than direct computation in near-singular cases.
 */
  chol_dec(M2,p,p);
  for (i=0; i<p; i++)
    for (j=0; j<i; j++)
    { M2[j*p+i] = M2[i*p+j];
      M2[i*p+j] = 0.0;
    }
  for (i=0; i<p; i++) jacob_solve(&des->xtwx,&M2[i*p]);
  for (i=0; i<p; i++)
  { for (j=0; j<p; j++)
    { M12[i*p+j] = 0;
      for (k=0; k<p; k++)
        M12[i*p+j] += M2[k*p+i]*M2[k*p+j]; /* ith column of covariance */
    }
  }
  if ((fam(sp)==TDEN) && (link(sp)==LIDENT))
    multmatscal(M12,1/SQR(des->smwt),p*p);

/* this computes the influence function as des->f1[0]. */
  unitvec(des->f1,0,des->p);
  jacob_solve(&des->xtwx,des->f1);
}

/* local_df computes:
 *   tr[0] = trace(W)
 *   tr[1] = trace(W*W)
 *   tr[2] = trace( M1^{-1} M2 )
 *   tr[3] = trace( M1^{-1} M3 )
 *   tr[4] = trace( (M1^{-1} M2)^2 )
 *   tr[5] = var(theta-hat).
 */
void local_df(lfd,sp,des,tr)
lfdata *lfd;
smpar *sp;
design *des;
double *tr;
{ int i, ii, j, p;
  double *m2, *V, ww, link[LLEN];

  tr[0] = tr[1] = tr[2] = tr[3] = tr[4] = tr[5] = 0.0;
  m2 = des->V; V = des->P; p = des->p;

  vmat(lfd,sp,des,m2,V);  /* M = X^T W^2 V X  tr0=sum(W) tr1=sum(W*W) */
  tr[0] = des->tr0;
  tr[1] = des->tr1;
  tr[2] = m_trace(m2,p);   /* tr (XTWVX)^{-1}(XTW^2VX) */

  unitvec(des->f1,0,p);
  jacob_solve(&des->xtwx,des->f1);
  for (i=0; i<p; i++)
    for (j=0; j<p; j++)
    { tr[4] += m2[i*p+j]*m2[j*p+i];  /* tr(M^2) */
      tr[5] += des->f1[i]*V[i*p+j]*des->f1[j]; /* var(thetahat) */
  }
  tr[5] = sqrt(tr[5]);

  setzero(m2,p*p);
  for (i=0; i<des->n; i++)
  { ii = des->ind[i];
    stdlinks(link,lfd,sp,ii,fitv(des,ii),robscale);
    ww = wght(des,ii)*wght(des,ii)*wght(des,ii)*link[ZDDLL];
    addouter(m2,d_xi(des,ii),d_xi(des,ii),p,ww);
  }
  for (i=0; i<p; i++)
  { jacob_solve(&des->xtwx,&m2[i*p]);
    tr[3] += m2[i*(p+1)];
  }

  return;
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
/*
 *  Routines for computing weight diagrams.
 *     wdiag(lf,des,lx,deg,ty,exp)
 *  Must locfit() first, unless ker==WPARM and has par. comp.
 *  
 */

#include "locf.h"

static double *wd;
extern double robscale;
void nnresproj(lfd,sp,des,u,m,p)
lfdata *lfd;
smpar *sp;
design *des;
double *u;
int m, p;
{ int i, ii, j;
  double link[LLEN];
  setzero(des->f1,p);
  for (j=0; j<m; j++)
  { ii = des->ind[j];
    stdlinks(link,lfd,sp,ii,fitv(des,ii),robscale);
    for (i=0; i<p; i++) des->f1[i] += link[ZDDLL]*d_xij(des,j,ii)*u[j];
  }
  jacob_solve(&des->xtwx,des->f1);
  for (i=0; i<m; i++)
  { ii = des->ind[i];
    u[i] -= innerprod(des->f1,d_xi(des,ii),p)*wght(des,ii);
  }
}

void wdexpand(l,n,ind,m)
double *l;
int *ind, n, m;
{ int i, j, t;
  double z;
  for (j=m; j<n; j++) { l[j] = 0.0; ind[j] = -1; }
  j = m-1;
  while (j>=0)
  { if (ind[j]==j) j--;
    else
    { i = ind[j];
      z = l[j]; l[j] = l[i]; l[i] = z;
      t = ind[j]; ind[j] = ind[i]; ind[i] = t;
      if (ind[j]==-1) j--;
    }
  }

/*  for (i=n-1; i>=0; i--)
  { l[i] = ((j>=0) && (ind[j]==i)) ? l[j--] : 0.0; } */
}

int wdiagp(lfd,sp,des,lx,pc,dv,deg,ty,exp)
lfdata *lfd;
smpar *sp;
design *des;
paramcomp *pc;
deriv *dv;
double *lx;
int deg, ty, exp;
{ int i, j, p, nd;
  double *l1;

  p = des->p;

  fitfun(lfd,sp,des->xev,pc->xbar,des->f1,dv);
  if (exp)
  { jacob_solve(&pc->xtwx,des->f1);
    for (i=0; i<lfd->n; i++)
      lx[i] = innerprod(des->f1,d_xi(des,des->ind[i]),p);
    return(lfd->n);
  }
  jacob_hsolve(&pc->xtwx,des->f1);
  for (i=0; i<p; i++) lx[i] = des->f1[i];

  nd = dv->nd;
  dv->nd = nd+1;
  if (deg>=1)
    for (i=0; i<lfd->d; i++)
    { dv->deriv[nd] = i;
      l1 = &lx[(i+1)*p];
      fitfun(lfd,sp,des->xev,pc->xbar,l1,dv);
      jacob_hsolve(&pc->xtwx,l1);
    }

  dv->nd = nd+2;
  if (deg>=2)
    for (i=0; i<lfd->d; i++)
    { dv->deriv[nd] = i;
      for (j=0; j<lfd->d; j++)
      { dv->deriv[nd+1] = j;
        l1 = &lx[(i*lfd->d+j+lfd->d+1)*p];
        fitfun(lfd,sp,des->xev,pc->xbar,l1,dv);
        jacob_hsolve(&pc->xtwx,l1);
    } }
  dv->nd = nd;
  return(p);
}

int wdiag(lfd,sp,des,lx,dv,deg,ty,exp)
lfdata *lfd;
smpar *sp;
design *des;
deriv *dv;
double *lx;
int deg, ty, exp;
/* deg=0: l(x) only.
   deg=1: l(x), l'(x)
   deg=2: l(x), l'(x), l''(x)
   ty = 1: e1 (X^T WVX)^{-1} X^T W        -- hat matrix
   ty = 2: e1 (X^T WVX)^{-1} X^T WV^{1/2} -- scb's
*/
{ double w, *X, *lxd, *lxdd, wdd, wdw, *ulx, link[LLEN], h;
  double dfx[MXDIM], hs[MXDIM];
  int i, ii, j, k, l, m, d, p, nd;

  h = des->h;
  nd = dv->nd;
  wd = des->wd;
  d = lfd->d; p = des->p; X = d_x(des);
  ulx = des->res;
  m = des->n;
  for (i=0; i<d; i++) hs[i] = h*lfd->sca[i];
  if (deg>0)
  { lxd = &lx[m];
    setzero(lxd,m*d);
    if (deg>1)
    { lxdd = &lxd[d*m];
      setzero(lxdd,m*d*d);
  } }

  if (nd>0) fitfun(lfd,sp,des->xev,des->xev,des->f1,dv); /* c(0) */
    else unitvec(des->f1,0,p);
  jacob_solve(&des->xtwx,des->f1);   /* c(0) (X^TWX)^{-1} */
  for (i=0; i<m; i++)
  { ii = des->ind[i];
    lx[i] = innerprod(des->f1,&X[ii*p],p); /* c(0)(XTWX)^{-1}X^T */
    if (deg>0)
    { wd[i] = Wd(dist(des,ii)/h,ker(sp));
      for (j=0; j<d; j++)
      { dfx[j] = datum(lfd,j,ii)-des->xev[j];
        lxd[j*m+i] = lx[i]*wght(des,ii)*weightd(dfx[j],lfd->sca[j],
          d,ker(sp),kt(sp),h,lfd->sty[j],dist(des,ii));
             /* c(0) (XTWX)^{-1}XTW' */
      }
      if (deg>1)
      { wdd = Wdd(dist(des,ii)/h,ker(sp));
        for (j=0; j<d; j++)
          for (k=0; k<d; k++)
          { w = (dist(des,ii)==0) ? 0 : h/dist(des,ii);
            w = wdd * (des->xev[k]-datum(lfd,k,ii)) * (des->xev[j]-datum(lfd,j,ii))
                  * w*w / (hs[k]*hs[k]*hs[j]*hs[j]);
            if (j==k) w += wd[i]/(hs[j]*hs[j]);
            lxdd[(j*d+k)*m+i] = lx[i]*w;
              /* c(0)(XTWX)^{-1}XTW'' */
          }
      }
    }
    lx[i] *= wght(des,ii);
  }

  dv->nd = nd+1;
  if (deg==2)
  { for (i=0; i<d; i++)
    { dv->deriv[nd] = i;
      fitfun(lfd,sp,des->xev,des->xev,des->f1,dv);
      for (k=0; k<m; k++)
      { ii = des->ind[i];
        stdlinks(link,lfd,sp,ii,fitv(des,ii),robscale);
        for (j=0; j<p; j++)
          des->f1[j] -= link[ZDDLL]*lxd[i*m+k]*X[ii*p+j];
        /* c'(x)-c(x)(XTWX)^{-1}XTW'X */
      }
      jacob_solve(&des->xtwx,des->f1); /* (...)(XTWX)^{-1} */
      for (j=0; j<m; j++)
      { ii = des->ind[j];
        ulx[j] = innerprod(des->f1,&X[ii*p],p); /* (...)XT */
      }
      for (j=0; j<d; j++)
        for (k=0; k<m; k++)
        { ii = des->ind[k];
          dfx[j] = datum(lfd,j,ii)-des->xev[j];
          wdw = wght(des,ii)*weightd(dfx[j],lfd->sca[j],d,ker(sp),
            kt(sp),h,lfd->sty[j],dist(des,ii));
          lxdd[(i*d+j)*m+k] += ulx[k]*wdw;
          lxdd[(j*d+i)*m+k] += ulx[k]*wdw;
        } /* + 2(c'-c(XTWX)^{-1}XTW'X)(XTWX)^{-1}XTW' */
    }
    for (j=0; j<d*d; j++) nnresproj(lfd,sp,des,&lxdd[j*m],m,p);
        /* * (I-X(XTWX)^{-1} XTW */
  }
  if (deg>0)
  { for (j=0; j<d; j++) nnresproj(lfd,sp,des,&lxd[j*m],m,p);
      /* c(0)(XTWX)^{-1}XTW'(I-X(XTWX)^{-1}XTW) */
    for (i=0; i<d; i++)
    { dv->deriv[nd]=i;
      fitfun(lfd,sp,des->xev,des->xev,des->f1,dv);
      jacob_solve(&des->xtwx,des->f1);
      for (k=0; k<m; k++)
      { ii = des->ind[k];
        for (l=0; l<p; l++)
          lxd[i*m+k] += des->f1[l]*X[ii*p+l]*wght(des,ii);
      } /* add c'(0)(XTWX)^{-1}XTW */
    }
  }

  dv->nd = nd+2;
  if (deg==2)
  { for (i=0; i<d; i++)
    { dv->deriv[nd]=i;
      for (j=0; j<d; j++)
      { dv->deriv[nd+1]=j;
        fitfun(lfd,sp,des->xev,des->xev,des->f1,dv);
        jacob_solve(&des->xtwx,des->f1);
        for (k=0; k<m; k++)
        { ii = des->ind[k];
          for (l=0; l<p; l++)
            lxdd[(i*d+j)*m+k] += des->f1[l]*X[ii*p+l]*wght(des,ii);
        } /* + c''(x)(XTWX)^{-1}XTW */
      }
    }
  }
  dv->nd = nd;

  k = 1+d*(deg>0)+d*d*(deg==2);

  if (exp) wdexpand(lx,lfd->n,des->ind,m);
 
  if (ty==1) return(m);
  for (i=0; i<m; i++)
  { ii = des->ind[i];
    stdlinks(link,lfd,sp,ii,fitv(des,ii),robscale);
    link[ZDDLL] = sqrt(fabs(link[ZDDLL]));
    for (j=0; j<k; j++) lx[j*m+i] *= link[ZDDLL];
  }
  return(m);
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
/*
 *  String  matching functions. For a given argument string, find
 *  the best match from an array of possibilities. Is there a library
 *  function somewhere to do something like this?
 *
 *  return values of -1 indicate failure/unknown string.
 */

#include "locf.h"

int ct_match(z1, z2)
char *z1, *z2;
{ int ct = 0;
  while (z1[ct]==z2[ct])
  { if (z1[ct]=='\0') return(ct+1);
    ct++;
  }
  return(ct);
}

int pmatch(z, strings, vals, n, def)
char *z, **strings;
int *vals, n, def;
{ int i, ct, best, best_ct;
  best = -1;
  best_ct = 0;

  for (i=0; i<n; i++)
  { ct = ct_match(z,strings[i]);
    if (ct==strlen(z)+1) return(vals[i]);
    if (ct>best_ct) { best = i; best_ct = ct; }
  }
  if (best==-1) return(def);
  return(vals[best]);
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
#include "locf.h"

int lf_maxit = 20;
int lf_debug = 0;
int lf_error = 0;

double s0, s1;
static lfdata *lf_lfd;
static design *lf_des;
static smpar   *lf_sp;
int lf_status;
int ident=0;
double lf_tol;
extern double robscale;

void lfdata_init(lfd)
lfdata *lfd;
{ int i;
  for (i=0; i<MXDIM; i++)
  { lfd->sty[i] = 0;
    lfd->sca[i] = 1.0;
    lfd->xl[i] = lfd->xl[i+MXDIM] = 0.0;
  }
  lfd->y = lfd->w = lfd->c = lfd->b = NULL;
  lfd->d = lfd->n = 0;
}

void smpar_init(sp,lfd)
smpar *sp;
lfdata *lfd;
{ nn(sp)  = 0.7;
  fixh(sp)= 0.0;
  pen(sp) = 0.0;
  acri(sp)= ANONE;
  deg(sp) = deg0(sp) = 2;
  ubas(sp) = 0;
  kt(sp) = KSPH;
  ker(sp) = WTCUB;
  fam(sp) = 64+TGAUS;
  link(sp)= LDEFAU;
  npar(sp) = calcp(sp,lfd->d);
}

void deriv_init(dv)
deriv *dv;
{ dv->nd = 0;
}

int des_reqd(n,p)
int n, p;
{
  return(n*(p+5)+2*p*p+4*p + jac_reqd(p));
}
int des_reqi(n,p)
int n, p;
{ return(n+p);
}
 
void des_init(des,n,p)
design *des;
int n, p;
{ double *z;
  int k;

  if (n<=0) WARN(("des_init: n <= 0"));
  if (p<=0) WARN(("des_init: p <= 0"));

  if (des->des_init_id != DES_INIT_ID)
  { des->lwk = des->lind = 0;
    des->des_init_id = DES_INIT_ID;
  }

  k = des_reqd(n,p);
  if (k>des->lwk)
  { des->wk = (double *)calloc(k,sizeof(double));
    if ( des->wk == NULL ) {
      printf("Problem allocating memory for des->wk\n");fflush(stdout);
    }
    des->lwk = k;
  }
  z = des->wk;

  des->X = z; z += n*p;
  des->w = z; z += n;
  des->res=z; z += n;
  des->di =z; z += n;
  des->th =z; z += n;
  des->wd =z; z += n;
  des->V  =z; z += p*p;
  des->P  =z; z += p*p;
  des->f1 =z; z += p;
  des->ss =z; z += p;
  des->oc =z; z += p;
  des->cf =z; z += p;
 
  z = jac_alloc(&des->xtwx,p,z);
 
  k = des_reqi(n,p);
  if (k>des->lind)
  {
    des->ind = (int *)calloc(k,sizeof(int));
    if ( des->ind == NULL ) {
      printf("Problem allocating memory for des->ind\n");fflush(stdout);
    }
    des->lind = k;
  }
  des->fix = &des->ind[n];
  for (k=0; k<p; k++) des->fix[k] = 0;

  des->n = n; des->p = p;
  des->smwt = n;
  des->xtwx.p = p;                                                              
}

void deschk(des,n,p)
design *des;
int n, p;
{ WARN(("deschk deprecated - use des_init()"));
  des_init(des,n,p);
}

int likereg(coef, lk0, f1, Z)
double *coef, *lk0, *f1, *Z;
{ int i, ii, j, p;
  double lk, ww, link[LLEN], *X;

  if (lf_debug>2) mut_printf("  likereg: %8.5f\n",coef[0]);
  lf_status = LF_OK;
  lk = 0.0; p = lf_des->p;
  setzero(Z,p*p);
  setzero(f1,p);
  for (i=0; i<lf_des->n; i++)
  {
    ii = lf_des->ind[i];
    X = d_xi(lf_des,ii);
    fitv(lf_des,ii) = base(lf_lfd,ii)+innerprod(coef,X,p);
    lf_status = stdlinks(link,lf_lfd,lf_sp,ii,fitv(lf_des,ii),robscale);
    if (lf_status == LF_BADP)
    { *lk0 = -1.0e300;
      return(NR_REDUCE);
    }
    if (lf_error) lf_status = LF_ERR;
    if (lf_status != LF_OK) return(NR_BREAK);

    ww = wght(lf_des,ii);
    lk += ww*link[ZLIK];
    for (j=0; j<p; j++)
      f1[j] += X[j]*ww*link[ZDLL];
    addouter(Z, X, X, p, ww*link[ZDDLL]);
  }
  for (i=0; i<p; i++) if (lf_des->fix[i])
  { for (j=0; j<p; j++) Z[i*p+j] = Z[j*p+i] = 0.0;
    Z[i*p+i] = 1.0;
    f1[i] = 0.0;
  }

  if (lf_debug>4) prresp(coef,Z,p);
  if (lf_debug>3) mut_printf("  likelihood: %8.5f\n",lk);
  *lk0 = lf_des->llk = lk;

  lf_status = fami(lf_sp)->pcheck(lf_sp,lf_des,lf_lfd);
  switch(lf_status)
  { case LF_DONE: return(NR_BREAK);
    case LF_OOB:  return(NR_REDUCE);
    case LF_PF:   return(NR_REDUCE);
    case LF_NSLN: return(NR_BREAK);
  }

  return(NR_OK);
}

int reginit(lfd,des,sp)
lfdata *lfd;
design *des;
smpar *sp;
{ int i, ii;
  double sb, link[LLEN];
  s0 = s1 = sb = 0;
  for (i=0; i<des->n; i++)
  { ii = des->ind[i];
    links(base(lfd,ii),resp(lfd,ii),fami(sp),LINIT,link,cens(lfd,ii),prwt(lfd,ii),1.0);
    s1 += wght(des,ii)*link[ZDLL];
    s0 += wght(des,ii)*prwt(lfd,ii);
    sb += wght(des,ii)*prwt(lfd,ii)*base(lfd,ii);
  }
  if (s0==0) return(LF_NOPT); /* no observations with W>0 */
  setzero(des->cf,des->p);
  lf_tol = 1.0e-6*s0;
  switch(link(sp))
  { case LIDENT:
      des->cf[0] = (s1-sb)/s0;
      return(LF_OK);
    case LLOG:
      if (s1<=0.0)
      { des->cf[0] = -1000;
        return(LF_INFA);
      }
      des->cf[0] = log(s1/s0) - sb/s0;
      return(LF_OK);
    case LLOGIT:
      if (s1<=0.0)
      { des->cf[0] = -1000;
        return(LF_INFA);
      }
      if (s1>=s0)
      { des->cf[0] = 1000;
        return(LF_INFA);
      }
      des->cf[0] = logit(s1/s0)-sb/s0;
      return(LF_OK);
    case LINVER:
      if (s1<=0.0)
      { des->cf[0] = 1e100;
        return(LF_INFA);
      }
      des->cf[0] = s0/s1-sb/s0;
      return(LF_OK);
    case LSQRT:
      des->cf[0] = sqrt(s1/s0)-sb/s0;
      return(LF_OK);
    case LASIN:
      des->cf[0] = asin(sqrt(s1/s0))-sb/s0;
      return(LF_OK);
    default:
      LERR(("reginit: invalid link %d",link(sp)));
      return(LF_ERR);
  }
}

int lfinit(lfd,sp,des)
lfdata *lfd;
smpar *sp;
design *des;
{ int initstat;
  des->xtwx.sm = (deg0(sp)<deg(sp)) ? JAC_CHOL : JAC_EIGD;

  designmatrix(lfd,sp,des);
  setfamily(sp);
  initstat = fami(sp)->initial(lfd,des,sp);

  return(initstat);
}

void lfiter(lfd,sp,des,maxit)
lfdata *lfd;
smpar *sp;
design *des;
int maxit;
{ int err;
  if (lf_debug>1) mut_printf(" lfiter: %8.5f\n",des->cf[0]);

  lf_des = des;
  lf_lfd = lfd;
  lf_sp  = sp;

  max_nr(fami(sp)->like, des->cf, des->oc, des->res, des->f1,
    &des->xtwx, des->p, maxit, lf_tol, &err);
  switch(err)
  { case NR_OK: return;
    case NR_NCON:
      WARN(("max_nr not converged"));
      return;
    case NR_NDIV:
      WARN(("max_nr reduction problem"));
      return;
  }
  WARN(("max_nr return status %d",err));
}

int use_robust_scale(int tg)
{ if ((tg&64)==0) return(0); /* not quasi - no scale */
  if (((tg&128)==0) & (((tg&63)!=TROBT) & ((tg&63)!=TCAUC))) return(0);
  return(1);
}

/*
 * noit not really needed any more, since
 * gauss->pcheck returns LF_DONE, and likereg NR_BREAK
 * in gaussian case.
 * nb: 0/1: does local neighborhood and weights need computing?
 * cv: 0/1: is variance/covariance matrix needed?
 */
int locfit(lfd,des,sp,noit,nb,cv)
lfdata *lfd;
design *des;
smpar *sp;
int noit, nb, cv;
{ int i;

  if (des->xev==NULL)
  { LERR(("locfit: NULL evaluation point?"));
    return(246);
  }

  if (lf_debug>0)
  { mut_printf("locfit: ");
    for (i=0; i<lfd->d; i++) mut_printf(" %10.6f",des->xev[i]);
    mut_printf("\n");
  }

/* the 1e-12 avoids problems that can occur with roundoff */
  if (nb) nbhd(lfd,des,(int)(lfd->n*nn(sp)+1e-12),0,sp);

  lf_status = lfinit(lfd,sp,des);

  if (lf_status == LF_OK)
  { if (use_robust_scale(fam(sp)))
      lf_robust(lfd,sp,des,lf_maxit);
    else
    { if ((fam(sp)&63)==TQUANT)
        lfquantile(lfd,sp,des,lf_maxit);
      else
      { robscale = 1.0;
        lfiter(lfd,sp,des,lf_maxit);
      }
    }
  }

  if (lf_status == LF_DONE) lf_status = LF_OK;
  if (lf_status == LF_OOB) lf_status = LF_OK;

  if ((fam(sp)&63)==TDEN) /* convert from rate to density */
  { switch(link(sp))
    { case LLOG:
        des->cf[0] -= log(des->smwt);
        break;
      case LIDENT:
        multmatscal(des->cf,1.0/des->smwt,des->p);
        break;
      default: LERR(("Density adjustment; invalid link"));
    }
  }

  /* variance calculations, if requested */
  if (cv)
  { switch(lf_status)
    { case LF_PF:  /* for these cases, variance calc. would likely fail. */
      case LF_NOPT:
      case LF_NSLN:
      case LF_INFA:
      case LF_DEMP:
      case LF_XOOR:
      case LF_DNOP:
      case LF_BADP:
        des->llk = des->tr0 = des->tr1 = des->tr2 = 0.0;
        setzero(des->V,des->p*des->p);
        setzero(des->f1,des->p);
        break;
      default: lf_vcov(lfd,sp,des);
    }
  }

  return(lf_status);
}

void lf_status_msg(status)
int status;
{ switch(status)
{ case LF_OK: return;
  case LF_NCON: WARN(("locfit did not converge")); return;
  case LF_OOB: WARN(("parameters out of bounds")); return;
  case LF_PF: WARN(("perfect fit")); return;
  case LF_NOPT: WARN(("no points with non-zero weight")); return;
  case LF_NSLN: WARN(("no solution")); return;
  case LF_INFA: WARN(("initial value problem")); return;
  case LF_DEMP: WARN(("density estimate, empty integration region")); return;
  case LF_XOOR: WARN(("procv: fit point outside xlim region")); return;
  case LF_DNOP: WARN(("density estimation -- insufficient points in smoothing window")); return;
  case LF_BADP: WARN(("bad parameters")); return;
  default: WARN(("procv: unknown return code %d",status)); return;
} }
/*
 * Copyright 1996-2006 Catherine Loader.
 */
/*
 *   Compute minimax weights for local regression.
 */

#include "locf.h"
#define NR_EMPTY 834

int mmsm_ct;

static int debug=0;
#define CONVTOL 1.0e-8
#define SINGTOL 1.0e-10
#define NR_SINGULAR 100

static lfdata *mm_lfd;
static design *mm_des;
static double mm_gam, mmf, lb;
static int st;

double ipower(x,n) /* use for n not too large!! */
double x;
int n;
{ if (n==0) return(1.0);
  if (n<0) return(1/ipower(x,-n));
  return(x*ipower(x,n-1));
}

double setmmwt(des,a,gam)
design *des;
double *a, gam;
{ double ip, w0, w1, sw, wt;
  int i;
  sw = 0.0;
  for (i=0; i<mm_lfd->n; i++)
  { ip = innerprod(a,d_xi(des,i),des->p);
    wt = prwt(mm_lfd,i);
    w0 = ip - gam*des->wd[i];
    w1 = ip + gam*des->wd[i];
    wght(des,i) = 0.0;
    if (w0>0) { wght(des,i) = w0; sw += wt*w0*w0; }
    if (w1<0) { wght(des,i) = w1; sw += wt*w1*w1; }
  }
  return(sw/2-a[0]);
}

/* compute sum_{w!=0} AA^T; e1-sum wA  */
int mmsums(des,coef,f,z,J)
design *des;
double *coef, *f, *z;
jacobian *J;
{ int ct, i, j, p, sing;
  double *A;

mmsm_ct++;
  A = J->Z;
  *f = setmmwt(des,coef,mm_gam);

  p = des->p;
  setzero(A,p*p);
  setzero(z,p);
  z[0] = 1.0;
  ct = 0;

  for (i=0; i<mm_lfd->n; i++)
    if (wght(des,i)!=0.0)
    { addouter(A,d_xi(des,i),d_xi(des,i),p,prwt(mm_lfd,i));
      for (j=0; j<p; j++) z[j] -= prwt(mm_lfd,i)*wght(des,i)*d_xij(des,i,j);
      ct++;
    }
  if (ct==0) return(NR_EMPTY);

  J->st = JAC_RAW;
  J->p = p;
  jacob_dec(J,JAC_EIGD);

  sing = 0;
  for (i=0; i<p; i++) sing |= (J->Z[i*p+i]<SINGTOL);
  if ((debug) & (sing)) mut_printf("SINGULAR!!!!\n");

  return((sing) ? NR_SINGULAR : NR_OK);
}

int descenddir(des,coef,dlt,f,af)
design *des;
double *coef, *dlt, *f;
int af;
{ int i, p;
  double f0, *oc;

  if (debug) mut_printf("descenddir: %8.5f %8.5f\n",dlt[0],dlt[1]);

  f0 = *f;
  oc = des->oc;
  p = des->p;
  memcpy(oc,coef,p*sizeof(double));

  for (i=0; i<p; i++) coef[i] = oc[i]+lb*dlt[i];
  st = mmsums(des,coef,f,des->f1,&des->xtwx);

  if (*f>f0) /* halve till we drop */
  { while (*f>f0)
    { lb = lb/2.0;
      for (i=0; i<p; i++) coef[i] = oc[i]+lb*dlt[i];
      st = mmsums(des,coef,f,des->f1,&des->xtwx);
    }
    return(st);
  }

  if (!af) return(st);

  /* double */
  while (*f<f0)
  { f0 = *f;
    lb *= 2.0;
    for (i=0; i<p; i++) coef[i] = oc[i]+lb*dlt[i];
    st = mmsums(des,coef,f,des->f1,&des->xtwx);
  }

  lb /= 2.0;
  for (i=0; i<p; i++) coef[i] = oc[i]+lb*dlt[i];
  st = mmsums(des,coef,f,des->f1,&des->xtwx);

  return(st);
}

int mm_initial(des)
design *des;
{ double *dlt;

  dlt = des->ss;

  setzero(des->cf,des->p);
  st = mmsums(des,des->cf,&mmf,des->f1,&des->xtwx);

  setzero(dlt,des->p);
  dlt[0] = 1;
  lb = 1.0;
  st = descenddir(des,des->cf,dlt,&mmf,1);
  return(st);
}

void getsingdir(des,dlt)
design *des;
double *dlt;
{ double f, sw, c0;
  int i, j, p, sd;

  sd = -1; p = des->p;
  setzero(dlt,p);
  for (i=0; i<p; i++) if (des->xtwx.Z[i*p+i]<SINGTOL) sd = i;
  if (sd==-1)
  { mut_printf("getsingdir: nonsing?\n");
    return;
  }
  if (des->xtwx.dg[sd]>0)
    for (i=0; i<p; i++) dlt[i] = des->xtwx.Q[p*i+sd]*des->xtwx.dg[i];
  else
  { dlt[sd] = 1.0;
  }

  c0 = innerprod(dlt,des->f1,p);
  if (c0<0) for (i=0; i<p; i++) dlt[i] = -dlt[i];
}

void mmax(coef, old_coef, delta, J, p, maxit, tol, err)
double *coef, *old_coef, *delta, tol;
int p, maxit, *err;
jacobian *J;
{ double old_f, lambda;
  int i, j;

  *err = NR_OK;
 
  for (j=0; j<maxit; j++)
  { memcpy(old_coef,coef,p*sizeof(double));
    old_f = mmf;

    if (st == NR_SINGULAR)
    {
      getsingdir(mm_des,delta);
      st = descenddir(mm_des,coef,delta,&mmf,1);
    }
    if (st == NR_EMPTY)
    { 
      setzero(delta,p);
      delta[0] = 1.0;
      st = descenddir(mm_des,coef,delta,&mmf,1);
    }
    if (st == NR_OK)
    { 
      lb = 1.0;
      jacob_solve(J,mm_des->f1);
      memcpy(delta,mm_des->f1,p*sizeof(double));
      st = descenddir(mm_des,coef,delta,&mmf,0);
    }

    if ((j>0) & (fabs(mmf-old_f)<tol)) return;
  }
  WARN(("findab not converged"));
  *err = NR_NCON;
  return;
}

double findab(gam)
double gam;
{ double sl;
  int i, p, nr_stat;

  if (debug) mut_printf("  findab: gam %8.5f\n",gam);
  mm_gam = gam;
  p = mm_des->p;
  lb = 1.0;
  st = mm_initial(mm_des);

    mmax(mm_des->cf, mm_des->oc, mm_des->ss,
       &mm_des->xtwx, p, lf_maxit, CONVTOL, &nr_stat);

  sl = 0.0;
  for (i=0; i<mm_lfd->n; i++) sl += fabs(wght(mm_des,i))*mm_des->wd[i];

  if (debug) mut_printf("  sl %8.5f  gam %8.5f    %8.5f %d\n", sl,gam,sl-gam,nr_stat);
  return(sl-gam);
}

double weightmm(coef,di,ff,gam)
double *coef, di, *ff, gam;
{ double y1, y2, ip;
  ip = innerprod(ff,coef,mm_des->p);
  y1 = ip-gam*di; if (y1>0) return(y1/ip);
  y2 = ip+gam*di; if (y2<0) return(y2/ip);
  return(0.0);
}

double minmax(lfd,des,sp)
lfdata *lfd;
design *des;
smpar *sp;
{ double h, u[MXDIM], gam;
  int i, j, m, d1, p1, err_flag;

  if (debug) mut_printf("minmax: x %8.5f\n",des->xev[0]);
  mm_lfd = lfd;
  mm_des = des;

mmsm_ct = 0;
  d1 = deg(sp)+1;
  p1 = factorial(d1);
  for (i=0; i<lfd->n; i++)
  { for (j=0; j<lfd->d; j++) u[j] = datum(lfd,j,i);
    des->wd[i] = sp->nn/p1*ipower(dist(des,i),d1);
    des->ind[i] = i;
    fitfun(lfd, sp, u,des->xev,d_xi(des,i),NULL);
  }

/* find gamma (i.e. solve eqn 13.17 from book), using the secant method.
 * As a side effect, this finds the other minimax coefficients.
 * Note that 13.17 is rewritten as
 *   g2 = sum |l_i(x)| (||xi-x||^(p+1) M/(s*(p+1)!))
 * where g2 = gamma * s * (p+1)! / M. The gam variable below is g2.
 * The smoothing parameter is sp->nn == M/s.
 */
  gam = solve_secant(findab, 0.0, 0.0,1.0, 0.0000001, BDF_EXPRIGHT, &err_flag);

/*
 * Set the smoothing weights, in preparation for the actual fit.
 */
  h = 0.0; m = 0;
  for (i=0; i<lfd->n; i++)
  { wght(des,i) = weightmm(des->cf, des->wd[i],d_xi(des,i),gam);
    if (wght(des,i)>0)
    { if (dist(des,i)>h) h = dist(des,i);
      des->ind[m] = i;
      m++;
    }
  }
  des->n = m;
  return(h);
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
/*
 *
 *  Defines the weight functions and related quantities used
 *  in LOCFIT.
 */

#include "locf.h"

/*
 * convert kernel and kernel type strings to numeric codes.
 */
#define NWFUNS 13
static char *wfuns[NWFUNS] = {
  "rectangular", "epanechnikov", "bisquare",    "tricube",
  "triweight",   "gaussian",     "triangular",  "ququ",
  "6cub",        "minimax",      "exponential", "maclean", "parametric" };
static int wvals[NWFUNS] = { WRECT, WEPAN, WBISQ, WTCUB,
  WTRWT, WGAUS, WTRIA, WQUQU, W6CUB, WMINM, WEXPL, WMACL, WPARM };
int lfkernel(char *z)
{ return(pmatch(z, wfuns, wvals, NWFUNS, WTCUB));
}

#define NKTYPE 5
static char *ktype[NKTYPE] = { "spherical", "product", "center", "lm", "zeon" };
static int   kvals[NKTYPE] = { KSPH, KPROD, KCE, KLM, KZEON };
int lfketype(char *z)
{ return(pmatch(z, ktype, kvals, NKTYPE, KSPH));
}

/* The weight functions themselves.  Used everywhere. */
double W(u,ker)
double u;
int ker;
{ u = fabs(u);
  switch(ker)
  { case WRECT: return((u>1) ? 0.0 : 1.0);
    case WEPAN: return((u>1) ? 0.0 : 1-u*u);
    case WBISQ: if (u>1) return(0.0);
                u = 1-u*u; return(u*u);
    case WTCUB: if (u>1) return(0.0);
                u = 1-u*u*u; return(u*u*u);
    case WTRWT: if (u>1) return(0.0);
                u = 1-u*u; return(u*u*u);
    case WQUQU: if (u>1) return(0.0);
                u = 1-u*u; return(u*u*u*u);
    case WTRIA: if (u>1) return(0.0);
                return(1-u);
    case W6CUB: if (u>1) return(0.0);
                u = 1-u*u*u; u = u*u*u; return(u*u);
    case WGAUS: return(exp(-SQR(GFACT*u)/2.0));
    case WEXPL: return(exp(-EFACT*u));
    case WMACL: return(1/((u+1.0e-100)*(u+1.0e-100)));
    case WMINM: LERR(("WMINM in W"));
                return(0.0);
    case WPARM: return(1.0);
  }
  LERR(("W(): Unknown kernel %d\n",ker));
  return(1.0);
}

int iscompact(ker)
int ker;
{ if ((ker==WEXPL) | (ker==WGAUS) | (ker==WMACL) | (ker==WPARM)) return(0);
  return(1);
}

double weightprod(lfd,u,h,ker)
lfdata *lfd;
double *u, h;
int ker;
{ int i;
  double sc, w;
  w = 1.0;
  for (i=0; i<lfd->d; i++)
  { sc = lfd->sca[i];
    switch(lfd->sty[i])
    { case STLEFT:
        if (u[i]>0) return(0.0);
        w *= W(-u[i]/(h*sc),ker);
        break;
      case STRIGH:
        if (u[i]<0) return(0.0);
        w *= W(u[i]/(h*sc),ker);
        break;
      case STANGL:
        w *= W(2*fabs(sin(u[i]/(2*sc)))/h,ker);
        break;
      case STCPAR:
        break;
      default:
        w *= W(fabs(u[i])/(h*sc),ker);
    }
    if (w==0.0) return(w);
  }
  return(w);
}

double weightsph(lfd,u,h,ker, hasdi,di)
lfdata *lfd;
double *u, h, di;
int ker, hasdi;
{ int i;

  if (!hasdi) di = rho(u,lfd->sca,lfd->d,KSPH,lfd->sty);

  for (i=0; i<lfd->d; i++)
  { if ((lfd->sty[i]==STLEFT) && (u[i]>0.0)) return(0.0);
    if ((lfd->sty[i]==STRIGH) && (u[i]<0.0)) return(0.0);
  }
  if (h==0) return((di==0.0) ? 1.0 : 0.0);

  return(W(di/h,ker));
}

double weight(lfd,sp,x,t,h, hasdi,di)
lfdata *lfd;
smpar *sp;
double *x, *t, h, di;
int hasdi;
{ double u[MXDIM];
  int i;
  for (i=0; i<lfd->d; i++) u[i] = (t==NULL) ? x[i] : x[i]-t[i];
  switch(kt(sp))
  { case KPROD: return(weightprod(lfd,u,h,ker(sp)));
    case KSPH:  return(weightsph(lfd,u,h,ker(sp), hasdi,di));
  }
  LERR(("weight: unknown kernel type %d",kt(sp)));
  return(1.0);
}

double sgn(x)
double x;
{ if (x>0) return(1.0);
  if (x<0) return(-1.0);
  return(0.0);
}

double WdW(u,ker) /* W'(u)/W(u) */
double u;
int ker;
{ double eps=1.0e-10;
  if (ker==WGAUS) return(-GFACT*GFACT*u);
  if (ker==WPARM) return(0.0);
  if (fabs(u)>=1) return(0.0);
  switch(ker)
  { case WRECT: return(0.0);
    case WTRIA: return(-sgn(u)/(1-fabs(u)+eps));
    case WEPAN: return(-2*u/(1-u*u+eps));
    case WBISQ: return(-4*u/(1-u*u+eps));
    case WTRWT: return(-6*u/(1-u*u+eps));
    case WTCUB: return(-9*sgn(u)*u*u/(1-u*u*fabs(u)+eps));
    case WEXPL: return((u>0) ? -EFACT : EFACT);
  }
  LERR(("WdW: invalid kernel"));
  return(0.0);
}

/* deriv. weights .. spherical, product etc
   u, sc, sty needed only in relevant direction
   Acutally, returns (d/dx W(||x||/h) ) / W(.)
*/
double weightd(u,sc,d,ker,kt,h,sty,di)
double u, sc, h, di;
int d, ker, kt, sty;
{ if (sty==STANGL)
  { if (kt==KPROD)
      return(-WdW(2*sin(u/(2*sc)),ker)*cos(u/(2*sc))/(h*sc));
    if (di==0.0) return(0.0);
    return(-WdW(di/h,ker)*sin(u/sc)/(h*sc*di));
  }
  if (sty==STCPAR) return(0.0);
  if (kt==KPROD)
    return(-WdW(u/(h*sc),ker)/(h*sc));
  if (di==0.0) return(0.0);
  return(-WdW(di/h,ker)*u/(h*di*sc*sc));
}

double weightdd(u,sc,d,ker,kt,h,sty,di,i0,i1)
double *u, *sc, h, di;
int d, ker, kt, i0, i1, *sty;
{ double w;
  w = 1;
  if (kt==KPROD)
  {
    w = WdW(u[i0]/(h*sc[i0]),ker)*WdW(u[i1]/(h*sc[i1]),ker)/(h*h*sc[i0]*sc[i1]);
  }
  return(0.0);
}

/* Derivatives W'(u)/u.
   Used in simult. conf. band computations,
   and kernel density bandwidth selectors. */
double Wd(u,ker)
double u;
int ker;
{ double v;
  if (ker==WGAUS) return(-SQR(GFACT)*exp(-SQR(GFACT*u)/2));
  if (ker==WPARM) return(0.0);
  if (fabs(u)>1) return(0.0);
  switch(ker)
  { case WEPAN: return(-2.0);
    case WBISQ: return(-4*(1-u*u));
    case WTCUB: v = 1-u*u*u;
                return(-9*v*v*u);
    case WTRWT: v = 1-u*u;
                return(-6*v*v);
    default: LERR(("Invalid kernel %d in Wd",ker));
  }
  return(0.0);
}

/* Second derivatives W''(u)-W'(u)/u.
   used in simult. conf. band computations in >1 dimension. */
double Wdd(u,ker)
double u;
int ker;
{ double v;
  if (ker==WGAUS) return(SQR(u*GFACT*GFACT)*exp(-SQR(u*GFACT)/2));
  if (ker==WPARM) return(0.0);
  if (u>1) return(0.0);
  switch(ker)
  { case WBISQ: return(12*u*u);
    case WTCUB: v = 1-u*u*u;
                return(-9*u*v*v+54*u*u*u*u*v);
    case WTRWT: return(24*u*u*(1-u*u));
    default: LERR(("Invalid kernel %d in Wdd",ker));
  }
  return(0.0);
}

/* int u1^j1..ud^jd W(u) du.
   Used for local log-linear density estimation.
   Assume all j_i are even.
   Also in some bandwidth selection.
*/
double wint(d,j,nj,ker)
int d, *j, nj, ker;
{ double I, z;
  int k, dj;
  dj = d;
  for (k=0; k<nj; k++) dj += j[k];
  switch(ker) /* int_0^1 u^(dj-1) W(u)du  */
  { case WRECT: I = 1.0/dj; break;
    case WEPAN: I = 2.0/(dj*(dj+2)); break;
    case WBISQ: I = 8.0/(dj*(dj+2)*(dj+4)); break;
    case WTCUB: I = 162.0/(dj*(dj+3)*(dj+6)*(dj+9)); break;
    case WTRWT: I = 48.0/(dj*(dj+2)*(dj+4)*(dj+6)); break;
    case WTRIA: I = 1.0/(dj*(dj+1)); break;
    case WQUQU: I = 384.0/(dj*(dj+2)*(dj+4)*(dj+6)*(dj+8)); break;
    case W6CUB: I = 524880.0/(dj*(dj+3)*(dj+6)*(dj+9)*(dj+12)*(dj+15)*(dj+18)); break;
    case WGAUS: switch(d)
                { case 1: I = S2PI/GFACT; break;
                  case 2: I = 2*PI/(GFACT*GFACT); break;
                  default: I = exp(d*log(S2PI/GFACT)); /* for nj=0 */
                }
                for (k=0; k<nj; k++) /* deliberate drop */
                  switch(j[k])
                  { case 4: I *= 3.0/(GFACT*GFACT);
                    case 2: I /= GFACT*GFACT;
                  }
                return(I);
    case WEXPL: I = factorial(dj-1)/ipower(EFACT,dj); break;
    default: LERR(("Unknown kernel %d in exacint",ker));
  }
  if ((d==1) && (nj==0)) return(2*I); /* common case quick */
  z = (d-nj)*LOGPI/2-mut_lgammai(dj);
  for (k=0; k<nj; k++) z += mut_lgammai(j[k]+1);
  return(2*I*exp(z));
}

/* taylor series expansion of weight function around x.
   0 and 1 are common arguments, so are worth programming
   as special cases.
   Used in density estimation.
*/
int wtaylor(f,x,ker)
double *f, x;
int ker;
{ double v;
  switch(ker)
  { case WRECT:
      f[0] = 1.0;
      return(1);
    case WEPAN:
      f[0] = 1-x*x; f[1] = -2*x; f[2] = -1;
      return(3);
    case WBISQ:
      v = 1-x*x;
      f[0] = v*v;   f[1] = -4*x*v; f[2] = 4-6*v;
      f[3] = 4*x;   f[4] = 1;
      return(5);
    case WTCUB:
      if (x==1.0)
      { f[0] = f[1] = f[2] = 0; f[3] = -27; f[4] = -81; f[5] = -108;
        f[6] = -81; f[7] = -36; f[8] = -9; f[9] = -1; return(10); }
      if (x==0.0)
      { f[1] = f[2] = f[4] = f[5] = f[7] = f[8] = 0;
        f[0] = 1; f[3] = -3; f[6] = 3; f[9] = -1; return(10); }
      v = 1-x*x*x;
      f[0] = v*v*v; f[1] = -9*v*v*x*x; f[2] = x*v*(27-36*v);
      f[3] = -27+v*(108-84*v);         f[4] = -3*x*x*(27-42*v);
      f[5] = x*(-108+126*v);           f[6] = -81+84*v;
      f[7] = -36*x*x; f[8] = -9*x;     f[9] = -1;
      return(10);
    case WTRWT:
      v = 1-x*x;
      f[0] = v*v*v; f[1] = -6*x*v*v; f[2] = v*(12-15*v);
      f[3] = x*(20*v-8); f[4] = 15*v-12; f[5] = -6; f[6] = -1;
      return(7);
    case WTRIA:
      f[0] = 1-x; f[1] = -1;
      return(2);
    case WQUQU:
      v = 1-x*x;
      f[0] = v*v*v*v; f[1] = -8*x*v*v*v; f[2] = v*v*(24-28*v);
      f[3] = v*x*(56*v-32); f[4] = (70*v-80)*v+16; f[5] = x*(32-56*v);
      f[6] = 24-28*v; f[7] = 8*x; f[8] = 1;
      return(9);
    case W6CUB:
      v = 1-x*x*x;
      f[0] = v*v*v*v*v*v;
      f[1] = -18*x*x*v*v*v*v*v;
      f[2] = x*v*v*v*v*(135-153*v);
      f[3] = v*v*v*(-540+v*(1350-816*v));
      f[4] = x*x*v*v*(1215-v*(4050-v*3060));
      f[5] = x*v*(-1458+v*(9234+v*(-16254+v*8568)));
      f[6] = 729-v*(10206-v*(35154-v*(44226-v*18564)));
      f[7] = x*x*(4374-v*(30132-v*(56862-v*31824)));
      f[8] = x*(12393-v*(61479-v*(92664-v*43758)));
      f[9] = 21870-v*(89100-v*(115830-v*48620));
      f[10]= x*x*(26730-v*(69498-v*43758));
      f[11]= x*(23814-v*(55458-v*31824));
      f[12]= 15849-v*(34398-v*18564);
      f[13]= x*x*(7938-8568*v);
      f[14]= x*(2970-3060*v);
      f[15]= 810-816*v;
      f[16]= 153*x*x;
      f[17]= 18*x;
      f[18]= 1;
      return(19);
  }
  LERR(("Invalid kernel %d in wtaylor",ker));
  return(0);
}

/* convolution int W(x)W(x+v)dx.
   used in kde bandwidth selection.
*/
double Wconv(v,ker)
double v;
int ker;
{ double v2;
  switch(ker)
  { case WGAUS: return(SQRPI/GFACT*exp(-SQR(GFACT*v)/4));
    case WRECT:
      v = fabs(v);
      if (v>2) return(0.0);
      return(2-v);
    case WEPAN:
      v = fabs(v);
      if (v>2) return(0.0);
      return((2-v)*(16+v*(8-v*(16-v*(2+v))))/30);
    case WBISQ:
      v = fabs(v);
      if (v>2) return(0.0);
      v2 = 2-v;
      return(v2*v2*v2*v2*v2*(16+v*(40+v*(36+v*(10+v))))/630);
  }
  LERR(("Wconv not implemented for kernel %d",ker));
  return(0.0);
}

/* derivative of Wconv.
   1/v d/dv int W(x)W(x+v)dx
   used in kde bandwidth selection.
*/
double Wconv1(v,ker)
double v;
int ker;
{ double v2;
  v = fabs(v);
  switch(ker)
  { case WGAUS: return(-0.5*SQRPI*GFACT*exp(-SQR(GFACT*v)/4));
    case WRECT:
      if (v>2) return(0.0);
      return(1.0);
    case WEPAN:
      if (v>2) return(0.0);
      return((-16+v*(12-v*v))/6);
    case WBISQ:
      if (v>2) return(0.0);
      v2 = 2-v;
      return(-v2*v2*v2*v2*(32+v*(64+v*(24+v*3)))/210);
  }
  LERR(("Wconv1 not implemented for kernel %d",ker));
  return(0.0);
}

/* 4th derivative of Wconv.
   used in kde bandwidth selection (BCV, SJPI, GKK)
*/
double Wconv4(v,ker)
double v;
int ker;
{ double gv;
  switch(ker)
  { case WGAUS:
      gv = GFACT*v;
      return(exp(-SQR(gv)/4)*GFACT*GFACT*GFACT*(12-gv*gv*(12-gv*gv))*SQRPI/16);
  }
  LERR(("Wconv4 not implemented for kernel %d",ker));
  return(0.0);
}

/* 5th derivative of Wconv.
   used in kde bandwidth selection (BCV method only)
*/
double Wconv5(v,ker) /* (d/dv)^5 int W(x)W(x+v)dx */
double v;
int ker;
{ double gv;
  switch(ker)
  { case WGAUS:
      gv = GFACT*v;
      return(-exp(-SQR(gv)/4)*GFACT*GFACT*GFACT*GFACT*gv*(60-gv*gv*(20-gv*gv))*SQRPI/32);
  }
  LERR(("Wconv5 not implemented for kernel %d",ker));
  return(0.0);
}

/* 6th derivative of Wconv.
   used in kde bandwidth selection (SJPI)
*/
double Wconv6(v,ker)
double v;
int ker;
{ double gv, z;
  switch(ker)
  { case WGAUS:
      gv = GFACT*v;
      gv = gv*gv;
      z = exp(-gv/4)*(-120+gv*(180-gv*(30-gv)))*0.02769459142;
      gv = GFACT*GFACT;
      return(z*gv*gv*GFACT);
  }
  LERR(("Wconv6 not implemented for kernel %d",ker));
  return(0.0);
}

/* int W(v)^2 dv / (int v^2 W(v) dv)^2
   used in some bandwidth selectors
*/
double Wikk(ker,deg)
int ker, deg;
{ switch(deg)
  { case 0:
    case 1: /* int W(v)^2 dv / (int v^2 W(v) dv)^2 */
      switch(ker)
      { case WRECT: return(4.5);
        case WEPAN: return(15.0);
        case WBISQ: return(35.0);
        case WGAUS: return(0.2820947918*GFACT*GFACT*GFACT*GFACT*GFACT);
        case WTCUB: return(34.152111046847892);   /* 59049 / 1729 */
        case WTRWT: return(66.083916083916080);   /* 9450/143 */
      }
    case 2:
    case 3: /* 4!^2/8*int(W1^2)/int(v^4W1)^2
               W1=W*(n4-v^2n2)/(n0n4-n2n2) */
      switch(ker)
      { case WRECT: return(11025.0);
        case WEPAN: return(39690.0);
        case WBISQ: return(110346.9231);
        case WGAUS: return(14527.43412);
        case WTCUB: return(126500.5904);
        case WTRWT: return(254371.7647);
      }
  }
  LERR(("Wikk not implemented for kernel %d, deg %d",ker,deg));
  return(0.0);
}
