/*
 * Copyright 1996-2006 Catherine Loader.
 */

#include "mex.h"
/*
 * Copyright 1996-2006 Catherine Loader.
 */
/*
 *   Copyright (c) 1998-2006 Catherine Loader.
 *   This file contains functions to compute the constants
 *   appearing in the tube formula.
 */

#include <stdio.h>
#include <math.h>
#include "tube.h"

static double *fd, *ft;
static int globm, (*wdf)(), use_covar, kap_terms;

int k0_reqd(d,n,uc)
int d, n, uc;
{ int m;
  m = d*(d+1)+1;
  if (uc) return(2*m*m);
     else return(2*n*m);
}

void assignk0(z,d,n) /* z should be n*(2*d*d+2*d+2); */
double *z;
int d, n;
{ ft = z; z += n*(d*(d+1)+1);
  fd = z; z += n*(d*(d+1)+1);
}

/* Residual projection of y to the columns of A,
 * (I - A(R^TR)^{-1}A^T)y
 * R should be from the QR-decomp. of A.
 */
void rproject(y,A,R,n,p)
double *y, *A, *R;
int n, p;
{ double v[1+TUBE_MXDIM];
  int i, j;

  for (i=0; i<p; i++) v[i] = innerprod(&A[i*n],y,n);
  qrsolv(R,v,n,p);
  for (i=0; i<n; i++)
    for (j=0; j<p; j++)
      y[i] -= A[j*n+i]*v[j];
}

double k2c(lij,A,m,dd,d)
double *lij, *A;
int m, d, dd;
{ int i, j, k, l;
  double sum, *bk, v[TUBE_MXDIM];

  for (i=0; i<dd*d; i++)
    chol_hsolve(fd,&lij[i*m],m,dd+1);
  for (i=0; i<dd*d; i++)
    for (j=0; j<dd*d; j++)
      lij[i*m+j+d+1] -= innerprod(&lij[i*m],&lij[j*m],dd+1);

  sum = 0;
  for (i=0; i<dd; i++)
    for (j=0; j<i; j++)
    { bk = &lij[i*d*m + j*d + d+1];
      for (k=0; k<dd; k++)
      { v[0] = 0;
        for (l=0; l<dd; l++) v[l+1] = bk[k*m+l];
        chol_solve(fd,v,m,dd+1);
        for (l=0; l<dd; l++) bk[k*m+l] = v[l+1];
      }
      for (k=0; k<dd; k++)
      { v[0] = 0;
        for (l=0; l<dd; l++) v[l+1] = bk[l*m+k];
        chol_solve(fd,v,m,dd+1);
        for (l=0; l<dd; l++) bk[l*m+k] = v[l+1];
      }
      sum += bk[i*m+j] - bk[j*m+i];
    }
  return(sum*fd[0]*fd[0]);
}

double k2x(lij,A,m,d,dd)
double *lij, *A;
int m, d, dd;
{ int i, j, k;
  double s, v[1+TUBE_MXDIM], *ll;

/* residual projections of lij onto A = [l,l1,...,ld] */
  for (i=0; i<d; i++)
    for (j=i; j<d; j++)
    { ll = &lij[(i*dd+j)*m];
      rproject(ll,A,fd,m,d+1);
      if (i!=j) memcpy(&lij[(j*dd+i)*m],ll,m*sizeof(double));
    }

/* compute lij[j][i] = e_i^T (A^T A)^{-1} B_j^T */
  for (k=0; k<m; k++)
    for (j=0; j<d; j++)
    { v[0] = 0;
      for (i=0; i<d; i++) v[i+1] = lij[(j*dd+i)*m+k];
      qrsolv(fd,v,m,d+1);
      for (i=0; i<d; i++) lij[(j*dd+i)*m+k] = v[i+1];
    } 

/* finally, add up to get the kappa2 term */
  s = 0;
  for (j=0; j<d; j++)
    for (k=0; k<j; k++)
      s += innerprod(&lij[(j*dd+j)*m],&lij[(k*dd+k)*m],m)
         - innerprod(&lij[(j*dd+k)*m],&lij[(k*dd+j)*m],m);

  return(s*fd[0]*fd[0]);
}

void d2c(ll,nn,li,ni,lij,nij,M,m,dd,d)
double *ll, *nn, *li, *ni, *lij, *nij, *M;
int m, dd, d;
{ int i, j, k, l, t, u, v, w;
  double z;

  for (i=0; i<dd; i++)
    for (j=0; j<dd; j++)
    { for (k=0; k<d; k++)
      { for (l=0; l<d; l++)
        { z = M[i*d+k]*M[j*d+l];
          if (z != 0.0)
          { nij[(i*d+j)*m] += z*lij[(k*d+l)*m];
            for (t=0; t<d; t++) /* need d, not dd here */
              for (u=0; u<d; u++)
                nij[(i*d+j)*m+t+1] += z*M[t*d+u]*lij[(k*d+l)*m+u+1];
            for (t=0; t<dd; t++)
              for (u=0; u<dd; u++)
              { for (v=0; v<d; v++)
                  for (w=0; w<d; w++)
                    nij[(i*d+j)*m+(t*d+u)+d+1] +=
                      z*M[t*d+v]*M[u*d+w]*lij[(k*d+l)*m+(v*d+w)+d+1];
                for (v=0; v<d; v++)
                  nij[(i*d+j)*m+(t*d+u)+d+1] += z*M[(v+1)*d*d+t*d+u]*lij[(k*d+l)*m+v+1];
              }
          }
        }

        z = M[(k+1)*d*d+i*d+j];
        if (z!=0.0)
        { nij[(i*d+j)*m] += z*li[k*m];
          for (t=0; t<d; t++)
            for (u=0; u<d; u++)
              nij[(i*d+j)*m+t+1] += z*M[t*d+u]*li[k*m+u+1];
          for (t=0; t<dd; t++)
            for (u=0; u<dd; u++)
            { for (v=0; v<d; v++)
                for (w=0; w<d; w++)
                  nij[(i*d+j)*m+(t*d+u)+d+1] += z*M[t*d+v]*M[u*d+w]*lij[(v*d+w)*m+k+1];
              for (v=0; v<d; v++)
                nij[(i*d+j)*m+(t*d+u)+d+1] += z*M[(v+1)*d*d+t*d+u]*li[k*m+v+1];
            }
        }
      }
    }
}

void d2x(li,lij,nij,M,m,dd,d)
double *li, *lij, *nij, *M;
int m, dd, d;
{ int i, j, k, l, z;
  double t;
  for (i=0; i<dd; i++)
    for (j=0; j<dd; j++)
    { for (k=0; k<d; k++)
      { for (l=0; l<d; l++)
        { t = M[i*d+k] * M[j*d+l];
          if (t != 0.0)
          { for (z=0; z<m; z++)
              nij[(i*d+j)*m+z] += t*lij[(k*d+l)*m+z];
          }
        }
        t = M[(k+1)*d*d+i*d+j];
        if (t!=0.0)
          for (z=0; z<m; z++)
            nij[(i*d+j)*m+z] += t*li[k*m+z];
      }
    }
}

int k0x(x,d,kap,M)
double *x, *kap, *M;
int d;
{ double det, *lij, *nij, z;
  int j, m, r;

  r = 1 + ((d>=2) & (kap_terms >= 3));
  m = globm = wdf(x,ft,r);

  memcpy(fd,ft,m*(d+1)*sizeof(double));
  if (use_covar) chol_dec(fd,m,d+1);
            else qr(fd,m,d+1,NULL);

  det = 1;
  for (j=1; j<=d; j++)
    det *= fd[j*(m+1)]/fd[0];
  kap[0] = det;
  if (kap_terms == 1) return(1);
  kap[1] = 0.0;
  if ((kap_terms == 2) | (d<=1)) return(2);

  lij = &ft[(d+1)*m];
  nij = &fd[(d+1)*m];
  memcpy(nij,lij,m*d*d*sizeof(double));
  z = (use_covar) ?  k2c(nij,ft,m,d,d) : k2x(nij,ft,m,d,d);
  kap[2] = z*det;
  if ((kap_terms == 3) | (d==2)) return(3);

  kap[3] = 0;
  return(4);
}

void d1c(li,ni,m,d,M)
double *li, *ni, *M;
int m, d;
{ int i, j, k, l;
  double t;

  fd[0] = ft[0];
  for (i=0; i<d; i++)
  { t = 0;
    for (j=0; j<d; j++) t += M[i*d+j]*li[j*m];
    fd[i+1] = ni[i*m] = t;

    for (j=0; j<d; j++)
    { t = 0;
      for (k=0; k<d; k++)
        for (l=0; l<d; l++)
          t += li[k*m+l+1] * M[i*d+k] * M[j*d+l];
      ni[i*m+j+1] = t;
    }
  }
}

void d1x(li,ni,m,d,M)
double *li, *ni, *M;
int m, d;
{ int i, j, k;
  memcpy(fd,ft,m*sizeof(double));
  setzero(ni,m*d);
  for (j=0; j<d; j++)
    for (k=0; k<d; k++)
      if (M[j*d+k]!=0)
        for (i=0; i<m; i++) ni[j*m+i] += M[j*d+k]*li[k*m+i];
}

int l1x(x,d,lap,M)
double *x, *lap, *M;
int d;
{ double det, sumcj, *u, v[TUBE_MXDIM];
  double *ll, *li, *lij, *ni, *nij;
  int i, j, m;
  if (kap_terms<=1) return(0);
  m = globm;
  li  = &ft[m]; lij = &ft[(d+1)*m];
  ni  = &fd[m]; nij = &fd[(d+1)*m];
  setzero(ni,m*d);
  setzero(nij,m*d*d);

  if (use_covar) d1c(li,ni,m,d,M);
            else d1x(li,ni,m,d,M);

/* the last (d+1) columns of nij are free, use for an extra copy of ni */
  ll = &fd[d*d*m];
  u = &ll[d*m];
  if (use_covar)
    memcpy(u,&ni[(d-1)*m],d*sizeof(double)); /* cov(ld, (l,l1,...ld-1)) */
  else
    memcpy(ll,fd,(d+1)*m*sizeof(double));

  if (use_covar) chol_dec(fd,m,d+1);
            else qr(fd,m,d+1,NULL);
  det = 1;
  for (j=1; j<d; j++)
    det *= fd[(m+1)*j]/fd[0];
  lap[0] = det;
  if ((kap_terms==2) | (d<=1)) return(1);

  sumcj = 0.0;
  if (use_covar)
  { d2c(ft,fd,li,ni,lij,nij,M,m,d-1,d);
    chol_solve(fd,u,m,d);
    for (i=0; i<d-1; i++)
    { v[0] = 0;
      for (j=0; j<d-1; j++)
        v[j+1] = nij[(i*d+j)*m+d] - innerprod(u,&nij[(i*d+j)*m],d);
      chol_solve(fd,v,m,d);
      sumcj -= v[i+1];
    }
  }
  else
  { d2x(li,lij,nij,M,m,d-1,d);
    rproject(u,ll,fd,m,d);
    for (i=0; i<d-1; i++)
    { v[0] = 0;
      for (j=0; j<d-1; j++) v[j+1] = innerprod(&nij[(i*d+j)*m],u,m);
      qrsolv(fd,v,m,d);
      sumcj -= v[i+1];
    }
  }

  lap[1] = sumcj*det*fd[0]/fd[(m+1)*d];
  if ((kap_terms==3) | (d==2)) return(2);

  if (use_covar) lap[2] = k2c(nij,ll,m,d-1,d)*det;
            else lap[2] = k2x(nij,ll,m,d-1,d)*det;
  return(3);
}

int m0x(x,d,m0,M)
double *x, *m0, *M;
int d;
{ double det, *li, *ni, *lij, *nij, *ll, *u1, *u2;
  double om, so, co, sumcj, v[TUBE_MXDIM];
  int m, i, j;
  
  if ((kap_terms<=2) | (d<=1)) return(0);

  m = globm;
  li  = &ft[m]; lij = &ft[(d+1)*m];
  ni  = &fd[m]; nij = &fd[(d+1)*m];
  setzero(ni,m*d); 
  setzero(nij,m*d*d);

  if (use_covar) d1c(li,ni,m,d,M);
            else d1x(li,ni,m,d,M);

/* the last (d+1) columns of nij are free, use for an extra copy of ni */
  ll = &fd[d*d*m];
  u1 = &ll[d*m];
  u2 = &ll[(d-1)*m];
  if (use_covar)
  { memcpy(u1,&ni[(d-1)*m],d*sizeof(double));
    memcpy(u2,&ni[(d-2)*m],d*sizeof(double));
  }
  else
    memcpy(ll,fd,(d+1)*m*sizeof(double));

  if (use_covar) chol_dec(fd,m,d+1);
            else qr(fd,m,d+1,NULL);
  det = 1;
  for (j=1; j<d-1; j++)
    det *= fd[j*(m+1)]/fd[0];
  om = atan2(fd[d*(m+1)],-fd[d*(m+1)-1]);
  m0[0] = det*om;
  if ((kap_terms==3) | (d==2)) return(1);

  so = sin(om)/fd[d*(m+1)];
  co = (1-cos(om))/fd[(d-1)*(m+1)];
  sumcj = 0.0;
  if (use_covar)
  { d2c(ft,fd,li,ni,lij,nij,M,m,d-2,d);
    chol_solve(fd,u1,m,d);
    chol_solve(fd,u2,m,d-1);
    for (i=0; i<d-2; i++)
    { v[0] = 0;
      for (j=0; j<d-2; j++)
        v[j+1] =
          so*(nij[(i*d+j)*m+d]-innerprod(u1,&nij[(i*d+j)*m],d))
         +co*(nij[(i*d+j)*m+d-1]-innerprod(u2,&nij[(i*d+j)*m],d-1));
      qrsolv(fd,v,m,d-1);
      sumcj -= v[i+1];
    }
  }
  else
  { d2x(li,lij,nij,M,m,d-2,d);
    rproject(u1,ll,fd,m,d);
    rproject(u2,ll,fd,m,d-1); /* now, u1, u2 are unnormalized n1*, n2* */
    for (i=0; i<m; i++)
      u1[i] = so*u1[i] + co*u2[i];  /* for n1*, n2* */
    for (i=0; i<d-2; i++)
    { v[0] = 0;
      for (j=0; j<d-2; j++)
        v[j+1] = innerprod(&nij[(i*d+j)*m],u1,m);
      qrsolv(fd,v,m,d-1);
      sumcj -= v[i+1];
    }
  }

  m0[1] = sumcj*det*fd[0];
  return(2);
}

int n0x(x,d,n0,M)
double *x, *n0, *M;
int d;
{ double det, *li, *ni, *a0, *a1, *a2;
  int j, m;

  if ((kap_terms <= 3) | (d <= 2)) return(0);

  m = globm;
  li  = &ft[m];
  ni  = &fd[m];

  if (use_covar) d1c(li,ni,m,d,M);
            else d1x(li,ni,m,d,M);

  det = 1;
  if (use_covar) chol_dec(fd,m,d+1);
            else qr(fd,m,d+1,NULL);
  for (j=1; j<d-2; j++)
    det *= fd[j*(m+1)]/fd[0];

  a0 = &ni[(d-3)*m+d-2];
  a1 = &ni[(d-2)*m+d-2];
  a2 = &ni[(d-1)*m+d-2];

  a0[0] = a1[1]*a2[2];
  a0[1] =-a1[0]*a2[2];
  a0[2] = a1[0]*a2[1]-a1[1]*a2[0];
  a1[0] = 0;
  a1[1] = a2[2];
  a1[2] =-a2[1];
  a2[0] = a2[1] = 0.0; a2[2] = 1.0;
  rn3(a0); rn3(a1);
  n0[0] = det*sptarea(a0,a1,a2);
  return(1);
}

int kodf(ll,ur,mg,kap,lap)
double *ll, *ur, *kap, *lap;
int *mg;
{ double x[1], *l0, *l1, t, sum;
  int i, j, n;

  sum = 0.0;
  for (i=0; i<=mg[0]; i++)
  { if (i&1) { l1 = fd; l0 = ft; }
        else { l1 = ft; l0 = fd; }
    x[0] = ll[0] + (ur[0]-ll[0])*i/mg[0];
    n = wdf(x,l0,0);

    t = sqrt(innerprod(l0,l0,n));
    for (j=0; j<n; j++) l0[j] /= t;

    if (i>0)
    { t = 0.0;
      for (j=0; j<n; j++) t += (l1[j]-l0[j])*(l1[j]-l0[j]);
      sum += 2*asin(sqrt(t)/2);
    }
  }
  kap[0] = sum;
  if (kap_terms<=1) return(1);
  kap[1] = 0.0;
  lap[0] = 2.0;
  return(2);
}

int tube_constants(f,d,m,ev,mg,fl,kap,wk,terms,uc)
double *fl, *kap, *wk;
int d, m, ev, *mg, (*f)(), terms, uc;
{ int aw, deb=0;
  double k0[4], l0[3], m0[2], n0[1], z[TUBE_MXDIM];

  wdf = f;

  aw = (wk==NULL);
  if (aw) {
    wk = (double *)calloc(k0_reqd(d,m,uc),sizeof(double));
    if ( wk == NULL ) {
      printf("Problem allocating memory for wk\n");fflush(stdout);
    }
  }
  assignk0(wk,d,m);

  k0[0] = k0[1] = k0[2] = k0[3] = 0.0;
  l0[0] = l0[1] = l0[2] = 0.0;
  m0[0] = m0[1] = 0.0;
  n0[0] = 0.0;

  use_covar = uc;
  kap_terms = terms;
  if ((kap_terms <=0) | (kap_terms >= 5))
    mut_printf("Warning: terms = %2d\n",kap_terms);

  switch(ev)
  {
    case IMONTE:
      monte(k0x,fl,&fl[d],d,k0,mg[0]);
      break;
    case ISPHERIC:
      if (d==2) integ_disc(k0x,l1x,fl,k0,l0,mg);
      if (d==3) integ_sphere(k0x,l1x,fl,k0,l0,mg);
      break;
    case ISIMPSON:
      if (use_covar) simpson4(k0x,l1x,m0x,n0x,fl,&fl[d],d,k0,l0,m0,n0,mg,z);
                else simpson4(k0x,l1x,m0x,n0x,fl,&fl[d],d,k0,l0,m0,n0,mg,z);
      break;
    case IDERFREE:
      kodf(fl,&fl[d],mg,k0,l0);
      break;
    default:
      mut_printf("Unknown integration type in tube_constants().\n");
  }

  if (deb>0)
  { mut_printf("constants:\n");
    mut_printf("  k0: %8.5f %8.5f %8.5f %8.5f\n",k0[0],k0[1],k0[2],k0[3]);
    mut_printf("  l0: %8.5f %8.5f %8.5f\n",l0[0],l0[1],l0[2]);
    mut_printf("  m0: %8.5f %8.5f\n",m0[0],m0[1]);
    mut_printf("  n0: %8.5f\n",n0[0]);
    if (d==2) mut_printf("  check: %8.5f\n",(k0[0]+k0[2]+l0[1]+m0[0])/(2*PI));
    if (d==3) mut_printf("  check: %8.5f\n",(l0[0]+l0[2]+m0[1]+n0[0])/(4*PI));
  }

  if (aw) free(wk);

  kap[0] = k0[0];
  if (kap_terms==1) return(1);
  kap[1] = l0[0]/2;
  if ((kap_terms==2) | (d==1)) return(2);
  kap[2] = (k0[2]+l0[1]+m0[0])/(2*PI);
  if ((kap_terms==3) | (d==2)) return(3);
  kap[3] = (l0[2]+m0[1]+n0[0])/(4*PI);
  return(4);
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
/*
 *   Copyright (c) 1998-2006 Catherine Loader.
 *
 *   Computes the critical values from constants kappa0 etc
 *   and significance level.
 */

#include <math.h>
#include "tube.h"

#define LOGPI 1.144729885849400174143427

/* area(d) = 2 pi^(d/2) / Gamma(d/2)
 *         = surface area of unit sphere in R^d
 */
static double A[10] =
  { 1, /* d=0, whatever */
    2,
    6.2831853071795864770, /* 2*pi */
    12.566370614359172954, /* 4*pi */
    19.739208802178717238, /* 2*pi^2 */
    26.318945069571622985, /* 8/3*pi^2 */
    31.006276680299820177, /* pi^3 */
    33.073361792319808190, /* 16/15*pi^3 */
    32.469697011334145747, /* 1/3*pi^4 */
    29.686580124648361825  /* 32/105*pi^4 */
  };

double area(d)
int d;
{ if (d<10) return(A[d]);
  return(2*exp(d*LOGPI/2.0-mut_lgammai(d)));
}

double tailp_uniform(c,k0,m,d,s,n)
double c, *k0, n;
int m, d, s;
{ int i;
  double p;
  p = 0.0;
  for (i=0; i<m; i++) if (k0[i] != 0.0)
    p += k0[i] * ibeta(1-c*c,(n-d+i-1)/2.0,(d+1-i)/2.0) / area(d+1-i);
  return( (s==TWO_SIDED) ? 2*p : p );
}

double tailp_gaussian(c,k0,m,d,s,n)
double c, *k0, n;
int m, d, s;
{ int i;
  double p;
  p = 0.0;
  for (i=0; i<m; i++) if (k0[i] != 0.0)
    p += k0[i] * (1-pchisq(c*c,(double) d+1-i)) / area(d+1-i);
  return( (s==TWO_SIDED) ? 2*p : p );
}

double tailp_tprocess(c,k0,m,d,s,n)
double c, *k0, n;
int m, d, s;
{ int i;
  double p;
  p = 0.0;
  for (i=0; i<m; i++) if (k0[i] != 0.0)
    p += k0[i] * (1-pf(c*c/(d+1-i),(double) d+1-i, n)) / area(d+1-i);
  return( (s==TWO_SIDED) ? 2*p : p );
}

double taild_uniform(c,k0,m,d,s,n)
double c, *k0, n;
int m, d, s;
{ int i;
  double p;
  p = 0.0;
  for (i=0; i<m; i++) if (k0[i] != 0.0)
    p += k0[i] * 2*c*dbeta(1-c*c,(n-d+i-1)/2.0,(d+1-i)/2.0,0) / area(d+1-i);
  return( (s==TWO_SIDED) ? 2*p : p );
}

double taild_gaussian(c,k0,m,d,s,n)
double c, *k0, n;
int m, d, s;
{ int i;
  double p;
  p = 0.0;
  for (i=0; i<m; i++) if (k0[i] != 0.0)
    p += k0[i] * 2*c*dchisq(c*c,(double) d+1-i,0) / area(d+1-i);
  return( (s==TWO_SIDED) ? 2*p : p );
}

double taild_tprocess(c,k0,m,d,s,n)
double c, *k0, n;
int m, d, s;
{ int i;
  double p;
  p = 0.0;
  for (i=0; i<m; i++) if (k0[i] != 0.0)
    p += k0[i] * 2*c*df(c*c/(d+1-i),(double) d+1-i, n,0) / ((d+1-i)*area(d+1-i));
  return( (s==TWO_SIDED) ? 2*p : p );
}

double tailp(c,k0,m,d,s,nu, process)
double c, *k0, nu;
int m, d, s, process;
{ switch(process)
  { case UNIF:  return(tailp_uniform(c,k0,m,d,s,nu));
    case GAUSS: return(tailp_gaussian(c,k0,m,d,s,nu));
    case TPROC: return(tailp_tprocess(c,k0,m,d,s,nu));
  }
  mut_printf("taild: unknown process.\n");
  return(0.0);
}

double taild(c,k0,m,d,s,nu, process)
double c, *k0, nu;
int m, d, s, process;
{ switch(process)
  { case UNIF:  return(taild_uniform(c,k0,m,d,s,nu));
    case GAUSS: return(taild_gaussian(c,k0,m,d,s,nu));
    case TPROC: return(taild_tprocess(c,k0,m,d,s,nu));
  }
  mut_printf("taild: unknown process.\n");
  return(0.0);
}

double critval(alpha,k0,m,d,s,nu,process)
double alpha, *k0, nu;
int m, d, s, process;
{ double c, cn, c0, c1, tp, td;
  int j, maxit;
  double (*tpf)(), (*tdf)();

  maxit = 20;
  if (m<0)
  { mut_printf("critval: no terms?\n");
    return(2.0);
  }
  if (m>d+1) m = d+1;
  if ((alpha<=0) | (alpha>=1))
  { mut_printf("critval: invalid alpha %8.5f\n",alpha);
    return(2.0);
  }
  if (alpha>0.5)
    mut_printf("critval: A mighty large tail probability alpha=%8.5f\n",alpha);
  if (m==0) { d = 0; k0[0] = 1; m = 1; }

  switch(process)
  { case UNIF:
      c = 0.5; c0 = 0.0; c1 = 1.0;
      tpf = tailp_uniform;
      tdf = taild_uniform;
      break;
    case GAUSS:
      c = 2.0; c0 = 0.0; c1 = 0.0;
      tpf = tailp_gaussian;
      tdf = taild_gaussian;
      break;
    case TPROC:
      c = 2.0; c0 = 0.0; c1 = 0.0;
      tpf = tailp_tprocess;
      tdf = taild_tprocess;
      break;
    default:
      mut_printf("critval: unknown process.\n");
      return(0.0);
  }

  for (j=0; j<maxit; j++)
  { tp = tpf(c,k0,m,d,s,nu)-alpha;
    td = tdf(c,k0,m,d,s,nu);
    if (tp>0) c0 = c;
    if (tp<0) c1 = c;
    cn = c + tp/td;
    if (cn<c0) cn = (c+c0)/2;
    if ((c1>0.0) && (cn>c1)) cn = (c+c1)/2;
    c = cn;
    if (fabs(tp/alpha)<1.0e-10) return(c);
  }
  return(c);
}
