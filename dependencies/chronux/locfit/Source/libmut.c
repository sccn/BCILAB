/*
 * Copyright 1996-2006 Catherine Loader.
 */

#include "mex.h"
/*
 * Copyright 1996-2006 Catherine Loader.
 */
#include <math.h>
#include "mut.h"

/* stirlerr(n) = log(n!) - log( sqrt(2*pi*n)*(n/e)^n ) */

#define S0 0.083333333333333333333       /* 1/12 */
#define S1 0.00277777777777777777778     /* 1/360 */
#define S2 0.00079365079365079365079365  /* 1/1260 */
#define S3 0.000595238095238095238095238 /* 1/1680 */
#define S4 0.0008417508417508417508417508 /* 1/1188 */

/*
  error for 0, 0.5, 1.0, 1.5, ..., 14.5, 15.0.
*/
static double sferr_halves[31] = {
0.0, /* n=0 - wrong, place holder only */
0.1534264097200273452913848,  /* 0.5 */
0.0810614667953272582196702,  /* 1.0 */
0.0548141210519176538961390,  /* 1.5 */
0.0413406959554092940938221,  /* 2.0 */
0.03316287351993628748511048, /* 2.5 */
0.02767792568499833914878929, /* 3.0 */
0.02374616365629749597132920, /* 3.5 */
0.02079067210376509311152277, /* 4.0 */
0.01848845053267318523077934, /* 4.5 */
0.01664469118982119216319487, /* 5.0 */
0.01513497322191737887351255, /* 5.5 */
0.01387612882307074799874573, /* 6.0 */
0.01281046524292022692424986, /* 6.5 */
0.01189670994589177009505572, /* 7.0 */
0.01110455975820691732662991, /* 7.5 */
0.010411265261972096497478567, /* 8.0 */
0.009799416126158803298389475, /* 8.5 */
0.009255462182712732917728637, /* 9.0 */
0.008768700134139385462952823, /* 9.5 */
0.008330563433362871256469318, /* 10.0 */
0.007934114564314020547248100, /* 10.5 */
0.007573675487951840794972024, /* 11.0 */
0.007244554301320383179543912, /* 11.5 */
0.006942840107209529865664152, /* 12.0 */
0.006665247032707682442354394, /* 12.5 */
0.006408994188004207068439631, /* 13.0 */
0.006171712263039457647532867, /* 13.5 */
0.005951370112758847735624416, /* 14.0 */
0.005746216513010115682023589, /* 14.5 */
0.005554733551962801371038690  /* 15.0 */
};

double stirlerr(n)
double n;
{ double nn;

  if (n<15.0)
  { nn = 2.0*n;
    if (nn==(int)nn) return(sferr_halves[(int)nn]);
    return(mut_lgamma(n+1.0) - (n+0.5)*log((double)n)+n - HF_LG_PIx2);
  }

  nn = (double)n;
  nn = nn*nn;
  if (n>500) return((S0-S1/nn)/n);
  if (n>80) return((S0-(S1-S2/nn)/nn)/n);
  if (n>35) return((S0-(S1-(S2-S3/nn)/nn)/nn)/n);
  return((S0-(S1-(S2-(S3-S4/nn)/nn)/nn)/nn)/n);
}

double bd0(x,np)
double x, np;
{ double ej, s, s1, v;
  int j;
  if (fabs(x-np)<0.1*(x+np))
  {
    s = (x-np)*(x-np)/(x+np);
    v = (x-np)/(x+np);
    ej = 2*x*v; v = v*v;
    for (j=1; ;++j)
    { ej *= v;
      s1 = s+ej/((j<<1)+1);
      if (s1==s) return(s1);
      s = s1;
    }
  }
  return(x*log(x/np)+np-x);
}

/*
  Raw binomial probability calculation.
  (1) This has both p and q arguments, when one may be represented
      more accurately than the other (in particular, in df()).
  (2) This should NOT check that inputs x and n are integers. This
      should be done in the calling function, where necessary.
  (3) Does not check for 0<=p<=1 and 0<=q<=1 or NaN's. Do this in
      the calling function.
*/
double dbinom_raw(x,n,p,q,give_log)
double x, n, p, q;
int give_log;
{ double f, lc;

  if (p==0.0) return((x==0) ? D_1 : D_0);
  if (q==0.0) return((x==n) ? D_1 : D_0);

  if (x==0)
  { lc = (p<0.1) ? -bd0(n,n*q) - n*p : n*log(q);
    return( DEXP(lc) );
  }

  if (x==n)
  { lc = (q<0.1) ? -bd0(n,n*p) - n*q : n*log(p);
    return( DEXP(lc) );
  }

  if ((x<0) | (x>n)) return( D_0 );

  lc = stirlerr(n) - stirlerr(x) - stirlerr(n-x)
         - bd0(x,n*p) - bd0(n-x,n*q);
  f = (PIx2*x*(n-x))/n;

  return( FEXP(f,lc) );
}

double dbinom(x,n,p,give_log)
int x, n;
double p;
int give_log;
{ 
  if ((p<0) | (p>1) | (n<0)) return(INVALID_PARAMS);
  if (x<0) return( D_0 );
  
  return( dbinom_raw((double)x,(double)n,p,1-p,give_log) );
}

/*
  Poisson probability  lb^x exp(-lb) / x!.
  I don't check that x is an integer, since other functions
  that call dpois_raw() (i.e. dgamma) may use a fractional
  x argument.
*/
double dpois_raw(x,lambda,give_log)
int give_log;
double x, lambda;
{
  if (lambda==0) return( (x==0) ? D_1 : D_0 );
  if (x==0) return( DEXP(-lambda) );
  if (x<0) return( D_0 );

  return(FEXP( PIx2*x, -stirlerr(x)-bd0(x,lambda) ));
}

double dpois(x,lambda,give_log)
int x, give_log;
double lambda;
{
  if (lambda<0) return(INVALID_PARAMS);
  if (x<0) return( D_0 );

  return( dpois_raw((double)x,lambda,give_log) );
}

double dbeta(x,a,b,give_log)
double x, a, b;
int give_log;
{ double f, p;

  if ((a<=0) | (b<=0)) return(INVALID_PARAMS);
  if ((x<=0) | (x>=1)) return(D_0);

  if (a<1)
  { if (b<1)                                    /* a<1, b<1 */
    { f = a*b/((a+b)*x*(1-x));
      p = dbinom_raw(a,a+b,x,1-x,give_log);
    }
    else                                        /* a<1, b>=1 */
    { f = a/x;
      p = dbinom_raw(a,a+b-1,x,1-x,give_log);
    }
  }
  else
  { if (b<1)                                    /* a>=1, b<1 */
    { f = b/(1-x);
      p = dbinom_raw(a-1,a+b-1,x,1-x,give_log);
    }
    else                                        /* a>=1, b>=1 */
    { f = a+b-1;
      p = dbinom_raw(a-1,(a-1)+(b-1),x,1-x,give_log);
    }
  }

  return( (give_log) ? p + log(f) : p*f );
}

/*
 *   To evaluate the F density, write it as a Binomial probability
 *   with p = x*m/(n+x*m). For m>=2, use the simplest conversion.
 *   For m<2, (m-2)/2<0 so the conversion will not work, and we must use
 *   a second conversion. Note the division by p; this seems unavoidable
 *   for m < 2, since the F density has a singularity as x (or p) -> 0.
 */
double df(x,m,n,give_log)
double x, m, n;
int give_log;
{ double p, q, f, dens;

  if ((m<=0) | (n<=0)) return(INVALID_PARAMS);
  if (x <= 0.0) return(D_0);

  f = 1.0/(n+x*m);
  q = n*f;
  p = x*m*f;

  if (m>=2)
  { f = m*q/2;
    dens = dbinom_raw((m-2)/2.0, (m+n-2)/2.0, p, q, give_log);
  }
  else
  { f = m*m*q / (2*p*(m+n));
    dens = dbinom_raw(m/2.0, (m+n)/2.0, p, q, give_log);
  }

  return((give_log) ? log(f)+dens : f*dens);
}

/*
 * Gamma density,
 *                  lb^r x^{r-1} exp(-lb*x)
 *      p(x;r,lb) = -----------------------
 *                          (r-1)!
 *
 * If USE_SCALE is defined below, the lb argument will be interpreted
 * as a scale parameter (i.e. replace lb by 1/lb above). Otherwise,
 * it is interpreted as a rate parameter, as above.
 */

/* #define USE_SCALE */

double dgamma(x,r,lambda,give_log)
int give_log;
double x, r, lambda;
{ double pr;

  if ((r<=0) | (lambda<0)) return(INVALID_PARAMS);
  if (x<=0.0) return( D_0 );

#ifdef USE_SCALE
  lambda = 1.0/lambda;
#endif

  if (r<1)
  { pr = dpois_raw(r,lambda*x,give_log);
    return( (give_log) ?  pr + log(r/x) : pr*r/x );
  }

  pr = dpois_raw(r-1.0,lambda*x,give_log);
  return( (give_log) ? pr + log(lambda) : lambda*pr);
}

double dchisq(x, df, give_log)
double x, df;
int give_log;
{
 return(dgamma(x, df/2.0,
  0.5
  ,give_log));
/*
#ifdef USE_SCALE
  2.0
#else
  0.5
#endif
  ,give_log));
*/
}

/*
 * Given a sequence of r successes and b failures, we sample n (\le b+r)
 * items without replacement. The hypergeometric probability is the
 * probability of x successes:
 *
 *                dbinom(x,r,p) * dbinom(n-x,b,p)
 *   p(x;r,b,n) = ---------------------------------
 *                          dbinom(n,r+b,p)
 *
 * for any p. For numerical stability, we take p=n/(r+b); with this choice,
 * the denominator is not exponentially small.
 */
double dhyper(x,r,b,n,give_log)
int x, r, b, n, give_log;
{ double p, q, p1, p2, p3;

  if ((r<0) | (b<0) | (n<0) | (n>r+b))
    return( INVALID_PARAMS );

  if (x<0) return(D_0);

  if (n==0) return((x==0) ? D_1 : D_0);

  p = ((double)n)/((double)(r+b));
  q = ((double)(r+b-n))/((double)(r+b));

  p1 = dbinom_raw((double)x,(double)r,p,q,give_log);
  p2 = dbinom_raw((double)(n-x),(double)b,p,q,give_log);
  p3 = dbinom_raw((double)n,(double)(r+b),p,q,give_log);

  return( (give_log) ? p1 + p2 - p3 : p1*p2/p3 );
}

/*
  probability of x failures before the nth success.
*/
double dnbinom(x,n,p,give_log)
double n, p;
int x, give_log;
{ double prob, f;

  if ((p<0) | (p>1) | (n<=0)) return(INVALID_PARAMS);

  if (x<0) return( D_0 );

  prob = dbinom_raw(n,x+n,p,1-p,give_log);
  f = n/(n+x);

  return((give_log) ? log(f) + prob : f*prob);
}

double dt(x, df, give_log)
double x, df;
int give_log;
{ double t, u, f;

  if (df<=0.0) return(INVALID_PARAMS);

  /*
     exp(t) = Gamma((df+1)/2) /{ sqrt(df/2) * Gamma(df/2) }
            = sqrt(df/2) / ((df+1)/2) * Gamma((df+3)/2) / Gamma((df+2)/2).
     This form leads to a computation that should be stable for all
     values of df, including df -> 0 and df -> infinity.
  */
  t = -bd0(df/2.0,(df+1)/2.0) + stirlerr((df+1)/2.0) - stirlerr(df/2.0);

  if (x*x>df)
    u = log( 1+ x*x/df ) * df/2;
  else
    u = -bd0(df/2.0,(df+x*x)/2.0) + x*x/2.0;

  f = PIx2*(1+x*x/df);

  return( FEXP(f,t-u) );
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
/*
 * Provides mut_erf() and mut_erfc() functions. Also mut_pnorm().
 * Had too many problems with erf()'s built into math libraries
 * (when they existed). Hence the need to write my own...
 *
 * Algorithm from Walter Kr\"{a}mer, Frithjof Blomquist.
 * "Algorithms with Guaranteed Error Bounds for the Error Function
 * and Complementary Error Function"
 * Preprint 2000/2, Bergische Universt\"{a}t GH Wuppertal
 * http://www.math.uni-wuppertal.de/wrswt/preprints/prep_00_2.pdf
 *
 * Coded by Catherine Loader, September 2006.
 */

#include "mut.h"

double erf1(double x)  /* erf; 0 < x < 0.65) */
{ double p[5] = {1.12837916709551256e0,  /* 2/sqrt(pi) */
                 1.35894887627277916e-1,
                 4.03259488531795274e-2,
                 1.20339380863079457e-3,
                 6.49254556481904354e-5};
  double q[5] = {1.00000000000000000e0,
                 4.53767041780002545e-1,
                 8.69936222615385890e-2,
                 8.49717371168693357e-3,
                 3.64915280629351082e-4};
  double x2, p4, q4;
  x2 = x*x;
  p4 = p[0] + p[1]*x2 + p[2]*x2*x2 + p[3]*x2*x2*x2 + p[4]*x2*x2*x2*x2;
  q4 = q[0] + q[1]*x2 + q[2]*x2*x2 + q[3]*x2*x2*x2 + q[4]*x2*x2*x2*x2;
  return(x*p4/q4);
}

double erf2(double x) /* erfc; 0.65 <= x < 2.2 */
{ double p[6] = {9.99999992049799098e-1,
                 1.33154163936765307e0,
                 8.78115804155881782e-1,
                 3.31899559578213215e-1,
                 7.14193832506776067e-2,
                 7.06940843763253131e-3};
  double q[7] = {1.00000000000000000e0,
                 2.45992070144245533e0,
                 2.65383972869775752e0,
                 1.61876655543871376e0,
                 5.94651311286481502e-1,
                 1.26579413030177940e-1,
                 1.25304936549413393e-2};
  double x2, p5, q6;
  x2 = x*x;
  p5 = p[0] + p[1]*x + p[2]*x2 + p[3]*x2*x + p[4]*x2*x2 + p[5]*x2*x2*x;
  q6 = q[0] + q[1]*x + q[2]*x2 + q[3]*x2*x + q[4]*x2*x2 + q[5]*x2*x2*x + q[6]*x2*x2*x2;
  return( exp(-x2)*p5/q6 );
}

double erf3(double x) /* erfc; 2.2 < x <= 6 */
{ double p[6] = {9.99921140009714409e-1,
                 1.62356584489366647e0,
                 1.26739901455873222e0,
                 5.81528574177741135e-1,
                 1.57289620742838702e-1,
                 2.25716982919217555e-2};
  double q[7] = {1.00000000000000000e0,
                 2.75143870676376208e0,
                 3.37367334657284535e0,
                 2.38574194785344389e0,
                 1.05074004614827206e0,
                 2.78788439273628983e-1,
                 4.00072964526861362e-2};
  double x2, p5, q6;
  x2 = x*x;
  p5 = p[0] + p[1]*x + p[2]*x2 + p[3]*x2*x + p[4]*x2*x2 + p[5]*x2*x2*x;
  q6 = q[0] + q[1]*x + q[2]*x2 + q[3]*x2*x + q[4]*x2*x2 + q[5]*x2*x2*x + q[6]*x2*x2*x2;
  return( exp(-x2)*p5/q6 );
}

double erf4(double x)  /* erfc; x > 6.0 */
{ double p[5] = {5.64189583547756078e-1,
                 8.80253746105525775e0,
                 3.84683103716117320e1,
                 4.77209965874436377e1,
                 8.08040729052301677e0};
  double q[5] = {1.00000000000000000e0,
                 1.61020914205869003e1,
                 7.54843505665954743e1,
                 1.12123870801026015e2,
                 3.73997570145040850e1};
  double x2, p4, q4;
  if (x>26.5432) return(0.0);
  x2 = 1.0/(x*x);
  p4 = p[0] + p[1]*x2 + p[2]*x2*x2 + p[3]*x2*x2*x2 + p[4]*x2*x2*x2*x2;
  q4 = q[0] + q[1]*x2 + q[2]*x2*x2 + q[3]*x2*x2*x2 + q[4]*x2*x2*x2*x2;
  return(exp(-x*x)*p4/(x*q4));
}

double mut_erfc(double x)
{ if (x<0.0) return(2.0-mut_erfc(-x));
  if (x==0.0) return(1.0);
  if (x<0.65) return(1.0-erf1(x));
  if (x<2.2) return(erf2(x));
  if (x<6.0) return(erf3(x));
  return(erf4(x));
}

double mut_erf(double x)
{ 
  if (x<0.0) return(-mut_erf(-x));
  if (x==0.0) return(0.0);
  if (x<0.65) return(erf1(x));
  if (x<2.2) return(1.0-erf2(x));
  if (x<6.0) return(1.0-erf3(x));
  return(1.0-erf4(x));
}

double mut_pnorm(double x)
{ if (x<0.0) return(mut_erfc(-x/SQRT2)/2);
  return((1.0+mut_erf(x/SQRT2))/2);
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
#include "mut.h"

static double lookup_gamma[21] = {
  0.0, /* place filler */
  0.572364942924699971,  /* log(G(0.5)) = log(sqrt(pi)) */
  0.000000000000000000,  /* log(G(1)) = log(0!) */
 -0.120782237635245301,  /* log(G(3/2)) = log(sqrt(pi)/2)) */
  0.000000000000000000,  /* log(G(2)) = log(1!) */
  0.284682870472919181,  /* log(G(5/2)) = log(3sqrt(pi)/4) */
  0.693147180559945286,  /* log(G(3)) = log(2!) */
  1.200973602347074287,  /* etc */
  1.791759469228054957,
  2.453736570842442344,
  3.178053830347945752,
  3.957813967618716511,
  4.787491742782045812,
  5.662562059857141783,
  6.579251212010101213,
  7.534364236758732680,
  8.525161361065414667,
  9.549267257300996903,
 10.604602902745250859,
 11.689333420797268559,
 12.801827480081469091 };

/*
 * coefs are B(2n)/(2n(2n-1))  2n(2n-1) =
 *  2n  B(2n) 2n(2n-1) coef
 *   2  1/6     2      1/12
 *   4 -1/30   12     -1/360
 *   6  1/42   30      1/1260
 *   8 -1/30   56     -1/1680
 *  10  5/66   90      1/1188
 *  12 -691/2730 132   691/360360
 */

double mut_lgamma(double x)
{ double f, z, x2, se;
  int ix;

/* lookup table for common values.
 */
  ix = (int)(2*x);
  if (((ix>=1) & (ix<=20)) && (ix==2*x)) return(lookup_gamma[ix]);

  f = 1.0;
  while (x <= 15)
  { f *= x;
    x += 1.0;
  }

  x2 = 1.0/(x*x);
  z = (x-0.5)*log(x) - x + HF_LG_PIx2;
  se = (13860 - x2*(462 - x2*(132 - x2*(99 - 140*x2))))/(166320*x);

  return(z + se - log(f));
}

double mut_lgammai(int i)  /* log(Gamma(i/2)) for integer i */
{ if (i>20) return(mut_lgamma(i/2.0));
  return(lookup_gamma[i]);
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
/*
 * A is a n*p matrix, find the cholesky decomposition
 * of the first p rows. In most applications, will want n=p.
 *
 * chol_dec(A,n,p)       computes the decomoposition R'R=A.
 *                       (note that R is stored in the input A).
 * chol_solve(A,v,n,p)   computes (R'R)^{-1}v
 * chol_hsolve(A,v,n,p)  computes (R')^{-1}v
 * chol_isolve(A,v,n,p)  computes (R)^{-1}v
 * chol_qf(A,v,n,p)      computes ||(R')^{-1}v||^2.
 * chol_mult(A,v,n,p)    computes (R'R)v
 *
 * The solve functions assume that A is already decomposed.
 * chol_solve(A,v,n,p) is equivalent to applying chol_hsolve()
 * and chol_isolve() in sequence.
 */

#include <math.h>
#include "mut.h"

void chol_dec(A,n,p)
double *A;
int n, p;
{ int i, j, k;

  for (j=0; j<p; j++)
  { k = n*j+j;
    for (i=0; i<j; i++) A[k] -= A[n*j+i]*A[n*j+i];
    if (A[k]<=0)
    { for (i=j; i<p; i++) A[n*i+j] = 0.0; }
    else
    { A[k] = sqrt(A[k]);
      for (i=j+1; i<p; i++)
      { for (k=0; k<j; k++)
          A[n*i+j] -= A[n*i+k]*A[n*j+k];
        A[n*i+j] /= A[n*j+j];
      }
    }
  }
  for (j=0; j<p; j++)
    for (i=j+1; i<p; i++) A[n*j+i] = 0.0;
}

int chol_solve(A,v,n,p)
double *A, *v;
int n, p;
{ int i, j;

  for (i=0; i<p; i++)
  { for (j=0; j<i; j++) v[i] -= A[i*n+j]*v[j];
    v[i] /= A[i*n+i];
  }
  for (i=p-1; i>=0; i--)
  { for (j=i+1; j<p; j++) v[i] -= A[j*n+i]*v[j];
    v[i] /= A[i*n+i];
  }
  return(p);
}

int chol_hsolve(A,v,n,p)
double *A, *v;
int n, p;
{ int i, j;

  for (i=0; i<p; i++)
  { for (j=0; j<i; j++) v[i] -= A[i*n+j]*v[j];
    v[i] /= A[i*n+i];
  }
  return(p);
}

int chol_isolve(A,v,n,p)
double *A, *v;
int n, p;
{ int i, j;

  for (i=p-1; i>=0; i--)
  { for (j=i+1; j<p; j++) v[i] -= A[j*n+i]*v[j];
    v[i] /= A[i*n+i];
  }
  return(p);
}

double chol_qf(A,v,n,p)
double *A, *v;
int n, p;
{ int i, j;
  double sum;
 
  sum = 0.0;
  for (i=0; i<p; i++)
  { for (j=0; j<i; j++) v[i] -= A[i*n+j]*v[j];
    v[i] /= A[i*n+i];
    sum += v[i]*v[i];
  }
  return(sum);
}

int chol_mult(A,v,n,p)
double *A, *v;
int n, p;
{ int i, j;
  double sum;
  for (i=0; i<p; i++)
  { sum = 0;
    for (j=i; j<p; j++) sum += A[j*n+i]*v[j];
    v[i] = sum;
  }
  for (i=p-1; i>=0; i--)
  { sum = 0;
    for (j=0; j<=i; j++) sum += A[i*n+j]*v[j];
    v[i] = sum;
  }
  return(1);
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
#include <stdio.h>
#include <math.h>
#include "mut.h"
#define E_MAXIT 20
#define E_TOL 1.0e-10
#define SQR(x) ((x)*(x))

double e_tol(D,p)
double *D;
int p;
{ double mx;
  int i;
  if (E_TOL <= 0.0) return(0.0);
  mx = D[0];
  for (i=1; i<p; i++) if (D[i*(p+1)]>mx) mx = D[i*(p+1)];
  return(E_TOL*mx);
}

void eig_dec(X,P,d)
double *X, *P;
int d;
{ int i, j, k, iter, ms;
  double c, s, r, u, v;

  for (i=0; i<d; i++)
    for (j=0; j<d; j++) P[i*d+j] = (i==j);

  for (iter=0; iter<E_MAXIT; iter++)
  { ms = 0;
    for (i=0; i<d; i++)
      for (j=i+1; j<d; j++)
        if (SQR(X[i*d+j]) > 1.0e-15*fabs(X[i*d+i]*X[j*d+j]))
        { c = (X[j*d+j]-X[i*d+i])/2;
          s = -X[i*d+j];
          r = sqrt(c*c+s*s);
          c /= r;
          s = sqrt((1-c)/2)*(2*(s>0)-1);
          c = sqrt((1+c)/2);
          for (k=0; k<d; k++)
          { u = X[i*d+k]; v = X[j*d+k];
            X[i*d+k] = u*c+v*s;
            X[j*d+k] = v*c-u*s;
          }
          for (k=0; k<d; k++)
          { u = X[k*d+i]; v = X[k*d+j];
            X[k*d+i] = u*c+v*s;
            X[k*d+j] = v*c-u*s;
          }
          X[i*d+j] = X[j*d+i] = 0.0;
          for (k=0; k<d; k++)
          { u = P[k*d+i]; v = P[k*d+j];
            P[k*d+i] = u*c+v*s;
            P[k*d+j] = v*c-u*s;
          }
          ms = 1;
        }
    if (ms==0) return;
  }
  mut_printf("eig_dec not converged\n");
}

int eig_solve(J,x)
jacobian *J;
double *x;
{ int d, i, j, rank;
  double *D, *P, *Q, *w;
  double tol;

  D = J->Z;
  P = Q = J->Q;
  d = J->p;
  w = J->wk;

  tol = e_tol(D,d);

  rank = 0;
  for (i=0; i<d; i++)
  { w[i] = 0.0;
    for (j=0; j<d; j++) w[i] += P[j*d+i]*x[j];
  }
  for (i=0; i<d; i++)
    if (D[i*d+i]>tol)
    { w[i] /= D[i*(d+1)];
      rank++;
    }
  for (i=0; i<d; i++)
  { x[i] = 0.0;
    for (j=0; j<d; j++) x[i] += Q[i*d+j]*w[j];
  }
  return(rank);
}

int eig_hsolve(J,v)
jacobian *J;
double *v;
{ int i, j, p, rank;
  double *D, *Q, *w;
  double tol;

  D = J->Z;
  Q = J->Q;
  p = J->p;
  w = J->wk;

  tol = e_tol(D,p);
  rank = 0;

  for (i=0; i<p; i++)
  { w[i] = 0.0;
    for (j=0; j<p; j++) w[i] += Q[j*p+i]*v[j];
  }
  for (i=0; i<p; i++)
  { if (D[i*p+i]>tol)
    { v[i] = w[i]/sqrt(D[i*(p+1)]);
      rank++;
    }
    else v[i] = 0.0;
  }
  return(rank);
}

int eig_isolve(J,v)
jacobian *J;
double *v;
{ int i, j, p, rank;
  double *D, *Q, *w;
  double tol;

  D = J->Z;
  Q = J->Q;
  p = J->p;
  w = J->wk;

  tol = e_tol(D,p);
  rank = 0;

  for (i=0; i<p; i++)
  { if (D[i*p+i]>tol)
    { v[i] = w[i]/sqrt(D[i*(p+1)]);
      rank++;
    }
    else v[i] = 0.0;
  }

  for (i=0; i<p; i++)
  { w[i] = 0.0;
    for (j=0; j<p; j++) w[i] += Q[i*p+j]*v[j];
  }

  return(rank);
}

double eig_qf(J,v)
jacobian *J;
double *v;
{ int i, j, p;
  double sum, tol;

  p = J->p;
  sum = 0.0;
  tol = e_tol(J->Z,p);

  for (i=0; i<p; i++)
    if (J->Z[i*p+i]>tol)
    { J->wk[i] = 0.0;
      for (j=0; j<p; j++) J->wk[i] += J->Q[j*p+i]*v[j];
      sum += J->wk[i]*J->wk[i]/J->Z[i*p+i];
    }
  return(sum);
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
/*
 *  Integrate a function f over a circle or disc.
 */

#include "mut.h"

void setM(M,r,s,c,b)
double *M, r, s, c;
int b;
{ M[0] =-r*s; M[1] = r*c;
  M[2] = b*c; M[3] = b*s;
  M[4] =-r*c; M[5] = -s;
  M[6] = -s;  M[7] = 0.0;
  M[8] =-r*s; M[9] = c;
  M[10]=   c; M[11]= 0.0;
}

void integ_circ(f,r,orig,res,mint,b)
int (*f)(), mint, b;
double r, *orig, *res;
{ double y, x[2], theta, tres[MXRESULT], M[12], c, s;
  int i, j, nr;
  
  y = 0;
  for (i=0; i<mint; i++)
  { theta = 2*PI*(double)i/(double)mint;
    c = cos(theta); s = sin(theta);
    x[0] = orig[0]+r*c;
    x[1] = orig[1]+r*s;

    if (b!=0)
    { M[0] =-r*s; M[1] = r*c;
      M[2] = b*c; M[3] = b*s;
      M[4] =-r*c; M[5] = -s;
      M[6] = -s;  M[7] = 0.0;
      M[8] =-r*s; M[9] = c;
      M[10]=   c; M[11]= 0.0;
    }

    nr = f(x,2,tres,M);
    if (i==0) setzero(res,nr);
    for (j=0; j<nr; j++) res[j] += tres[j];
  }
  y = 2 * PI * ((b==0)?r:1.0) / mint;
  for (j=0; j<nr; j++) res[j] *= y;
}

void integ_disc(f,fb,fl,res,resb,mg)
int (*f)(), (*fb)(), *mg;
double *fl, *res, *resb;
{ double x[2], y, r, tres[MXRESULT], *orig, rmin, rmax, theta, c, s, M[12];
  int ct, ctb, i, j, k, nr, nrb, w;

  orig = &fl[2];
  rmax = fl[1];
  rmin = fl[0];
  y = 0.0;
  ct = ctb = 0;

  for (j=0; j<mg[1]; j++)
  { theta = 2*PI*(double)j/(double)mg[1];
    c = cos(theta); s = sin(theta);
    for (i= (rmin>0) ? 0 : 1; i<=mg[0]; i++)
    { r = rmin + (rmax-rmin)*i/mg[0];
      w = (2+2*(i&1)-(i==0)-(i==mg[0]));
      x[0] = orig[0] + r*c;
      x[1] = orig[1] + r*s;
      nr = f(x,2,tres,NULL);
      if (ct==0) setzero(res,nr);
      for (k=0; k<nr; k++) res[k] += w*r*tres[k];
      ct++;
      if (((i==0) | (i==mg[0])) && (fb!=NULL))
      { setM(M,r,s,c,1-2*(i==0));
        nrb = fb(x,2,tres,M);
        if (ctb==0) setzero(resb,nrb);
        ctb++;
        for (k=0; k<nrb; k++) resb[k] += tres[k];
      }
    }
  }


/*  for (i= (rmin>0) ? 0 : 1; i<=mg[0]; i++)
  {
    r = rmin + (rmax-rmin)*i/mg[0];
    w = (2+2*(i&1)-(i==0)-(i==mg[0]));

    for (j=0; j<mg[1]; j++)
    { theta = 2*PI*(double)j/(double)mg[1];
      c = cos(theta); s = sin(theta);
      x[0] = orig[0] + r*c;
      x[1] = orig[1] + r*s;
      nr = f(x,2,tres,NULL);
      if (ct==0) setzero(res,nr);
      ct++;
      for (k=0; k<nr; k++) res[k] += w*r*tres[k];

      if (((i==0) | (i==mg[0])) && (fb!=NULL))
      { setM(M,r,s,c,1-2*(i==0));
        nrb = fb(x,2,tres,M);
        if (ctb==0) setzero(resb,nrb);
        ctb++;
        for (k=0; k<nrb; k++) resb[k] += tres[k];
      }
    }
  } */
  for (j=0; j<nr; j++) res[j] *= 2*PI*(rmax-rmin)/(3*mg[0]*mg[1]);
  for (j=0; j<nrb; j++) resb[j] *= 2*PI/mg[1];
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
/*
 *  Multivariate integration of a vector-valued function
 *  using Monte-Carlo method.
 *
 *  uses drand48() random number generator. Does not seed.
 */

#include <stdlib.h>
#include "mut.h"
extern void setzero();

static double M[(1+MXIDIM)*MXIDIM*MXIDIM];

void monte(f,ll,ur,d,res,n)
int (*f)(), d, n;
double *ll, *ur, *res;
{
  int i, j, nr;
#ifdef WINDOWS
  mut_printf("Sorry, Monte-Carlo Integration not enabled.\n");
  for (i=0; i<n; i++) res[i] = 0.0;
#else
  double z, x[MXIDIM], tres[MXRESULT];

srand48(234L);

  for (i=0; i<n; i++)
  { for (j=0; j<d; j++) x[j] = ll[j] + (ur[j]-ll[j])*drand48();
    nr = f(x,d,tres,NULL);
    if (i==0) setzero(res,nr);
    for (j=0; j<nr; j++) res[j] += tres[j];
  }

  z = 1;
  for (i=0; i<d; i++) z *= (ur[i]-ll[i]);
  for (i=0; i<nr; i++) res[i] *= z/n;
#endif
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
/*
 *  Multivariate integration of a vector-valued function
 *  using Simpson's rule.
 */

#include <math.h>
#include "mut.h"
extern void setzero();

static double M[(1+MXIDIM)*MXIDIM*MXIDIM];

/* third order corners */
void simp3(fd,x,d,resd,delta,wt,i0,i1,mg,ct,res2,index)
int (*fd)(), d, wt, i0, i1, *mg, ct, *index;
double *x, *resd, *delta, *res2;
{ int k, l, m, nrd;
  double zb;

  for (k=i1+1; k<d; k++) if ((index[k]==0) | (index[k]==mg[k]))
  {
    setzero(M,d*d);
    m = 0; zb = 1.0;
    for (l=0; l<d; l++)
      if ((l!=i0) & (l!=i1) & (l!=k))
      { M[m*d+l] = 1.0;
        m++;
        zb *= delta[l];
      }
    M[(d-3)*d+i0] = (index[i0]==0) ? -1 : 1;
    M[(d-2)*d+i1] = (index[i1]==0) ? -1 : 1;
    M[(d-1)*d+k] = (index[k]==0) ? -1 : 1;
    nrd = fd(x,d,res2,M);
    if ((ct==0) & (i0==0) & (i1==1) & (k==2)) setzero(resd,nrd);
    for (l=0; l<nrd; l++)
      resd[l] += wt*zb*res2[l];
  }
}

/* second order corners */
void simp2(fc,fd,x,d,resc,resd,delta,wt,i0,mg,ct,res2,index)
int (*fc)(), (*fd)(), d, wt, i0, *mg, ct, *index;
double *x, *resc, *resd, *delta, *res2;
{ int j, k, l, nrc;
  double zb;
  for (j=i0+1; j<d; j++) if ((index[j]==0) | (index[j]==mg[j]))
  { setzero(M,d*d);
    l = 0; zb = 1;
    for (k=0; k<d; k++) if ((k!=i0) & (k!=j))
    { M[l*d+k] = 1.0;
      l++;
      zb *= delta[k];
    }
    M[(d-2)*d+i0] = (index[i0]==0) ? -1 : 1;
    M[(d-1)*d+j] = (index[j]==0) ? -1 : 1;
    nrc = fc(x,d,res2,M);
    if ((ct==0) & (i0==0) & (j==1)) setzero(resc,nrc);
    for (k=0; k<nrc; k++) resc[k] += wt*zb*res2[k];
       
    if (fd!=NULL)
      simp3(fd,x,d,resd,delta,wt,i0,j,mg,ct,res2,index);
  }
}

/* first order boundary */
void simp1(fb,fc,fd,x,d,resb,resc,resd,delta,wt,mg,ct,res2,index)
int (*fb)(), (*fc)(), (*fd)(), d, wt, *mg, ct, *index;
double *x, *resb, *resc, *resd, *delta, *res2;
{ int i, j, k, nrb;
  double zb;
  for (i=0; i<d; i++) if ((index[i]==0) | (index[i]==mg[i]))
  { setzero(M,(1+d)*d*d);
    k = 0;
    for (j=0; j<d; j++) if (j!=i)
    { M[k*d+j] = 1;
      k++;
    }
    M[(d-1)*d+i] = (index[i]==0) ? -1 : 1;
    nrb = fb(x,d,res2,M);
    zb = 1;
    for (j=0; j<d; j++) if (i!=j) zb *= delta[j];
        if ((ct==0) && (i==0))
          for (j=0; j<nrb; j++) resb[j] = 0.0;
    for (j=0; j<nrb; j++) resb[j] += wt*zb*res2[j];

    if (fc!=NULL)
      simp2(fc,fd,x,d,resc,resd,delta,wt,i,mg,ct,res2,index);
  }
}

void simpson4(f,fb,fc,fd,ll,ur,d,res,resb,resc,resd,mg,res2)
int (*f)(), (*fb)(), (*fc)(), (*fd)(), d, *mg;
double *ll, *ur, *res, *resb, *resc, *resd, *res2;
{ int ct, i, j, nr, wt, index[MXIDIM];
  double x[MXIDIM], delta[MXIDIM], z;

  for (i=0; i<d; i++)
  { index[i] = 0;
    x[i] = ll[i];
    if (mg[i]&1) mg[i]++;
    delta[i] = (ur[i]-ll[i])/(3*mg[i]);
  }
  ct = 0;

  while(1)
  { wt = 1;
    for (i=0; i<d; i++)
      wt *= (4-2*(index[i]%2==0)-(index[i]==0)-(index[i]==mg[i]));
    nr = f(x,d,res2,NULL);
    if (ct==0) setzero(res,nr);
    for (i=0; i<nr; i++) res[i] += wt*res2[i];

    if (fb!=NULL)
      simp1(fb,fc,fd,x,d,resb,resc,resd,delta,wt,mg,ct,res2,index);

    /* compute next grid point */
    for (i=0; i<d; i++)
    { index[i]++;
      if (index[i]>mg[i])
      { index[i] = 0;
        x[i] = ll[i];
        if (i==d-1) /* done */
        { z = 1.0;
          for (j=0; j<d; j++) z *= delta[j];
          for (j=0; j<nr; j++) res[j] *= z;
          return;
        }
      }
      else
      { x[i] = ll[i] + 3*delta[i]*index[i];
        i = d;
      }
    }
    ct++;
  }
}

void simpsonm(f,ll,ur,d,res,mg,res2)
int (*f)(), d, *mg;
double *ll, *ur, *res, *res2;
{ simpson4(f,NULL,NULL,NULL,ll,ur,d,res,NULL,NULL,NULL,mg,res2);
}

double simpson(f,l0,l1,m)
double (*f)(), l0, l1;
int m;
{ double x, sum;
  int i;
  sum = 0;
  for (i=0; i<=m; i++)
  { x = ((m-i)*l0 + i*l1)/m;
    sum += (2+2*(i&1)-(i==0)-(i==m)) * f(x);
  }
  return( (l1-l0) * sum / (3*m) );
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
#include "mut.h"

static double *res, *resb, *orig, rmin, rmax;
static int ct0;

void sphM(M,r,u)
double *M, r, *u;
{ double h, u1[3], u2[3];

  /* set the orthogonal unit vectors. */
  h = sqrt(u[0]*u[0]+u[1]*u[1]);
  if (h<=0)
  { u1[0] = u2[1] = 1.0;
    u1[1] = u1[2] = u2[0] = u2[2] = 0.0;
  }
  else
  { u1[0] = u[1]/h; u1[1] = -u[0]/h; u1[2] = 0.0;
    u2[0] = u[2]*u[0]/h; u2[1] = u[2]*u[1]/h; u2[2] = -h;
  }

  /* parameterize the sphere as r(cos(t)cos(v)u + sin(t)u1 + cos(t)sin(v)u2).
   * first layer of M is (dx/dt, dx/dv, dx/dr) at t=v=0.
   */
  M[0] = r*u1[0]; M[1] = r*u1[1]; M[2] = r*u1[2];
  M[3] = r*u2[0]; M[4] = r*u2[1]; M[5] = r*u2[2];
  M[6] = u[0]; M[7] = u[1]; M[8] = u[2];

  /* next layers are second derivative matrix of components of x(r,t,v).
   * d^2x/dt^2 = d^2x/dv^2 = -ru;    d^2x/dtdv = 0;
   * d^2x/drdt = u1; d^2x/drdv = u2; d^2x/dr^2 = 0.
   */

  M[9] = M[13] = -r*u[0];
  M[11]= M[15] = u1[0];
  M[14]= M[16] = u2[0];
  M[10]= M[12] = M[17] = 0.0;

  M[18]= M[22] = -r*u[1];
  M[20]= M[24] = u1[1];
  M[23]= M[25] = u2[1];
  M[19]= M[21] = M[26] = 0.0;

  M[27]= M[31] = -r*u[1];
  M[29]= M[33] = u1[1];
  M[32]= M[34] = u2[1];
  M[28]= M[30] = M[35] = 0.0;

}

double ip3(a,b)
double *a, *b;
{ return(a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
}

void rn3(a)
double *a;
{ double s;
  s = sqrt(ip3(a,a));
  a[0] /= s; a[1] /= s; a[2] /= s;
}

double sptarea(a,b,c)
double *a, *b, *c;
{ double ea, eb, ec, yab, yac, ybc, sab, sac, sbc;
  double ab[3], ac[3], bc[3], x1[3], x2[3];

  ab[0] = a[0]-b[0]; ab[1] = a[1]-b[1]; ab[2] = a[2]-b[2];
  ac[0] = a[0]-c[0]; ac[1] = a[1]-c[1]; ac[2] = a[2]-c[2];
  bc[0] = b[0]-c[0]; bc[1] = b[1]-c[1]; bc[2] = b[2]-c[2];
 
  yab = ip3(ab,a); yac = ip3(ac,a); ybc = ip3(bc,b);

  x1[0] = ab[0] - yab*a[0]; x2[0] = ac[0] - yac*a[0];
  x1[1] = ab[1] - yab*a[1]; x2[1] = ac[1] - yac*a[1];
  x1[2] = ab[2] - yab*a[2]; x2[2] = ac[2] - yac*a[2];
  sab = ip3(x1,x1); sac = ip3(x2,x2);
  ea = acos(ip3(x1,x2)/sqrt(sab*sac));

  x1[0] = ab[0] + yab*b[0]; x2[0] = bc[0] - ybc*b[0];
  x1[1] = ab[1] + yab*b[1]; x2[1] = bc[1] - ybc*b[1];
  x1[2] = ab[2] + yab*b[2]; x2[2] = bc[2] - ybc*b[2];
  sbc = ip3(x2,x2);
  eb = acos(ip3(x1,x2)/sqrt(sab*sbc));

  x1[0] = ac[0] + yac*c[0]; x2[0] = bc[0] + ybc*c[0];
  x1[1] = ac[1] + yac*c[1]; x2[1] = bc[1] + ybc*c[1];
  x1[2] = ac[2] + yac*c[2]; x2[2] = bc[2] + ybc*c[2];
  ec = acos(ip3(x1,x2)/sqrt(sac*sbc));

/*
 * Euler's formula is a+b+c-PI, except I've cheated...
 * a=ea, c=ec, b=PI-eb, which is more stable.
 */
  return(ea+ec-eb);
}

void li(x,f,fb,mint,ar)
double *x, ar;
int (*f)(), (*fb)(), mint;
{ int i, j, nr, nrb, ct1, w;
  double u[3], r, M[36];
  double sres[MXRESULT], tres[MXRESULT];

/* divide mint by 2, and force to even (Simpson's rule...)
 * to make comparable with rectangular interpretation of mint
 */
  mint <<= 1;
  if (mint&1) mint++;

  ct1 = 0;
  for (i= (rmin==0) ? 1 : 0; i<=mint; i++)
  {
    r = rmin + (rmax-rmin)*i/mint;
    w = 2+2*(i&1)-(i==0)-(i==mint);
    u[0] = orig[0]+x[0]*r;
    u[1] = orig[1]+x[1]*r;
    u[2] = orig[2]+x[2]*r;
    nr = f(u,3,tres,NULL);
    if (ct1==0) setzero(sres,nr);
    for (j=0; j<nr; j++)
      sres[j] += w*r*r*tres[j];
    ct1++;

    if ((fb!=NULL) && (i==mint)) /* boundary */
    { sphM(M,rmax,x);
      nrb = fb(u,3,tres,M);
      if (ct0==0) for (j=0; j<nrb; j++) resb[j] = 0.0;
      for (j=0; j<nrb; j++)
        resb[j] += tres[j]*ar;
    }
  }

  if (ct0==0) for (j=0; j<nr; j++) res[j] = 0.0;
  ct0++;

  for (j=0; j<nr; j++)
    res[j] += sres[j] * ar * (rmax-rmin)/(3*mint);
}

void sphint(f,fb,a,b,c,lev,mint,cent)
double *a, *b, *c;
int (*f)(), (*fb)(), lev, mint, cent;
{ double x[3], ab[3], ac[3], bc[3], ar;
  int i;

  if (lev>1)
  { ab[0] = a[0]+b[0]; ab[1] = a[1]+b[1]; ab[2] = a[2]+b[2]; rn3(ab);
    ac[0] = a[0]+c[0]; ac[1] = a[1]+c[1]; ac[2] = a[2]+c[2]; rn3(ac);
    bc[0] = b[0]+c[0]; bc[1] = b[1]+c[1]; bc[2] = b[2]+c[2]; rn3(bc);
    lev >>= 1;
    if (cent==0)
    { sphint(f,fb,a,ab,ac,lev,mint,1);
      sphint(f,fb,ab,bc,ac,lev,mint,0);
    }
    else
    { sphint(f,fb,a,ab,ac,lev,mint,1);
      sphint(f,fb,b,ab,bc,lev,mint,1);
      sphint(f,fb,c,ac,bc,lev,mint,1);
      sphint(f,fb,ab,bc,ac,lev,mint,1);
    }
    return;
  }

  x[0] = a[0]+b[0]+c[0];
  x[1] = a[1]+b[1]+c[1];
  x[2] = a[2]+b[2]+c[2];
  rn3(x);
  ar = sptarea(a,b,c);

  for (i=0; i<8; i++)
  { if (i>0)
    { x[0] = -x[0];
      if (i%2 == 0) x[1] = -x[1];
      if (i==4) x[2] = -x[2];
    }
    switch(cent)
    { case 2: /* the reflection and its 120', 240' rotations */
        ab[0] = x[0]; ab[1] = x[2]; ab[2] = x[1]; li(ab,f,fb,mint,ar);
        ab[0] = x[2]; ab[1] = x[1]; ab[2] = x[0]; li(ab,f,fb,mint,ar);
        ab[0] = x[1]; ab[1] = x[0]; ab[2] = x[2]; li(ab,f,fb,mint,ar);
      case 1: /* and the 120' and 240' rotations */
        ab[0] = x[1]; ab[1] = x[2]; ab[2] = x[0]; li(ab,f,fb,mint,ar);
        ac[0] = x[2]; ac[1] = x[0]; ac[2] = x[1]; li(ac,f,fb,mint,ar);
      case 0: /* and the triangle itself. */
        li( x,f,fb,mint,ar);
    }
  }
}

void integ_sphere(f,fb,fl,Res,Resb,mg)
double *fl, *Res, *Resb;
int (*f)(), (*fb)(), *mg;
{ double a[3], b[3], c[3];

  a[0] = 1; a[1] = a[2] = 0;
  b[1] = 1; b[0] = b[2] = 0;
  c[2] = 1; c[0] = c[1] = 0;
  
  res = Res;
  resb=Resb;
  orig = &fl[2];
  rmin = fl[0];
  rmax = fl[1];

  ct0 = 0;
  sphint(f,fb,a,b,c,mg[1],mg[0],0);
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
/*
 * solving symmetric equations using the jacobian structure. Currently, three
 * methods can be used: cholesky decomposition, eigenvalues, eigenvalues on
 * the correlation matrix.
 *
 * jacob_dec(J,meth)    decompose the matrix, meth=JAC_CHOL, JAC_EIG, JAC_EIGD
 * jacob_solve(J,v)     J^{-1}v
 * jacob_hsolve(J,v)    (R')^{-1/2}v
 * jacob_isolve(J,v)    (R)^{-1/2}v
 * jacob_qf(J,v)         v' J^{-1} v
 * jacob_mult(J,v)      (R'R) v   (pres. CHOL only)
 * where for each decomposition, R'R=J, although the different decomp's will
 * produce different R's.
 *
 * To set up the J matrix:
 * first, allocate storage: jac_alloc(J,p,wk)
 *   where p=dimension of matrix, wk is a numeric vector of length
 *   jac_reqd(p) (or NULL, to allocate automatically).
 * now, copy the numeric values to J->Z (numeric vector with length p*p).
 *   (or, just set J->Z to point to the data vector. But remember this
 *   will be overwritten by the decomposition).
 * finally, set:
 *    J->st=JAC_RAW;
 *    J->p = p;
 *
 * now, call jac_dec(J,meth) (optional) and the solve functions as required.
 *
 */

#include "math.h"
#include "mut.h"

#define DEF_METH JAC_EIGD

int jac_reqd(int p) { return(2*p*(p+1)); }

double *jac_alloc(J,p,wk)
jacobian *J;
int p;
double *wk;
{ if (wk==NULL)
    wk = (double *)calloc(2*p*(p+1),sizeof(double));
     if ( wk == NULL ) {
      printf("Problem allocating memory for wk\n");fflush(stdout);
    }
  J->Z = wk; wk += p*p;
  J->Q = wk; wk += p*p;
  J->wk= wk; wk += p;
  J->dg= wk; wk += p;
  return(wk);
}

void jacob_dec(J, meth)
jacobian *J;
int meth;
{ int i, j, p;

  if (J->st != JAC_RAW) return;

  J->sm = J->st = meth;
  switch(meth)
  { case JAC_EIG:
      eig_dec(J->Z,J->Q,J->p);
      return;
    case JAC_EIGD:
      p = J->p;
      for (i=0; i<p; i++)
        J->dg[i] = (J->Z[i*(p+1)]<=0) ? 0.0 : 1/sqrt(J->Z[i*(p+1)]);
      for (i=0; i<p; i++)
        for (j=0; j<p; j++)
          J->Z[i*p+j] *= J->dg[i]*J->dg[j];
      eig_dec(J->Z,J->Q,J->p);
      J->st = JAC_EIGD;
      return;
    case JAC_CHOL:
      chol_dec(J->Z,J->p,J->p);
      return;
    default: mut_printf("jacob_dec: unknown method %d",meth);
  }
}

int jacob_solve(J,v) /* (X^T W X)^{-1} v */
jacobian *J;
double *v;
{ int i, rank;

  if (J->st == JAC_RAW) jacob_dec(J,DEF_METH);

  switch(J->st)
  { case JAC_EIG:
      return(eig_solve(J,v));
    case JAC_EIGD:
      for (i=0; i<J->p; i++) v[i] *= J->dg[i];
      rank = eig_solve(J,v);
      for (i=0; i<J->p; i++) v[i] *= J->dg[i];
      return(rank);
    case JAC_CHOL:
      return(chol_solve(J->Z,v,J->p,J->p));
  }
  mut_printf("jacob_solve: unknown method %d",J->st);
  return(0);
}

int jacob_hsolve(J,v) /*  J^{-1/2} v */
jacobian *J;
double *v;
{ int i;

  if (J->st == JAC_RAW) jacob_dec(J,DEF_METH);

  switch(J->st)
  { case JAC_EIG:
      return(eig_hsolve(J,v));
    case JAC_EIGD: /* eigenvalues on corr matrix */
      for (i=0; i<J->p; i++) v[i] *= J->dg[i];
      return(eig_hsolve(J,v));
    case JAC_CHOL:
      return(chol_hsolve(J->Z,v,J->p,J->p));
  }
  mut_printf("jacob_hsolve: unknown method %d\n",J->st);
  return(0);
}

int jacob_isolve(J,v) /*  J^{-1/2} v */
jacobian *J;
double *v;
{ int i, r;

  if (J->st == JAC_RAW) jacob_dec(J,DEF_METH);

  switch(J->st)
  { case JAC_EIG:
        return(eig_isolve(J,v));
    case JAC_EIGD: /* eigenvalues on corr matrix */
        r = eig_isolve(J,v);
        for (i=0; i<J->p; i++) v[i] *= J->dg[i];
	return(r);
    case JAC_CHOL:
      return(chol_isolve(J->Z,v,J->p,J->p));
  }
  mut_printf("jacob_hsolve: unknown method %d\n",J->st);
  return(0);
}

double jacob_qf(J,v)  /* vT J^{-1} v */
jacobian *J;
double *v;
{ int i;

  if (J->st == JAC_RAW) jacob_dec(J,DEF_METH);

  switch (J->st)
  { case JAC_EIG:
      return(eig_qf(J,v));
    case JAC_EIGD:
      for (i=0; i<J->p; i++) v[i] *= J->dg[i];
      return(eig_qf(J,v));
    case JAC_CHOL:
      return(chol_qf(J->Z,v,J->p,J->p));
    default:
      mut_printf("jacob_qf: invalid method\n");
      return(0.0);
  }
}

int jacob_mult(J,v)  /* J v */
jacobian *J;
double *v;
{
  if (J->st == JAC_RAW) jacob_dec(J,DEF_METH);
  switch (J->st)
  { case JAC_CHOL:
      return(chol_mult(J->Z,v,J->p,J->p));
    default:
      mut_printf("jacob_mult: invalid method\n");
      return(0);
  }
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
/*
 *  Routines for maximization of a one dimensional function f()
 *    over an interval [xlo,xhi]. In all cases. the flag argument
 *    controls the return:
 *      flag='x', the maximizer xmax is returned.
 *                otherwise, maximum f(xmax) is returned.
 *
 *  max_grid(f,xlo,xhi,n,flag)
 *    grid maximization of f() over [xlo,xhi] with n intervals.
 *
 *  max_golden(f,xlo,xhi,n,tol,err,flag)
 *    golden section maximization.
 *    If n>2, an initial grid search is performed with n intervals
 *      (this helps deal with local maxima).
 *    convergence criterion is |x-xmax| < tol.
 *    err is an error flag.
 *    if flag='x', return value is xmax.
 *       otherwise, return value is f(xmax).
 *
 *  max_quad(f,xlo,xhi,n,tol,err,flag)
 *    quadratic maximization.
 *
 *  max_nr()
 *    newton-raphson, handles multivariate case.
 *
 *  TODO: additional error checking, non-convergence stop.
 */

#include <math.h>
#include "mut.h"

#define max_val(a,b) ((flag=='x') ? a : b)

double max_grid(f,xlo,xhi,n,flag)
double (*f)(), xlo, xhi;
int n;
char flag;
{ int i, mi;
  double x, y, mx, my;
  for (i=0; i<=n; i++)
  { x = xlo + (xhi-xlo)*i/n;
    y = f(x);
    if ((i==0) || (y>my))
    { mx = x;
      my = y;
      mi = i;
    }
  }
  if (mi==0) return(max_val(xlo,my));
  if (mi==n) return(max_val(xhi,my));
  return(max_val(mx,my));
}

double max_golden(f,xlo,xhi,n,tol,err,flag)
double (*f)(), xhi, xlo, tol;
int n, *err;
char flag;
{ double dlt, x0, x1, x2, x3, y0, y1, y2, y3;
  *err = 0;

  if (n>2)
  { dlt = (xhi-xlo)/n;
    x0 = max_grid(f,xlo,xhi,n,'x');
    if (xlo<x0) xlo = x0-dlt;
    if (xhi>x0) xhi = x0+dlt;
  }

  x0 = xlo; y0 = f(xlo);
  x3 = xhi; y3 = f(xhi);
  x1 = gold_rat*x0 + (1-gold_rat)*x3; y1 = f(x1);
  x2 = gold_rat*x3 + (1-gold_rat)*x0; y2 = f(x2);

  while (fabs(x3-x0)>tol)
  { if ((y1>=y0) && (y1>=y2))
    { x3 = x2; y3 = y2;
      x2 = x1; y2 = y1;
      x1 = gold_rat*x0 + (1-gold_rat)*x3; y1 = f(x1);
    }
    else if ((y2>=y3) && (y2>=y1))
    { x0 = x1; y0 = y1;
      x1 = x2; y1 = y2;
      x2 = gold_rat*x3 + (1-gold_rat)*x0; y2 = f(x2);
    }
    else
    { if (y3>y0) { x0 = x2; y0 = y2; }
            else { x3 = x1; y3 = y1; }
      x1 = gold_rat*x0 + (1-gold_rat)*x3; y1 = f(x1);
      x2 = gold_rat*x3 + (1-gold_rat)*x0; y2 = f(x2);
    }
  }
  if (y0>=y1) return(max_val(x0,y0));
  if (y3>=y2) return(max_val(x3,y3));
  return((y1>y2) ? max_val(x1,y1) : max_val(x2,y2));
}

double max_quad(f,xlo,xhi,n,tol,err,flag)
double (*f)(), xhi, xlo, tol;
int n, *err;
char flag;
{ double x0, x1, x2, xnew, y0, y1, y2, ynew, a, b;
  *err = 0;

  if (n>2)
  { x0 = max_grid(f,xlo,xhi,n,'x');
    if (xlo<x0) xlo = x0-1.0/n;
    if (xhi>x0) xhi = x0+1.0/n;
  }

  x0 = xlo; y0 = f(x0);
  x2 = xhi; y2 = f(x2);
  x1 = (x0+x2)/2; y1 = f(x1);

  while (x2-x0>tol)
  {
    /* first, check (y0,y1,y2) is a peak. If not,
     * next interval is the halve with larger of (y0,y2).
     */
    if ((y0>y1) | (y2>y1))
    { 
      if (y0>y2) { x2 = x1; y2 = y1; }
            else { x0 = x1; y0 = y1; }
      x1 = (x0+x2)/2;
      y1 = f(x1);
    }
    else /* peak */
    { a = (y1-y0)*(x2-x1) + (y1-y2)*(x1-x0);
      b = ((y1-y0)*(x2-x1)*(x2+x1) + (y1-y2)*(x1-x0)*(x1+x0))/2;
      /* quadratic maximizer is b/a. But first check if a's too
       * small, since we may be close to constant.
       */
      if ((a<=0) | (b<x0*a) | (b>x2*a))
      { /* split the larger halve */
        xnew = ((x2-x1) > (x1-x0)) ? (x1+x2)/2 : (x0+x1)/2;
      }
      else
      { xnew = b/a;
        if (10*xnew < (9*x0+x1)) xnew = (9*x0+x1)/10;
        if (10*xnew > (9*x2+x1)) xnew = (9*x2+x1)/10;
        if (fabs(xnew-x1) < 0.001*(x2-x0))
        {
          if ((x2-x1) > (x1-x0))
            xnew = (99*x1+x2)/100;
          else
            xnew = (99*x1+x0)/100;
        }
      }
      ynew = f(xnew);
      if (xnew>x1)
      { if (ynew >= y1) { x0 = x1; y0 = y1; x1 = xnew; y1 = ynew; }
                   else { x2 = xnew; y2 = ynew; }
      }
      else
      { if (ynew >= y1) { x2 = x1; y2 = y1; x1 = xnew; y1 = ynew; }
                   else { x0 = xnew; y0 = ynew; }
      }
    }
  }
  return(max_val(x1,y1));
}

double max_nr(F, coef, old_coef, f1, delta, J, p, maxit, tol, err)
double *coef, *old_coef, *f1, *delta, tol;
int (*F)(), p, maxit, *err;
jacobian *J;
{ double old_f, f, lambda;
  int i, j, fr;
  double nc, nd, cut;
  int rank;

  *err = NR_OK;
  J->p = p;
  fr = F(coef, &f, f1, J->Z); J->st = JAC_RAW;

  for (i=0; i<maxit; i++)
  { memcpy(old_coef,coef,p*sizeof(double));
    old_f = f;
    rank = jacob_solve(J,f1);
    memcpy(delta,f1,p*sizeof(double));

    if (rank==0) /* NR won't move! */
      delta[0] = -f/f1[0];

    lambda = 1.0;

    nc = innerprod(old_coef,old_coef,p);
    nd = innerprod(delta, delta, p);
    cut = sqrt(nc/nd);
    if (cut>1.0) cut = 1.0;
    cut *= 0.0001;
    do
    { for (j=0; j<p; j++) coef[j] = old_coef[j] + lambda*delta[j];
      f = old_f - 1.0;
      fr = F(coef, &f, f1, J->Z); J->st = JAC_RAW;
      if (fr==NR_BREAK) return(old_f);

      lambda = (fr==NR_REDUCE) ? lambda/2 : lambda/10.0;
    } while ((lambda>cut) & (f <= old_f - 1.0e-3));

    if (f < old_f - 1.0e-3)
    { *err = NR_NDIV;
      return(f);
    }
    if (fr==NR_REDUCE) return(f);

    if (fabs(f-old_f) < tol) return(f);

  }
  *err = NR_NCON;
  return(f);
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
#include <math.h>
#include "mut.h"

/* qr decomposition of X (n*p organized by column).
 * Take w for the ride, if not NULL.
 */
void qr(X,n,p,w)
double *X, *w;
int n, p;
{ int i, j, k, mi;
  double c, s, mx, nx, t;

  for (j=0; j<p; j++)
  { mi = j;
    mx = fabs(X[(n+1)*j]);
    nx = mx*mx;

    /* find the largest remaining element in j'th column, row mi.
     * flip that row with row j.
     */
    for (i=j+1; i<n; i++)
    { nx += X[j*n+i]*X[j*n+i];
      if (fabs(X[j*n+i])>mx)
      { mi = i;
        mx = fabs(X[j*n+i]);
      }
    }
    for (i=j; i<p; i++)
    { t = X[i*n+j];
      X[i*n+j] = X[i*n+mi];
      X[i*n+mi] = t;
    }
    if (w!=NULL) { t = w[j]; w[j] = w[mi]; w[mi] = t; }

    /* want the diag. element -ve, so we do the `good' Householder reflect.
     */
    if (X[(n+1)*j]>0)
    { for (i=j; i<p; i++) X[i*n+j] = -X[i*n+j];
      if (w!=NULL) w[j] = -w[j];
    }

    nx = sqrt(nx);
    c = nx*(nx-X[(n+1)*j]);
    if (c!=0)
    { for (i=j+1; i<p; i++)
      { s = 0;
        for (k=j; k<n; k++)
          s += X[i*n+k]*X[j*n+k];
        s = (s-nx*X[i*n+j])/c;
        for (k=j; k<n; k++)
          X[i*n+k] -= s*X[j*n+k];
        X[i*n+j] += s*nx;
      }
      if (w != NULL)
      { s = 0;
        for (k=j; k<n; k++)
          s += w[k]*X[n*j+k];
        s = (s-nx*w[j])/c;
        for (k=j; k<n; k++)
          w[k] -= s*X[n*j+k];
        w[j] += s*nx;
      }
      X[j*n+j] = nx;
    }
  }
}

void qrinvx(R,x,n,p)
double *R, *x;
int n, p;
{ int i, j;
  for (i=p-1; i>=0; i--)
  { for (j=i+1; j<p; j++) x[i] -= R[j*n+i]*x[j];
    x[i] /= R[i*n+i];
  }
}

void qrtinvx(R,x,n,p)
double *R, *x;
int n, p;
{ int i, j;
  for (i=0; i<p; i++)
  { for (j=0; j<i; j++) x[i] -= R[i*n+j]*x[j];
    x[i] /= R[i*n+i];
  }
}

void qrsolv(R,x,n,p)
double *R, *x;
int n, p;
{ qrtinvx(R,x,n,p);
  qrinvx(R,x,n,p);
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
/*
 *  solve f(x)=c by various methods, with varying stability etc...
 *    xlo and xhi should be initial bounds for the solution.
 *    convergence criterion is |f(x)-c| < tol.
 *
 *  double solve_bisect(f,c,xmin,xmax,tol,bd_flag,err)
 *  double solve_secant(f,c,xmin,xmax,tol,bd_flag,err)
 *    Bisection and secant methods for solving of f(x)=c.
 *    xmin and xmax are starting values and bound for solution.
 *    tol = convergence criterion, |f(x)-c| < tol.
 *    bd_flag = if (xmin,xmax) doesn't bound a solution, what action to take?
 *      BDF_NONE returns error.
 *      BDF_EXPRIGHT increases xmax.
 *      BDF_EXPLEFT  decreases xmin.
 *    err = error flag.
 *    The (xmin,xmax) bound is not formally necessary for the secant method.
 *    But having such a bound vastly improves stability; the code performs
 *    a bisection step whenever the iterations run outside the bounds.
 *
 *  double solve_nr(f,f1,c,x0,tol,err)
 *    Newton-Raphson solution of f(x)=c.
 *    f1 = f'(x).
 *    x0 = starting value.
 *    tol = convergence criteria, |f(x)-c| < tol.
 *    err = error flag.
 *    No stability checks at present.
 *
 *  double solve_fp(f,x0,tol)
 *    fixed-point iteration to solve f(x)=x.
 *    x0 = starting value.
 *    tol = convergence criteria, stops when |f(x)-x| < tol.
 *    Convergence requires |f'(x)|<1 in neighborhood of true solution;
 *      f'(x) \approx 0 gives the fastest convergence.
 *    No stability checks at present.
 *
 *  TODO: additional error checking, non-convergence stop.
 */

#include <math.h>
#include "mut.h"

typedef struct {
  double xmin, xmax, x0, x1;
  double ymin, ymax, y0, y1;
} solvest;

int step_expand(f,c,sv,bd_flag)
double (*f)(), c;
solvest *sv;
int bd_flag;
{ double x, y;
  if (sv->ymin*sv->ymax <= 0.0) return(0);
  if (bd_flag == BDF_NONE)
  { mut_printf("invalid bracket\n");
    return(1); /* error */
  }
  if (bd_flag == BDF_EXPRIGHT)
  { while (sv->ymin*sv->ymax > 0)
    { x = sv->xmax + 2*(sv->xmax-sv->xmin);
      y = f(x) - c;
      sv->xmin = sv->xmax; sv->xmax = x;
      sv->ymin = sv->ymax; sv->ymax = y;
    }
    return(0);
  }
  if (bd_flag == BDF_EXPLEFT)
  { while (sv->ymin*sv->ymax > 0)
    { x = sv->xmin - 2*(sv->xmax-sv->xmin);
      y = f(x) - c;
      sv->xmax = sv->xmin; sv->xmin = x;
      sv->ymax = sv->ymin; sv->ymin = y;
    }
    return(0);
  }
  mut_printf("step_expand: unknown bd_flag %d.\n",bd_flag);
  return(1);
}

int step_addin(sv,x,y)
solvest *sv;
double x, y;
{ sv->x1 = sv->x0; sv->x0 = x;
  sv->y1 = sv->y0; sv->y0 = y;
  if (y*sv->ymin > 0)
  { sv->xmin = x;
    sv->ymin = y;
    return(0);
  }
  if (y*sv->ymax > 0)
  { sv->xmax = x;
    sv->ymax = y;
    return(0);
  }  
  if (y==0)
  { sv->xmin = sv->xmax = x;
    sv->ymin = sv->ymax = 0;
    return(0);
  }
  return(1);
}

int step_bisect(f,c,sv)
double (*f)(), c;
solvest *sv;
{ double x, y;
  x = sv->x0 = (sv->xmin + sv->xmax)/2;
  y = sv->y0 = f(x)-c;
  return(step_addin(sv,x,y));
}

double solve_bisect(f,c,xmin,xmax,tol,bd_flag,err)
double (*f)(), c, xmin, xmax, tol;
int bd_flag, *err;
{ solvest sv;
  int z;
  *err = 0;
  sv.xmin = xmin; sv.ymin = f(xmin)-c;
  sv.xmax = xmax; sv.ymax = f(xmax)-c;
  *err = step_expand(f,c,&sv,bd_flag);
  if (*err>0) return(sv.xmin);
  while(1) /* infinite loop if f is discontinuous */
  { z = step_bisect(f,c,&sv);
    if (z)
    { *err = 1;
      return(sv.x0);
    }
    if (fabs(sv.y0)<tol) return(sv.x0);
  }
}

int step_secant(f,c,sv)
double (*f)(), c;
solvest *sv;
{ double x, y;
  if (sv->y0==sv->y1) return(step_bisect(f,c,sv));
  x = sv->x0 + (sv->x1-sv->x0)*sv->y0/(sv->y0-sv->y1);
  if ((x<=sv->xmin) | (x>=sv->xmax)) return(step_bisect(f,c,sv));
  y = f(x)-c;
  return(step_addin(sv,x,y));
}

double solve_secant(f,c,xmin,xmax,tol,bd_flag,err)
double (*f)(), c, xmin, xmax, tol;
int bd_flag, *err;
{ solvest sv;
  int z;
  *err = 0;
  sv.xmin = xmin; sv.ymin = f(xmin)-c;
  sv.xmax = xmax; sv.ymax = f(xmax)-c;
  *err = step_expand(f,c,&sv,bd_flag);
  if (*err>0) return(sv.xmin);
  sv.x0 = sv.xmin; sv.y0 = sv.ymin;
  sv.x1 = sv.xmax; sv.y1 = sv.ymax;
  while(1) /* infinite loop if f is discontinuous */
  { z = step_secant(f,c,&sv);
    if (z)
    { *err = 1;
      return(sv.x0);
    }
    if (fabs(sv.y0)<tol) return(sv.x0);
  }
}

double solve_nr(f,f1,c,x0,tol,err)
double (*f)(), (*f1)(), c, x0, tol;
int *err;
{ double y;
  do
  { y = f(x0)-c;
    x0 -= y/f1(x0);
  } while (fabs(y)>tol);
  return(x0);
}

double solve_fp(f,x0,tol,maxit)
double (*f)(), x0, tol;
int maxit;
{ double x1;
  int i;
  for (i=0; i<maxit; i++)
  { x1 = f(x0);
    if (fabs(x1-x0)<tol) return(x1);
    x0 = x1;
  }
  return(x1); /* although it hasn't converged */
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
#include "mut.h"

void svd(x,p,q,d,mxit)  /* svd of square matrix */
double *x, *p, *q;
int d, mxit;
{ int i, j, k, iter, ms, zer;
  double r, u, v, cp, cm, sp, sm, c1, c2, s1, s2, mx;
  for (i=0; i<d; i++)
    for (j=0; j<d; j++) p[i*d+j] = q[i*d+j] = (i==j);
  for (iter=0; iter<mxit; iter++)
  { ms = 0;
    for (i=0; i<d; i++)
      for (j=i+1; j<d; j++)
      { s1 = fabs(x[i*d+j]);
        s2 = fabs(x[j*d+i]);
        mx = (s1>s2) ? s1 : s2;
        zer = 1;
        if (mx*mx>1.0e-15*fabs(x[i*d+i]*x[j*d+j]))
        { if (fabs(x[i*(d+1)])<fabs(x[j*(d+1)]))
          { for (k=0; k<d; k++)
            { u = x[i*d+k]; x[i*d+k] = x[j*d+k]; x[j*d+k] = u;
              u = p[k*d+i]; p[k*d+i] = p[k*d+j]; p[k*d+j] = u;
            }
            for (k=0; k<d; k++)
            { u = x[k*d+i]; x[k*d+i] = x[k*d+j]; x[k*d+j] = u;
              u = q[k*d+i]; q[k*d+i] = q[k*d+j]; q[k*d+j] = u;
            }
          }
          cp = x[i*(d+1)]+x[j*(d+1)];
          sp = x[j*d+i]-x[i*d+j];
          r = sqrt(cp*cp+sp*sp);
          if (r>0) { cp /= r; sp /= r; }
              else { cp = 1.0; zer = 0;}
          cm = x[i*(d+1)]-x[j*(d+1)];
          sm = x[i*d+j]+x[j*d+i];
          r = sqrt(cm*cm+sm*sm);
          if (r>0) { cm /= r; sm /= r; }
              else { cm = 1.0; zer = 0;}
          c1 = cm+cp;
          s1 = sm+sp;
          r = sqrt(c1*c1+s1*s1);
          if (r>0) { c1 /= r; s1 /= r; }
              else { c1 = 1.0; zer = 0;}
          if (fabs(s1)>ms) ms = fabs(s1);
          c2 = cm+cp;
          s2 = sp-sm;
          r = sqrt(c2*c2+s2*s2);
          if (r>0) { c2 /= r; s2 /= r; }
              else { c2 = 1.0; zer = 0;}
          for (k=0; k<d; k++)
          { u = x[i*d+k]; v = x[j*d+k];
            x[i*d+k] = c1*u+s1*v;
            x[j*d+k] = c1*v-s1*u;
            u = p[k*d+i]; v = p[k*d+j];
            p[k*d+i] = c1*u+s1*v;
            p[k*d+j] = c1*v-s1*u;
          }
          for (k=0; k<d; k++)
          { u = x[k*d+i]; v = x[k*d+j];
            x[k*d+i] = c2*u-s2*v;
            x[k*d+j] = s2*u+c2*v;
            u = q[k*d+i]; v = q[k*d+j];
            q[k*d+i] = c2*u-s2*v;
            q[k*d+j] = s2*u+c2*v;
          }
          if (zer) x[i*d+j] = x[j*d+i] = 0.0;
          ms = 1;
        }
      }
    if (ms==0) iter=mxit+10;
  }
  if (iter==mxit) mut_printf("Warning: svd not converged.\n");
  for (i=0; i<d; i++)
    if (x[i*d+i]<0)
    { x[i*d+i] = -x[i*d+i];
      for (j=0; j<d; j++) p[j*d+i] = -p[j*d+i];
    }
}

int svdsolve(x,w,P,D,Q,d,tol) /* original X = PDQ^T; comp. QD^{-1}P^T x */
double *x, *w, *P, *D, *Q, tol;
int d;
{ int i, j, rank;
  double mx;
  if (tol>0)
  { mx = D[0];
    for (i=1; i<d; i++) if (D[i*(d+1)]>mx) mx = D[i*(d+1)];
    tol *= mx;
  }
  rank = 0;
  for (i=0; i<d; i++)
  { w[i] = 0.0;
    for (j=0; j<d; j++) w[i] += P[j*d+i]*x[j];
  }
  for (i=0; i<d; i++)
    if (D[i*d+i]>tol)
    { w[i] /= D[i*(d+1)];
      rank++;
    }
  for (i=0; i<d; i++)
  { x[i] = 0.0;
    for (j=0; j<d; j++) x[i] += Q[i*d+j]*w[j];
  }
  return(rank);
}

void hsvdsolve(x,w,P,D,Q,d,tol) /* original X = PDQ^T; comp. D^{-1/2}P^T x */
double *x, *w, *P, *D, *Q, tol;
int d;
{ int i, j;
  double mx;
  if (tol>0)
  { mx = D[0];
    for (i=1; i<d; i++) if (D[i*(d+1)]>mx) mx = D[i*(d+1)];
    tol *= mx;
  }
  for (i=0; i<d; i++)
  { w[i] = 0.0;
    for (j=0; j<d; j++) w[i] += P[j*d+i]*x[j];
  }
  for (i=0; i<d; i++) if (D[i*d+i]>tol) w[i] /= sqrt(D[i*(d+1)]);
  for (i=0; i<d; i++) x[i] = w[i];
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
/*
 *   Includes some miscellaneous vector functions:
 *     setzero(v,p)           sets all elements of v to 0.
 *     unitvec(x,k,p)         sets x to k'th unit vector e_k.
 *     innerprod(v1,v2,p)     inner product.
 *     addouter(A,v1,v2,p,c)  A <- A + c * v_1 v2^T
 *     multmatscal(A,z,n)     A <- A*z
 *     matrixmultiply(A,B,C,m,n,p)    C(m*p) <- A(m*n) * B(n*p)
 *     transpose(x,m,n)       inline transpose
 *     m_trace(x,n)           trace
 *     vecsum(x,n)            sum elements.
 */

#include "mut.h"

void setzero(v,p)
double *v;
int p;
{ int i;
  for (i=0; i<p; i++) v[i] = 0.0;
}

void unitvec(x,k,p)
double *x;
int k, p;
{ setzero(x,p);
  x[k] = 1.0;
}

double innerprod(v1,v2,p)
double *v1, *v2;
int p;
{ int i;
  double s;
  s = 0;
  for (i=0; i<p; i++) s += v1[i]*v2[i];
  return(s);
}

void addouter(A,v1,v2,p,c)
double *A, *v1, *v2, c;
int p;
{ int i, j;
  for (i=0; i<p; i++)
    for (j=0; j<p; j++)
      A[i*p+j] += c*v1[i]*v2[j];
}

void multmatscal(A,z,n)
double *A, z;
int n;
{ int i;
  for (i=0; i<n; i++) A[i] *= z;
}

/* matrix multiply A (m*n) times B (n*p).
 * store in C (m*p).
 * all matrices stored by column.
 */
void matrixmultiply(A,B,C,m,n,p)
double *A, *B, *C;
int m, n, p;
{ int i, j, k, ij;
  for (i=0; i<m; i++)
    for (j=0; j<p; j++)
    { ij = j*m+i;
      C[ij] = 0.0;
      for (k=0; k<n; k++)
        C[ij] += A[k*m+i] * B[j*n+k];
    }
}

/*
 *  transpose() transposes an m*n matrix in place.
 *  At input, the matrix has n rows, m columns and
 *    x[0..n-1] is the is the first column.
 *  At output, the matrix has m rows, n columns and
 *    x[0..m-1] is the first column.
 */
void transpose(x,m,n)
double *x;
int m, n;
{ int t0, t, ti, tj;
  double z;
  for (t0=1; t0<m*n-2; t0++)
    { ti = t0%m; tj = t0/m;
      do
      { t = ti*n+tj;
        ti= t%m;
        tj= t/m;
      } while (t<t0);
      z = x[t];
      x[t] = x[t0];
      x[t0] = z;
    }
}

/* trace of an n*n square matrix. */
double m_trace(x,n)
double *x;
int n;
{ int i;
  double sum;
  sum = 0;
  for (i=0; i<n; i++)
    sum += x[i*(n+1)];
  return(sum);
}

double vecsum(v,n)
double *v;
int n;
{ int i;
  double sum;
  sum = 0.0;
  for (i=0; i<n; i++) sum += v[i];
  return(sum);
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
/*
  miscellaneous functions that may not be defined in the math
  libraries. The implementations are crude.
  mut_daws(x)   -- dawson's function
  mut_exp(x)    -- exp(x), but it won't overflow.

  where required, these must be #define'd in header files.

  also includes
  ptail(x)    -- exp(x*x/2)*int_{-\infty}^x exp(-u^2/2)du for x < -6.
  logit(x)    -- logistic function.
  expit(x)    -- inverse of logit.
  factorial(n)-- factorial
 */

#include "mut.h"

double mut_exp(x)
double x;
{ if (x>700.0) return(1.014232054735004e+304);
  return(exp(x));
}

double mut_daws(x)
double x;
{ static double val[] = {
  0,  0.24485619356002, 0.46034428261948, 0.62399959848185, 0.72477845900708,
      0.76388186132749, 0.75213621001998, 0.70541701910853, 0.63998807456541,
      0.56917098836654, 0.50187821196415, 0.44274283060424, 0.39316687916687,
      0.35260646480842, 0.31964847250685, 0.29271122077502, 0.27039629581340,
      0.25160207761769, 0.23551176224443, 0.22153505358518, 0.20924575719548,
      0.19833146819662, 0.18855782729305, 0.17974461154688, 0.17175005072385 };
  double h, f0, f1, f2, y, z, xx;
  int j, m;
  if (x<0) return(-mut_daws(-x));
  if (x>6)
  { /* Tail series: 1/x + 1/x^3 + 1.3/x^5 + 1.3.5/x^7 + ...  */
    y = z = 1/x;
    j = 0;
    while (((f0=(2*j+1)/(x*x))<1) && (y>1.0e-10*z))
    { y *= f0;
      z += y;
      j++;
    }
    return(z);
  }
  m = (int) (4*x);
  h = x-0.25*m;
  if (h>0.125)
  { m++;
    h = h-0.25;
  }
  xx = 0.25*m;
  f0 = val[m];
  f1 = 1-xx*f0;
  z = f0+h*f1;
  y = h;
  j = 2;
  while (fabs(y)>z*1.0e-10)
  { f2 = -(j-1)*f0-xx*f1;
    y *= h/j;
    z += y*f2;
    f0 = f1; f1 = f2;
    j++;
  }
  return(z);
}

double ptail(x) /* exp(x*x/2)*int_{-\infty}^x exp(-u^2/2)du for x < -6 */
double x;
{ double y, z, f0;
  int j;
  y = z = -1.0/x;
  j = 0;
  while ((fabs(f0= -(2*j+1)/(x*x))<1) && (fabs(y)>1.0e-10*z))
  { y *= f0;
    z += y;
    j++;
  }
  return(z);
}

double logit(x)
double x;
{ return(log(x/(1-x)));
}

double expit(x)
double x;
{ double u;
  if (x<0)
  { u = exp(x);
    return(u/(1+u));
  }
  return(1/(1+exp(-x)));
}

int factorial(n)
int n;
{ if (n<=1) return(1.0);
  return(n*factorial(n-1));
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
/*
 * Constrained maximization of a bivariate function.
 *   maxbvgrid(f,x,ll,ur,m0,m1)
 *     maximizes over a grid of m0*m1 points. Returns the maximum,
 *     and the maximizer through the array x. Usually this is a starter,
 *     to choose between local maxima, followed by other routines to refine.
 *
 *  maxbvstep(f,x,ymax,h,ll,ur,err)
 *     essentially multivariate bisection. A 3x3 grid of points is
 *     built around the starting value (x,ymax). This grid is moved
 *     around (step size h[0] and h[1] in the two dimensions) until
 *     the maximum is in the middle.  Then, the step size is halved.
 *     Usually, this will be called in a loop.
 *     The error flag is set if the maximum can't be centered in a
 *     reasonable number of steps.
 *
 *  maxbv(f,x,h,ll,ur,m0,m1,tol)
 *     combines the two previous functions. It begins with a grid search
 *     (if m0>0 and m1>0), followed by refinement. Refines until both h
 *     components are < tol.
 */
#include "mut.h"

#define max(a,b) ((a)>(b) ? (a) : (b))
#define min(a,b) ((a)<(b) ? (a) : (b))

double maxbvgrid(f,x,ll,ur,m0,m1,con)
double (*f)(), *x, *ll, *ur;
int m0, m1, *con;
{ int i, j, im, jm;
  double y, ymax;

  im = -1;
  for (i=0; i<=m0; i++)
  { x[0] = ((m0-i)*ll[0] + i*ur[0])/m0;
    for (j=0; j<=m1; j++)
    { x[1] = ((m1-j)*ll[1] + j*ur[1])/m1;
      y = f(x);
      if ((im==-1) || (y>ymax))
      { im = i; jm = j;
        ymax = y;
      }
    }
  }

  x[0] = ((m0-im)*ll[0] + im*ur[0])/m0;
  x[1] = ((m1-jm)*ll[1] + jm*ur[1])/m1;
  con[0] = (im==m0)-(im==0);
  con[1] = (jm==m1)-(jm==0);
  return(ymax);
}

double maxbvstep(f,x,ymax,h,ll,ur,err,con)
double (*f)(), *x, ymax, *h, *ll, *ur;
int *err, *con;
{ int i, j, ij, imax, steps, cts[2];
  double newx, X[9][2], y[9];

  imax =4; y[4] = ymax;

  for (i=(con[0]==-1)-1; i<2-(con[0]==1); i++)
    for (j=(con[1]==-1)-1; j<2-(con[1]==1); j++)
    { ij = 3*i+j+4;
      X[ij][0] = x[0]+i*h[0];
      if (X[ij][0] < ll[0]+0.001*h[0]) X[ij][0] = ll[0];
      if (X[ij][0] > ur[0]-0.001*h[0]) X[ij][0] = ur[0];
      X[ij][1] = x[1]+j*h[1];
      if (X[ij][1] < ll[1]+0.001*h[1]) X[ij][1] = ll[1];
      if (X[ij][1] > ur[1]-0.001*h[1]) X[ij][1] = ur[1];
      if (ij != 4)
      { y[ij] = f(X[ij]);
        if (y[ij]>ymax) { imax = ij; ymax = y[ij]; }
      }
    }

  steps = 0;
  cts[0] = cts[1] = 0;
  while ((steps<20) && (imax != 4))
  { steps++;
    if ((con[0]>-1) && ((imax/3)==0)) /* shift right */
    {
      cts[0]--;
      for (i=8; i>2; i--)
      { X[i][0] = X[i-3][0]; y[i] = y[i-3];
      }
      imax = imax+3;
      if (X[imax][0]==ll[0])
        con[0] = -1;
      else
      { newx = X[imax][0]-h[0];
        if (newx < ll[0]+0.001*h[0]) newx = ll[0];
        for (i=(con[1]==-1); i<3-(con[1]==1); i++)
        { X[i][0] = newx;
          y[i] = f(X[i]);
          if (y[i]>ymax) { ymax = y[i]; imax = i; }
        }
        con[0] = 0;
      }
    }

    if ((con[0]<1) && ((imax/3)==2)) /* shift left */
    {
      cts[0]++;
      for (i=0; i<6; i++)
      { X[i][0] = X[i+3][0]; y[i] = y[i+3];
      }
      imax = imax-3;
      if (X[imax][0]==ur[0])
        con[0] = 1;
      else
      { newx = X[imax][0]+h[0];
        if (newx > ur[0]-0.001*h[0]) newx = ur[0];
        for (i=6+(con[1]==-1); i<9-(con[1]==1); i++)
        { X[i][0] = newx;
          y[i] = f(X[i]);
          if (y[i]>ymax) { ymax = y[i]; imax = i; }
        }
        con[0] = 0;
      }
    }

    if ((con[1]>-1) && ((imax%3)==0)) /* shift up */
    {
      cts[1]--;
      for (i=9; i>0; i--) if (i%3 > 0)
      { X[i][1] = X[i-1][1]; y[i] = y[i-1];
      }
      imax = imax+1;
      if (X[imax][1]==ll[1])
        con[1] = -1;
      else
      { newx = X[imax][1]-h[1];
        if (newx < ll[1]+0.001*h[1]) newx = ll[1];
        for (i=3*(con[0]==-1); i<7-(con[0]==1); i=i+3)
        { X[i][1] = newx;
          y[i] = f(X[i]);
          if (y[i]>ymax) { ymax = y[i]; imax = i; }
        }
        con[1] = 0;
      }
    }

    if ((con[1]<1) && ((imax%3)==2)) /* shift down */
    {
      cts[1]++;
      for (i=0; i<9; i++) if (i%3 < 2)
      { X[i][1] = X[i+1][1]; y[i] = y[i+1];
      }
      imax = imax-1;
      if (X[imax][1]==ur[1])
        con[1] = 1;
      else
      { newx = X[imax][1]+h[1];
        if (newx > ur[1]-0.001*h[1]) newx = ur[1];
        for (i=2+3*(con[0]==-1); i<9-(con[0]==1); i=i+3)
        { X[i][1] = newx;
          y[i] = f(X[i]);
          if (y[i]>ymax) { ymax = y[i]; imax = i; }
        }
        con[1] = 0;
      }
    }
/* if we've taken 3 steps in one direction, try increasing the
 * corresponding h.
 */
    if ((cts[0]==-2) | (cts[0]==2))
    { h[0] = 2*h[0]; cts[0] = 0; }
    if ((cts[1]==-2) | (cts[1]==2))
    { h[1] = 2*h[1]; cts[1] = 0; }
  }

  if (steps==40)
    *err = 1;
  else
  {
    h[0] /= 2.0; h[1] /= 2.0;
    *err = 0;
  }

  x[0] = X[imax][0];
  x[1] = X[imax][1];
  return(y[imax]);
}

#define BQMmaxp 5

int boxquadmin(J,b0,p,x0,ll,ur)
jacobian *J;
double *b0, *x0, *ll, *ur;
int p;
{ double b[BQMmaxp], x[BQMmaxp], L[BQMmaxp*BQMmaxp], C[BQMmaxp*BQMmaxp], d[BQMmaxp];
  double f, fmin;
  int i, imin, m, con[BQMmaxp], rlx;

  if (p>BQMmaxp) mut_printf("boxquadmin: maxp is 5.\n");
  if (J->st != JAC_RAW) mut_printf("boxquadmin: must start with JAC_RAW.\n");

  m = 0;
  setzero(L,p*p);
  setzero(x,p);
  memcpy(C,J->Z,p*p*sizeof(double));
  for (i=0; i<p; i++) con[i] = 0;

  do
  {
/* first, keep minimizing and add constraints, one at a time.
 */
    do
    {
      matrixmultiply(C,x,b,p,p,1);
      for (i=0; i<p; i++) b[i] += b0[i];
      conquadmin(J,b,p,L,d,m);
/* if C matrix is not pd, don't even bother.
 * this relies on having used cholesky dec.
 */
if ((J->Z[0]==0.0) | (J->Z[3]==0.0)) return(1);
      fmin = 1.0;
      for (i=0; i<p; i++) if (con[i]==0)
      { f = 1.0;
        if (x0[i]+x[i]+b[i] < ll[i]) f = (ll[i]-x[i]-x0[i])/b[i];
        if (x0[i]+x[i]+b[i] > ur[i]) f = (ur[i]-x[i]-x0[i])/b[i];
        if (f<fmin) fmin = f;
        imin = i;
      }
      for (i=0; i<p; i++) x[i] += fmin*b[i];
      if (fmin<1.0)
      { L[m*p+imin] = 1;
        m++;
        con[imin] = (b[imin]<0) ? -1 : 1;
      }
    } while ((fmin < 1.0) & (m<p));

/* now, can I relax any constraints?
 * compute slopes at current point. Can relax if:
 *   slope is -ve on a lower boundary.
 *   slope is +ve on an upper boundary.
 */
    rlx = 0;
    if (m>0)
    { matrixmultiply(C,x,b,p,p,1);
      for (i=0; i<p; i++) b[i] += b0[i];
      for (i=0; i<p; i++)
      { if ((con[i]==-1)&& (b[i]<0)) { con[i] = 0; rlx = 1; }
        if ((con[i]==1) && (b[i]>0)) { con[i] = 0; rlx = 1; }
      }

      if (rlx) /* reconstruct the constraint matrix */
      { setzero(L,p*p); m = 0;
        for (i=0; i<p; i++) if (con[i] != 0)
        { L[m*p+i] = 1;
          m++;
        }
      }
    }
  } while (rlx);

  memcpy(b0,x,p*sizeof(double)); /* this is how far we should move from x0 */
  return(0);
}

double maxquadstep(f,x,ymax,h,ll,ur,err,con)
double (*f)(), *x, ymax, *h, *ll, *ur;
int *err, *con;
{ jacobian J;
  double b[2], c[2], d, jwork[12];
  double x0, x1, y0, y1, ym, h0, xl[2], xu[2], xi[2];
  int i, m;
  
  *err = 0;

/* first, can we relax any of the initial constraints?
 * if so, just do one step away from the boundary, and
 * return for restart.
 */
  for (i=0; i<2; i++)
    if (con[i] != 0)
    { xi[0] = x[0]; xi[1] = x[1];
      xi[i] = x[i]-con[i]*h[i];
      y0 = f(xi);
      if (y0>ymax)
      { memcpy(x,xi,2*sizeof(double));
        con[i] = 0;
        return(y0);
      }
    }

/* now, all initial constraints remain active.
 */

  m = 9;
  for (i=0; i<2; i++) if (con[i]==0)
  { m /= 3;
    xl[0] = x[0]; xl[1] = x[1];
    xl[i] = max(x[i]-h[i],ll[i]); y0 = f(xl);
    x0 = xl[i]-x[i]; y0 -= ymax;
    xu[0] = x[0]; xu[1] = x[1];
    xu[i] = min(x[i]+h[i],ur[i]); y1 = f(xu);
    x1 = xu[i]-x[i]; y1 -= ymax;
    if (x0*x1*(x1-x0)==0) { *err = 1; return(0.0); }
    b[i] = (x0*x0*y1-x1*x1*y0)/(x0*x1*(x0-x1));
    c[i] = 2*(x0*y1-x1*y0)/(x0*x1*(x1-x0));
    if (c[i] >= 0.0) { *err = 1; return(0.0); }
    xi[i] = (b[i]<0) ? xl[i] : xu[i];
  }
  else { c[i] = -1.0; b[i] = 0.0; } /* enforce initial constraints */

  if ((con[0]==0) && (con[1]==0))
  { x0 = xi[0]-x[0];
    x1 = xi[1]-x[1];
    ym = f(xi) - ymax - b[0]*x0 - c[0]*x0*x0/2 - b[1]*x1 - c[1]*x1*x1/2;
    d = ym/(x0*x1);
  }
  else d = 0.0;

/* now, maximize the quadratic.
 * y[4] + b0*x0 + b1*x1 + 0.5(c0*x0*x0 + c1*x1*x1 + 2*d*x0*x1)
 * -ve everything, to call quadmin.
 */
  jac_alloc(&J,2,jwork);
  J.Z[0] = -c[0];
  J.Z[1] = J.Z[2] = -d;
  J.Z[3] = -c[1];
  J.st = JAC_RAW;
  J.p = 2;
  b[0] = -b[0]; b[1] = -b[1];
  *err = boxquadmin(&J,b,2,x,ll,ur);
  if (*err) return(ymax);

/* test to see if this step successfully increases...
 */
  for (i=0; i<2; i++)
  { xi[i] = x[i]+b[i];
    if (xi[i]<ll[i]+1e-8*h[i]) xi[i] = ll[i];
    if (xi[i]>ur[i]-1e-8*h[i]) xi[i] = ur[i];
  }
  y1 = f(xi);
  if (y1 < ymax) /* no increase */
  { *err = 1;
    return(ymax);
  }

/* wonderful. update x, h, with the restriction that h can't decrease
 * by a factor over 10, or increase by over 2.
 */
  for (i=0; i<2; i++)
  { x[i] = xi[i];
    if (x[i]==ll[i]) con[i] = -1;
    if (x[i]==ur[i]) con[i] = 1;
    h0 = fabs(b[i]);
    h0 = min(h0,2*h[i]);
    h0 = max(h0,h[i]/10);
    h[i] = min(h0,(ur[i]-ll[i])/2.0); 
  }
  return(y1);
}

double maxbv(f,x,h,ll,ur,m0,m1,tol)
double (*f)(), *x, *h, *ll, *ur, tol;
int m0, m1;
{ double ymax;
  int err, con[2];

  con[0] = con[1] = 0;
  if ((m0>0) & (m1>0))
  {
    ymax = maxbvgrid(f,x,ll,ur,m0,m1,con);
    h[0] = (ur[0]-ll[0])/(2*m0);
    h[1] = (ur[1]-ll[1])/(2*m1);
  }
  else
  { x[0] = (ll[0]+ur[0])/2;
    x[1] = (ll[1]+ur[1])/2;
    h[0] = (ur[0]-ll[0])/2;
    h[1] = (ur[1]-ll[1])/2;
    ymax = f(x);
  }

  while ((h[0]>tol) | (h[1]>tol))
  { ymax = maxbvstep(f,x,ymax,h,ll,ur,&err,con);
    if (err) mut_printf("maxbvstep failure\n");
  }

  return(ymax);
}

double maxbvq(f,x,h,ll,ur,m0,m1,tol)
double (*f)(), *x, *h, *ll, *ur, tol;
int m0, m1;
{ double ymax;
  int err, con[2];

  con[0] = con[1] = 0;
  if ((m0>0) & (m1>0))
  {
    ymax = maxbvgrid(f,x,ll,ur,m0,m1,con);
    h[0] = (ur[0]-ll[0])/(2*m0);
    h[1] = (ur[1]-ll[1])/(2*m1);
  }
  else
  { x[0] = (ll[0]+ur[0])/2;
    x[1] = (ll[1]+ur[1])/2;
    h[0] = (ur[0]-ll[0])/2;
    h[1] = (ur[1]-ll[1])/2;
    ymax = f(x);
  }

  while ((h[0]>tol) | (h[1]>tol))
  { /* first, try a quadratric step */
    ymax = maxquadstep(f,x,ymax,h,ll,ur,&err,con);
    /* if the quadratic step fails, move the grid around */
    if (err)
    {
      ymax = maxbvstep(f,x,ymax,h,ll,ur,&err,con);
      if (err)
      { mut_printf("maxbvstep failure\n");
        return(ymax);
      }
    }
  }

  return(ymax);
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
#include "mut.h"

prf mut_printf = (prf)printf;

void mut_redirect(newprf)
prf newprf;
{ mut_printf = newprf;
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
/*
 * function to find order of observations in an array.
 *
 * mut_order(x,ind,i0,i1)
 *   x     array to find order of.
 *   ind   integer array of indexes.
 *   i0,i1 (integers) range to order.
 *
 * at output, ind[i0...i1] are permuted so that
 *   x[ind[i0]] <= x[ind[i0+1]] <= ... <= x[ind[i1]].
 * (with ties, output order of corresponding indices is arbitrary).
 * The array x is unchanged.
 *
 * Typically, if x has length n, then i0=0, i1=n-1 and
 * ind is (any permutation of) 0...n-1.
 */

#include "mut.h"

double med3(x0,x1,x2)
double x0, x1, x2;
{ if (x0<x1)
  { if (x2<x0) return(x0);
    if (x1<x2) return(x1);
    return(x2);
  }
/* x1 < x0 */
  if (x2<x1) return(x1);
  if (x0<x2) return(x0);
  return(x2);
}

void mut_order(x,ind,i0,i1)
double *x;
int *ind, i0, i1;
{ double piv;
  int i, l, r, z;

  if (i1<=i0) return;
  piv = med3(x[ind[i0]],x[ind[i1]],x[ind[(i0+i1)/2]]);
  l = i0; r = i0-1;

/* at each stage,
 * x[i0..l-1] < piv
 * x[l..r] = piv
 * x[r+1..i-1] > piv
 * then, decide where to put x[i].
 */
  for (i=i0; i<=i1; i++)
  { if (x[ind[i]]==piv)
    { r++;
      z = ind[i]; ind[i] = ind[r]; ind[r] = z;
    } 
    else if (x[ind[i]]<piv)
    { r++;
      z = ind[i]; ind[i] = ind[r]; ind[r] = ind[l]; ind[l] = z;
      l++;
    }
  }

  if (l>i0) mut_order(x,ind,i0,l-1);
  if (r<i1) mut_order(x,ind,r+1,i1);
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
#include "mut.h"

#define LOG_2           0.6931471805599453094172321214581765680755
#define IBETA_LARGE     1.0e30
#define IBETA_SMALL     1.0e-30
#define IGAMMA_LARGE    1.0e30
#define DOUBLE_EP     2.2204460492503131E-16

double ibeta(x, a, b)
double x, a, b;
{ int flipped = 0, i, k, count;
  double I = 0, temp, pn[6], ak, bk, next, prev, factor, val;
  if (x <= 0) return(0);
  if (x >= 1) return(1);
/* use ibeta(x,a,b) = 1-ibeta(1-x,b,z) */
  if ((a+b+1)*x > (a+1))
  { flipped = 1;
    temp = a;
    a = b;
    b = temp;
    x = 1 - x;
  }
  pn[0] = 0.0;
  pn[2] = pn[3] = pn[1] = 1.0;
  count = 1;
  val = x/(1.0-x);
  bk = 1.0;
  next = 1.0;
  do
  { count++;
    k = count/2;
    prev = next;
    if (count%2 == 0)
      ak = -((a+k-1.0)*(b-k)*val)/((a+2.0*k-2.0)*(a+2.0*k-1.0));
    else
      ak = ((a+b+k-1.0)*k*val)/((a+2.0*k)*(a+2.0*k-1.0));
    pn[4] = bk*pn[2] + ak*pn[0];
    pn[5] = bk*pn[3] + ak*pn[1];
    next = pn[4] / pn[5];
    for (i=0; i<=3; i++)
      pn[i] = pn[i+2];
    if (fabs(pn[4]) >= IBETA_LARGE)
      for (i=0; i<=3; i++)
        pn[i] /= IBETA_LARGE;
    if (fabs(pn[4]) <= IBETA_SMALL)
      for (i=0; i<=3; i++)
        pn[i] /= IBETA_SMALL;
  } while (fabs(next-prev) > DOUBLE_EP*prev);
  /* factor = a*log(x) + (b-1)*log(1-x);
  factor -= mut_lgamma(a+1) + mut_lgamma(b) - mut_lgamma(a+b); */
  factor = dbeta(x,a,b,1) + log(x/a);
  I = exp(factor) * next;
  return(flipped ? 1-I : I);
}

/*
 * Incomplete gamma function.
 * int_0^x u^{df-1} e^{-u} du / Gamma(df).
 */
double igamma(x, df)
double x, df;
{ double factor, term, gintegral, pn[6], rn, ak, bk;
  int i, count, k;
  if (x <= 0.0) return(0.0);

  if (df < 1.0)
    return( dgamma(x,df+1.0,1.0,0) + igamma(x,df+1.0) );

  factor = x * dgamma(x,df,1.0,0);
  /* factor = exp(df*log(x) - x - lgamma(df)); */

  if (x > 1.0 && x >= df)
  {
    pn[0] = 0.0;
    pn[2] = pn[1] = 1.0;
    pn[3] = x;
    count = 1;
    rn = 1.0 / x;
    do
    { count++;
      k = count / 2;
      gintegral = rn;
      if (count%2 == 0)
      { bk = 1.0;
        ak = (double)k - df;
      } else
      { bk = x;
        ak = (double)k;
      }
      pn[4] = bk*pn[2] + ak*pn[0];
      pn[5] = bk*pn[3] + ak*pn[1];
      rn = pn[4] / pn[5];
      for (i=0; i<4; i++)
        pn[i] = pn[i+2];
      if (pn[4] > IGAMMA_LARGE)
        for (i=0; i<4; i++)
          pn[i] /= IGAMMA_LARGE;
    } while (fabs(gintegral-rn) > DOUBLE_EP*rn);
    gintegral = 1.0 - factor*rn;
  }
  else
  { /*   For x<df, use the series
     *   dpois(df,x)*( 1 + x/(df+1) + x^2/((df+1)(df+2)) + ... )
     *   This could be slow if df large and x/df is close to 1.
     */
    gintegral = term = 1.0;
    rn = df;
    do
    { rn += 1.0;
      term *= x/rn;
      gintegral += term;
    } while (term > DOUBLE_EP*gintegral);
    gintegral *= factor/df;
  }
  return(gintegral);
}

double pf(q, df1, df2)
double q, df1, df2;
{ return(ibeta(q*df1/(df2+q*df1), df1/2, df2/2));
}
/*
 * Copyright 1996-2006 Catherine Loader.
 */
#include "mut.h"
#include <string.h>

/* quadmin: minimize the quadratic,
 *   2<x,b> + x^T A x.
 *   x = -A^{-1} b.
 *
 * conquadmin: min. subject to L'x = d   (m constraints)
 *   x = -A^{-1}(b+Ly)  (y = Lagrange multiplier)
 *   y = -(L'A^{-1}L)^{-1} (L'A^{-1}b)
 *   x = -A^{-1}b + A^{-1}L (L'A^{-1}L)^{-1} [(L'A^{-1})b + d]
 * (non-zero d to be added!!)
 *
 * qprogmin: min. subject to L'x >= 0.
 */

void quadmin(J,b,p)
jacobian *J;
double *b;
int p;
{ int i;
  jacob_dec(J,JAC_CHOL);
  i = jacob_solve(J,b);
  if (i<p) mut_printf("quadmin singular %2d %2d\n",i,p);
  for (i=0; i<p; i++) b[i] = -b[i];
}

/* project vector a (length n) onto
 * columns of X (n rows, m columns, organized by column).
 * store result in y.
 */
#define pmaxm 10
#define pmaxn 100
void project(a,X,y,n,m)
double *a, *X, *y;
int n, m;
{ double xta[pmaxm], R[pmaxn*pmaxm];
  int i;

  if (n>pmaxn) mut_printf("project: n too large\n");
  if (m>pmaxm) mut_printf("project: m too large\n");

  for (i=0; i<m; i++) xta[i] = innerprod(a,&X[i*n],n);
  memcpy(R,X,m*n*sizeof(double));
  qr(R,n,m,NULL);
  qrsolv(R,xta,n,m);

  matrixmultiply(X,xta,y,n,m,1);
}

void resproj(a,X,y,n,m)
double *a, *X, *y;
int n, m;
{ int i;
  project(a,X,y,n,m);
  for (i=0; i<n; i++) y[i] = a[i]-y[i];
}

/* x = -A^{-1}b + A^{-1}L (L'A^{-1}L)^{-1} [(L'A^{-1})b + d] */
void conquadmin(J,b,n,L,d,m)
jacobian *J;
double *b, *L, *d;
int m, n;
{ double bp[10], L0[100];
  int i, j;

  if (n>10) mut_printf("conquadmin: max. n is 10.\n");
  memcpy(L0,L,n*m*sizeof(double));
  jacob_dec(J,JAC_CHOL);
  for (i=0; i<m; i++) jacob_hsolve(J,&L[i*n]);
  jacob_hsolve(J,b);

  resproj(b,L,bp,n,m);

  jacob_isolve(J,bp);
  for (i=0; i<n; i++) b[i] = -bp[i];

  qr(L,n,m,NULL);
  qrsolv(L,d,n,m);
  for (i=0; i<n; i++)
  { bp[i] = 0;
    for (j=0; j<m; j++) bp[i] += L0[j*n+i]*d[j];
  }
  jacob_solve(J,bp);
  for (i=0; i<n; i++) b[i] += bp[i];
}

void qactivemin(J,b,n,L,d,m,ac)
jacobian *J;
double *b, *L, *d;
int m, n, *ac;
{ int i, nac;
  double M[100], dd[10];
  nac = 0;
  for (i=0; i<m; i++) if (ac[i]>0)
  { memcpy(&M[nac*n],&L[i*n],n*sizeof(double));
    dd[nac] = d[i];
    nac++;
  }
  conquadmin(J,b,n,M,dd,nac);
}

/* return 1 for full step; 0 if new constraint imposed. */
int movefrom(x0,x,n,L,d,m,ac)
double *x0, *x, *L, *d;
int n, m, *ac;
{ int i, imin;
  double c0, c1, lb, lmin;
  lmin = 1.0;
  for (i=0; i<m; i++) if (ac[i]==0)
  { c1 = innerprod(&L[i*n],x,n)-d[i];
    if (c1<0.0)
    { c0 = innerprod(&L[i*n],x0,n)-d[i];
      if (c0<0.0)
      { if (c1<c0) { lmin = 0.0; imin = 1; }
      }
      else
      { lb = c0/(c0-c1);
        if (lb<lmin) { lmin = lb; imin = i; }
      }
    }
  }
  for (i=0; i<n; i++)
    x0[i] = lmin*x[i]+(1-lmin)*x0[i];
  if (lmin==1.0) return(1);
  ac[imin] = 1;
  return(0);
}

int qstep(J,b,x0,n,L,d,m,ac,deac)
jacobian *J;
double *b, *x0, *L, *d;
int m, n, *ac, deac;
{ double x[10];
  int i;

  if (m>10) mut_printf("qstep: too many constraints.\n");
  if (deac)
  { for (i=0; i<m; i++) if (ac[i]==1)
    { ac[i] = 0;
      memcpy(x,b,n*sizeof(double));
      qactivemin(J,x,n,L,d,m,ac);
      if (innerprod(&L[i*n],x,n)>d[i]) /* deactivate this constraint; should rem. */
        i = m+10;
      else
        ac[i] = 1;
    }
    if (i==m) return(0); /* no deactivation possible */
  }

  do
  { if (!deac)
    { memcpy(x,b,n*sizeof(double));
      qactivemin(J,x,n,L,d,m,ac);
    }
    i = movefrom(x0,x,n,L,d,m,ac);

    deac = 0;
  } while (i==0);
  return(1);
}

/*
 * x0 is starting value; should satisfy constraints.
 * L is n*m constraint matrix.
 * ac is active constraint vector:
 *   ac[i]=0, inactive.
 *   ac[i]=1, active, but can be deactivated.
 *   ac[i]=2, active, cannot be deactivated.
 */

void qprogmin(J,b,x0,n,L,d,m,ac)
jacobian *J;
double *b, *x0, *L, *d;
int m, n, *ac;
{ int i;
  for (i=0; i<m; i++) if (ac[i]==0)
  { if (innerprod(&L[i*n],x0,n) < d[i]) ac[i] = 1; }
  jacob_dec(J,JAC_CHOL);
  qstep(J,b,x0,n,L,d,m,ac,0);
  while (qstep(J,b,x0,n,L,d,m,ac,1));
}

void qpm(A,b,x0,n,L,d,m,ac)
double *A, *b, *x0, *L, *d;
int *n, *m, *ac;
{ jacobian J;
  double wk[1000];
  jac_alloc(&J,*n,wk);
  memcpy(J.Z,A,(*n)*(*n)*sizeof(double));
  J.p = *n;
  J.st = JAC_RAW;
  qprogmin(&J,b,x0,*n,L,d,*m,ac);
}
