/*
 * Copyright 1996-2006 Catherine Loader.
 */
/*
 *   Headers for math utility functions.
 */

#ifndef I_MUT_H
#define I_MUT_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef int (*prf)();
extern prf mut_printf;
extern void mut_redirect();

typedef struct {
  double *Z;   /* jacobian matrix, length p*p          */
  double *Q;   /* eigenvalue matrix, length p*p        */
  double *wk;  /* work vector in eig_solve, length p   */
  double *dg;  /* diag vector in eigd, length p        */
  int p;       /* dimension */
  int st;      /* status    */
  int sm;      /* requested decomposition */
} jacobian;

/* m_jacob.c */
extern int jac_reqd();
extern double *jac_alloc();
extern void jacob_dec(),   chol_dec(),   eig_dec();
extern int  jacob_solve(), chol_solve(), eig_solve();
extern int  jacob_hsolve(),chol_hsolve(),eig_hsolve();
extern int  jacob_isolve(),chol_isolve(),eig_isolve();
extern double jacob_qf(),  chol_qf(),    eig_qf();
extern int  jacob_mult(),  chol_mult();

/* m_max.c */
extern double max_grid(), max_golden(), max_quad(), max_nr();

/* maxbv.c */
extern double maxbvgrid(), maxbvstep(), maxquadstep(), maxbv(), maxbvq();

/* quadmin.c */
extern void quadmin(), project(), resproj(), conquadmin();

/* m_qr.c */
extern void qr(), qrinvx(), qrtinvx(), qrsolv();

/* m_svd.c */
extern void svd(), hsvdsolve();
extern int svdsolve();

/* m_solve.c */
extern double solve_bisect(), solve_secant(), solve_nr(), solve_fp();

/* m_vector.c */
extern void setzero(), unitvec(), addouter();
extern void matrixmultiply(), multmatscal(), transpose();
extern double innerprod(), m_trace(), vecsum();

#define BDF_NONE  0
#define BDF_EXPLEFT  1
#define BDF_EXPRIGHT 2

/* return codes for functions optimized by max_nr */
#define NR_OK 0
#define NR_INVALID 1
#define NR_BREAK   2
#define NR_REDUCE  3
#define NR_NCON  10
#define NR_NDIV  11


/* jacobian status definitions */
#define JAC_RAW 0
#define JAC_CHOL 1
#define JAC_EIG  2
#define JAC_EIGD 3

/*  Numerical Integration Stuff
 */
#define MXRESULT 5
#define MXIDIM  10  /* max. dimension */
extern void simpsonm(), simpson4(), integ_disc(), integ_circ();
extern void integ_sphere(), monte(), rn3();
extern double simpson(), sptarea();

/*  Density, distribution stuff
 */

#ifndef PI
#define PI          3.141592653589793238462643
#endif
#define S2PI        2.506628274631000502415765       /* sqrt(2*pi) */
#define PIx2        6.283185307179586476925286       /* 2*pi */
#define SQRPI       1.772453850905516027298          /* sqrt(pi)   */
#define HF_LG_PIx2  0.918938533204672741780329736406 /* 0.5*log(2*pi) */
#define SQRT2       1.4142135623730950488            /* sqrt(2.0) */
#define gold_rat    0.6180339887498948482045870      /* (sqrt(5.0)-1)/2 */

#define LOG_ZERO -1e100
#define D_0 ((give_log) ? LOG_ZERO : 0.0)
#define D_1 ((give_log) ? 0.0 : 1.0)
#define DEXP(x)   ((give_log) ? (x) : exp(x))
#define FEXP(f,x) ((give_log) ? -0.5*log(f)+(x) : exp(x)/sqrt(f))

#define INVALID_PARAMS 0.0

extern double stirlerr(), bd0();
extern double dbinom_raw(), dpois_raw();
extern double dbinom(), dpois(), dnbinom(), dbeta(), dgamma(), dt(), df(), dhyper();
extern double dchisq();

extern double igamma(), ibeta();
extern double pf(), pchisq();
#define pchisq(x,df) igamma((x)/2.0,(df)/2.0)

/* order.c */
extern void mut_order();

/* erf.c */
extern double mut_erf(), mut_erfc(), mut_pnorm();

/* lgamma.c */
extern double mut_lgamma(), mut_lgammai();

/* math.c */
extern double mut_exp(), mut_daws(), ptail(), logit(), expit();
extern int factorial();
#endif  /* define I_MUT_H */
