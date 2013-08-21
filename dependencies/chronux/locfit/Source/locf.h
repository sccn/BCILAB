/*
 * Copyright 1996-2006 Catherine Loader.
 */
#ifndef I_LOCF_H
#define I_LOCF_H

#include "stdlib.h"
#include "ctype.h"
#include "mut.h"

#ifdef WARN
#undef WARN
#endif

#define LERR(args) {mut_printf("Error: "); mut_printf args; mut_printf("\n"); lf_error=1;}
#define WARN(args)  {mut_printf("Warning: "); mut_printf args; mut_printf("\n"); }

extern int lf_error;
#define LOGPI 1.144729885849400174143427
#define HUBERC 2.0
#define NOSLN 0.1278433
#define GFACT 2.5
#define EFACT 3.0
#define ISWAP(a,b) { int zz; zz = a; a = b; b = zz; }
#define MAX(a,b) (((a)>(b)) ? (a) : (b))
#define MIN(a,b) (((a)<(b)) ? (a) : (b))
#define SGN(x) (((x)>0) ? 1 : -1)
#define SQR(x) ((x)*(x))

extern int lf_error;

typedef struct {
  double *wk, *coef, *xbar, *f;
  jacobian xtwx;
  int lwk, haspc;
} paramcomp;
#define haspc(pc) ((pc)->haspc)

/*
 * MXDIM and MXDEG are maximum dimension and local polynomial
 * degree for Locfit. Note that some parts of the code may be
 * more restrictive.
 */

#define MXDIM 15
#define MXDEG 7

typedef struct {
  double *x[MXDIM];
  double *y;
  double *w;
  double *b;
  double *c;
  double sca[MXDIM];
  double xl[2*MXDIM];
  int n, d, ord;
  int sty[MXDIM];
} lfdata;
#define resp(lfd,i) (((lfd)->y==NULL) ? 0.0 : (lfd)->y[i])
#define base(lfd,i) (((lfd)->b==NULL) ? 0.0 : (lfd)->b[i])
#define prwt(lfd,i) (((lfd)->w==NULL) ? 1.0 : (lfd)->w[i])
#define cens(lfd,i) (((lfd)->c==NULL) ? 0 : (int)(lfd)->c[i])
#define datum(lfd,i,j) ((lfd)->x[i][j])
#define dvari(lfd,i)   ((lfd)->x[i])


/*
 *   The design structure used in Locfit, and associated macro definitions.
 */

typedef struct {
  int des_init_id;
  double *wk;
  int *ind;
  int lwk, lind;

  double *xev;                /* fitting point, length p          */
  double *X;                  /* design matrix, length n*p        */
  double *w, *di, *res, *th, *wd, h;
  double *V, *P;              /* matrices with length p*p         */
  double *f1, *ss, *oc, *cf;  /* work vectors, length p  */
  double llk, smwt;
  double tr0, tr1, tr2;       /* traces for local df computation */
  jacobian xtwx;     /* to store X'WVX and decomposition */
  int cfn[1+MXDIM], ncoef;
  int *fix;          /* integer vector for fixed coefficients. */
  int (*itype)();    /* density integration function     */
  int n, p;
} design;

#define cfn(des,i) (des->cfn[i])
#define d_x(des) ((des)->X)
#define d_xi(des,i) (&(des)->X[i*((des)->p)])
#define d_xij(des,i,j) ((des)->X[i*((des)->p)+j])
#define is_fixed(des,i) ((des)->fix[i]==1)
#define wght(des,i) ((des)->w[i])
#define dist(des,i) ((des)->di[i])
#define fitv(des,i) ((des)->th[i])
#define DES_INIT_ID 34988372

extern int des_reqd(), des_reqi();

typedef struct {
  int deflink, canlink, quasi, robust;
  int (*vallink)(), (*family)(), (*initial)(), (*like)(), (*pcheck)();
} family;
#define isquasi(fam) ((fam)->quasi)
#define isrobust(fam) ((fam)->robust)
extern int inllmix; /* flag needed to ensure correct behavior in llmix. */

typedef struct {
  double nn, fixh, adpen;
  int ker, kt;
  int deg, deg0, p;
  int acri;
  int fam, lin;
  family fami;
  int ubas;
  double (*vb)();
  void (*vbasis)();
} smpar;
#define nn(sp)   ((sp)->nn)
#define fixh(sp) ((sp)->fixh)
#define pen(sp)  ((sp)->adpen)
#define ker(sp)  ((sp)->ker)
#define kt(sp)   ((sp)->kt)
#define deg(sp)  ((sp)->deg)
#define deg0(sp) ((sp)->deg0)
#define npar(sp) ((sp)->p)
#define acri(sp) ((sp)->acri)
#define ubas(sp) ((sp)->ubas)
#define fam(sp)  ((sp)->fam)
#define fami(sp)  (&(sp)->fami)
#define link(sp) ((sp)->lin)

typedef struct {
  int deriv[MXDEG+2];
  int nd;
} deriv;

/*
 *  Criteria for adaptive local fitting  mi[MACRI]
 *  1: localized CP;  2: ICI (katkovnik);  3: curvature model index
 *  4: Increase bandwidth until locfit returns LF_OK
 */
#define ANONE 0
#define ACP  1
#define AKAT 2
#define AMDI 3
#define AOK  4

/*
 * weight functions mi[MKER].
 * see Table 3.1 or the function W() in weights.c for definitions.
 */
#define WRECT 1
#define WEPAN 2
#define WBISQ 3
#define WTCUB 4
#define WTRWT 5
#define WGAUS 6
#define WTRIA 7
#define WQUQU 8
#define W6CUB 9
#define WMINM 10
#define WEXPL 11
#define WMACL 12
#define WPARM 13

/*
 * type of multivariate weight function mi[MKT]
 * KSPH (spherical)  KPROD (product)
 * others shouldn't be used at present.
 */
#define KSPH   1
#define KPROD  2
#define KCE    3
#define KLM    4
#define KZEON  5

/*
 * Local likelihood family mi[MTG]
 * for quasi-likelihood, add 64.
 */
#define TNUL 0
#define TDEN 1
#define TRAT 2
#define THAZ 3
#define TGAUS 4
#define TLOGT 5
#define TPOIS 6
#define TGAMM 7
#define TGEOM 8
#define TCIRC 9
#define TROBT 10
#define TRBIN 11
#define TWEIB 12
#define TCAUC 13
#define TPROB 14
#define TQUANT 15

/*
 * Link functions mi[MLINK].
 * Mostly as in table 4.1 of the book.
 * LDEFAU and LCANON are used to select default and canonical
 * links respectively. LINIT shouldn't be selected by user...
 */
#define LINIT  0
#define LDEFAU 1
#define LCANON 2
#define LIDENT 3
#define LLOG   4
#define LLOGIT 5
#define LINVER 6
#define LSQRT  7
#define LASIN  8

/*
 * components of vector returned by the links() function
 * in family.c. ZLIK the likelihood; ZMEAN = estimated mean;
 * ZDLL = derivative of log-likelihood; ZDDLL = - second derivative
 */
#define LLEN  4
#define ZLIK  0
#define ZMEAN 1
#define ZDLL  2
#define ZDDLL 3

/*
 * return status for the locfit() function
 */
#define LF_OK   0
#define LF_DONE 1   /* done - forced break from iterations */
#define LF_OOB  2   /* out of bounds, or large unstable parameter */
#define LF_PF   3   /* perfect fit; interpolation; deviance=0 */
#define LF_NCON 4   /* not converged */
#define LF_NSLN 5   /* no solution - eg separation in binomial. */
#define LF_NOPT 6   /* no or insufficient points with non-zero wt */
#define LF_INFA 7   /* initial failure e.g. log(0) */
#define LF_DEMP 10  /* density -- empty integration region */
#define LF_XOOR 11  /* density -- fit point outside xlim region */
#define LF_DNOP 12  /* density version of 6 */
#define LF_BADP 81  /* bad parameters e.g. neg prob for binomial */
#define LF_LNK  82  /* invalid link */
#define LF_FAM  83  /* invalid family */
#define LF_ERR  99  /* error */

#define STANGL 4
#define STLEFT 5
#define STRIGH 6
#define STCPAR 7

/*
 * Integration type mi[MIT] for integration in
 * density estimation.
 */
#define INVLD 0
#define IDEFA 1
#define IMULT 2
#define IPROD 3
#define IMLIN 4
#define IHAZD 5
#define ISPHR 6
#define IMONT 7

/* density.c */
extern int densinit(), likeden(), deitype();
extern int fact[];
extern void prodintresp(), prresp();
extern int de_mint, de_itype, de_renorm;

/* dens_haz.c */
extern void haz_init();
extern int hazint();

/* dens_odi.c */
extern int onedint();
extern void recurint();

/* famquant.c */
extern void lfquantile();

/* family.c */
extern int lffamily(), lflink();
extern int links(), stdlinks(), defaultlink(), validlinks();
extern double b2(), b3(), b4(), lf_link(), invlink();
extern void setfamily();

/* lf_adap.c */
extern int alocfit(), lfacri();

/* lf_fitfun.c */
extern void fitfun(), makecfn(), designmatrix();
extern int calcp(), coefnumber();

/* lf_nbhd.c */
extern double kordstat(), rho();
extern void nbhd();

/* lf_robust.c */
extern double median();
extern void lf_robust();

/* lfstr.c */
extern int pmatch();

/* lf_vari.c */
extern void lf_vcov(), comp_vari(), local_df();

/* lf_wdiag.c */
extern int wdiag(), wdiagp();

/* locfit.c */
extern int locfit(), des_reqd(), des_reqi(), likereg();
extern int reginit();
extern void lfdata_init(), smpar_init(), deriv_init(), des_init(), lfiter();
extern int lf_maxit, lf_debug;
extern void lf_status_msg();

/* minmax.c */
extern double ipower(), minmax();

/* weight.c */
extern int lfkernel(), lfketype();
extern double W(), weight(), weightd(), Wd(), Wdd(), wint();
extern double Wconv(), Wconv1(), Wconv4(), Wconv5(), Wconv6(), Wikk();
extern int iscompact(), wtaylor();

#endif /* define I_LOCF_H */
