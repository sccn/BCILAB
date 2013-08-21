/*
 * Copyright 1996-2006 Catherine Loader.
 */
#ifndef I_LFEV_H
#define I_LFEV_H

#include "locf.h"
#include "tube.h"

/* uncomment the following line to allow plug-in modules. */
/*#define ALLOW_MODULES*/

#define DALP 0
#define DFXH 1

typedef struct {
  int ev;
  double *sv;
  double cut;
  double fl[2*MXDIM];
  int *iwk, *ce, *s, *lo, *hi;
  int liw, nce, ncm, maxk;
  int mg[MXDIM];
  void (*espec)();
} evstruc;
#define ev(evs)   ((evs)->ev)
#define cut(evs)  ((evs)->cut)
#define mk(evs)   ((evs)->maxk)
#define mg(evs)   ((evs)->mg)

typedef struct {
  double *xev, *wk, *coef, *nlx, *t0, *lik, *h, *deg, *kap;
  int lev, lwk;
  int d, dcor, geth, hasd;
  int nv, nvm;
} fitpt;
#define evp(fp)     ((fp)->xev)
#define evpt(fp,i)  (&(fp)->xev[(i)*(fp)->d])
#define evptx(fp,i,k) ((fp)->xev[(i)*(fp)->d+(k)])
#define llk(fp) ((fp)->kap[0])
#define df0(fp) ((fp)->kap[1])
#define df1(fp) ((fp)->kap[2])
#define rv(fp)  ((fp)->kap[3])
#define rsc(fp) ((fp)->kap[5])
#define dc(fp)  ((fp)->dcor)
#define geth(fp)  ((fp)->geth)

typedef struct {
  int (*procv)(), keepv, keepc, nopc, isset, np;
  void (*alloc)();
  void (*pproc)();
  double *params;
} module;
#define PROC_VERTEX(des,lf,i) (lf)->mdl.procv((des),(lf),(i))
#define MODPARAMS(lf) ((lf)->mdl.params)
#define MODPARAM(lf,i) ((lf)->mdl.params[i])
#define MODNPARAMS(lf) ((lf)->mdl.np)
#define PROCV(lf) (lf)->mdl.procv
#define ALLOC(lf) (lf)->mdl.alloc
#define PPROC(lf) (lf)->mdl.pproc
#define KEEPC(lf) (lf)->mdl.keepc
#define KEEPV(lf) (lf)->mdl.keepv
#define NOPC(lf)  (lf)->mdl.nopc

typedef struct {
  int       lf_init_id;
  lfdata    lfd;
  smpar     sp;
  evstruc   evs;
  fitpt     fp;
  deriv     dv;
  paramcomp pc;
  module    mdl;
  } lfit;
#define LF_INIT_ID 34897239
#define NOBS(lf) ((lf)->lfd.n)
#define NVAR(lf) ((lf)->lfd.d)
#define NPAR(lf) ((lf)->sp.p)

/*
 * VVEC(lf,i) is storage vector for the i'th item.
 * VVAL(lf,v,i) is the storage point for the i'th item when fitting at vertex v.
 *   should have 0 <= i < keepv; keepv defined when initializing module.
 */
#define VVEC(lf,i) (&(lf)->fp.wk[(i)*(lf)->fp.nvm])
#define VVAL(lf,v,i) ((lf)->fp.wk[(i)*(lf)->fp.nvm+(v)])
#define VECR(lf) ((lf)->kap)

/*
 *  mi[MGETH] codes
 *  scb(), pointwise codes are 71,...,75.
 *         add 10 for simultaneous codes.
 */
#define GSTD 0
#define GHAT 1
#define GKAP 2
#define GRBD 3
#define GAMF 4
#define GAMP 5
#define GLSC 6
#define GSMP 7
#define GMIX 8
#define GLM1 71
#define GLM2 72
#define GLM3 73
#define GLM4 74
#define GLDN 75

/* bandwidth criteria */
#define BGCV 1
#define BCP  2
#define BIND 3

/*
 * Evaluation structures
 * EFITP special for `interpolation' at fit points
 */
#define ENULL  0
#define ETREE  1
#define EPHULL 2
#define EDATA  3
#define EGRID  4
#define EKDTR  5
#define EKDCE  6
#define ECROS  7
#define EPRES  8
#define EXBAR  9
#define ENONE  10
#define ESPHR  11
#define EFITP  50
#define ESPEC  100

/*
 * For prediction functions, what to predict?
 * PCOEF -- coefficients        PT0   -- influence function
 * PNLX  -- ||l(x)||            PBAND -- bandwidth h(x)
 * PDEGR -- local poly. degree  PLIK  -- max. local likelihood
 * PRDF  -- local res. d.f.     PVARI -- ||l(x)||^2
 */
#define PCOEF 1
#define PT0   2
#define PNLX  3
#define PBAND 4
#define PDEGR 5
#define PLIK  6
#define PRDF  7
#define PVARI 8

/*
 *  Residual Types
 */
#define RDEV  1
#define RPEAR 2
#define RRAW  3
#define RLDOT 4
#define RDEV2 5
#define RLDDT 6
#define RFIT  7
#define RMEAN 8

/* band.c */
extern void band(), kdeselect(), kdecri(), bselect();

/* dens_int.c */
extern double dens_integrate();
extern void dens_renorm(), lforder();

/* ev_atree.c */
extern void atree_start(), atree_grow(), atree_guessnv();
extern double atree_int();

/* ev_interp.c */
extern double dointpoint(), cubintd();
extern double linear_interp(), cubic_interp(), rectcell_interp();
extern int exvval();
extern void exvvalpv(), hermite2();

/* ev_kdtre.c */
extern void kdtre_start(), kdtre_guessnv();
extern double kdtre_int();

/* ev_sphere.c */
extern void sphere_start(), sphere_guessnv();
extern double sphere_int();

/* ev_main.c */
extern void trchck(), guessnv(), lfit_alloc(), evstruc_alloc(), evstruc_init();
extern void dataf(), gridf(), crossf(), xbarf(), preset();
extern int findpt(), newsplit(), lfit_reqd(), evstruc_reqi();
extern int lfevstr();

/* ev_trian.c */
extern void triang_start(), triang_grow(), triang_guessnv();
extern double triang_int();

/* fitted.c */
extern void fitted();
extern double resid();
extern int restyp();

/* frend.c */
extern void ressumm();
extern double rss();

/* lf_dercor.c */
extern void dercor();

/* pcomp.c */
extern double addparcomp();
extern void compparcomp(), subparcomp(), subparcomp2(), pcchk();
extern int pc_reqd(), noparcomp();

/* preplot.c */
extern void preplot(), cpreplot();
extern int setpppoints(), ppwhat();

/* procv.c */
extern int procv_nov(), procv_var();

/* startlf.c */
extern void set_flim(), set_scales(), nstartlf(), startlf(), lfit_init();
extern void fitoptions(), clocfit(), endfit(), startmodule();
extern int nofit(), initmodule();

extern void initsimple(), initstd(), inithatm(), initgam(), initallcf();
extern void initlscv(), initrband(), initscb(), initkappa(), initvord();

/* modkappa.c */
extern int constants();

/* modrband.c */
extern void rband();

/* modstd.c */
extern int procvstd();
#endif /* define I_LFEV_H */
