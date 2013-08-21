extern int procv_llmix();

/* adap.c */
extern int alocfit();

/* band.c */
extern void band(), kdeselect(), kdecri();

/* density.c */
extern int densinit(), likeden();
extern int fact[];
extern void prodintresp(), prresp();
extern int de_mint, de_itype, de_renorm;

/* dens_haz.c */
extern void haz_init();
extern int hazint();

/* dens_int.c */
extern double dens_integrate();
extern void dens_renorm(), dens_lscv(), lforder();

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
extern void trchck(), guessnv(), lfit_alloc();
extern void dataf(), gridf(), crossf(), xbarf(), preset();
extern int findpt(), newsplit(), lfit_reqd(), lfit_reqi();

/* ev_trian.c */
extern void triang_start(), triang_grow(), triang_guessnv();
extern double triang_int();

/* family.c */
extern int links(), stdlinks(), defaultlink(), validlinks();
extern double b2(), b3(), b4(), lf_link(), invlink();
extern void setfamily();

/* fitted.c */
extern void fitted();

/* frend.c */
extern void ressumm();
extern double rss();

/* lf_dercor.c */
extern void dercor();

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
extern int lffamily(), lfkernel(), lfketype(), lflink();
extern int deitye(), lfevstr(), lfacri();
extern int ppwhat(), restyp();

/* lf_vari.c */
extern void lf_vcov(), comp_vari(), local_df();

/* locfit.c */
extern int locfit(), des_reqd(), des_reqi(), likereg();
extern int reginit();
extern void lfdata_init(), smpar_init(), deriv_init(), des_init(), lfiter();
extern int lf_maxit, lf_debug;

/* minmax.c */
extern double ipower(), minmax();

/* odint.c */
extern int onedint();
extern void recurint();

/* pcomp.c */
extern double addparcomp();
extern void compparcomp(), subparcomp(), subparcomp2(), pcchk();
extern int pc_reqd(), noparcomp();

/* preplot.c */
extern void preplot(), cpreplot();
extern int setpppoints();

/* procv.c */
extern int procvhatm(), procv(), procvraw(), procveraw(), procvvord(), calcp();

/* resid.c */
extern double resid();

/* scb.c */
extern void scb(), cscbsim();

/* scb_iface.c */
extern int constants();

/* simul.c */
extern void liksim(), scbsim(), scbmax(), regband(), rband();

/* startlf.c */
extern void set_flim(), set_scales(), nstartlf(), startlf(), lfit_init();
extern void fitoptions(), clocfit(), endfit(), startmodule();
extern int nofit(), initmodule();

/* strings.c */
extern int stm(), pmatch(), matchlf(), matchrt(), checkltor(), checkrtol();
extern void strip();

/* wdiag.c */
extern int wdiag(), wdiagp();

/* weight.c */
extern double W(), weight(), weightd(), Wd(), Wdd(), wint();
extern double Wconv(), Wconv1(), Wconv4(), Wconv5(), Wconv6(), Wikk();
extern int iscompact(), wtaylor();

extern void initsimple(), initstd(), inithatm(), initgam();
extern void initlscv(), initrband(), initscb(), initkappa();
