/*
 *   The design structure used in Locfit, and associated macro definitions.
 */

typedef struct {
  int des_init_id;
  double *wk;
  Sint *ind;
  int lwk, lind;

  double *xev;                /* fitting point, length p          */
  double *X;                  /* design matrix, length n*p        */
  double *w, *di, *res, *th, *wd, h;
  double *V, *P;              /* matrices with length p*p         */
  double *f1, *ss, *oc, *cf;  /* work vectors, length p  */
  double llk, smwt;
  jacobian xtwx;     /* to store X'WVX and decomposition */
  int cfn[1+MXDIM], ncoef;
  Sint *fix;         /* integer vector for fixed coefficients. */
  int (*itype)();    /* density integration function     */
  int n, p;
  int (*vfun)();     /* pointer to the vertex processing function. */
} design;

#define cfn(des,i) (des->cfn[i])
#define d_x(des) ((des)->X)
#define d_xi(des,i) (&(des)->X[i*((des)->p)])
#define d_xij(des,i,j) ((des)->X[i*((des)->p)+j])
#define is_fixed(des,i) ((des)->fix[i]==1)
#define DES_INIT_ID 34988372

extern int des_reqd(), des_reqi();
