#include "local.h"

int lf_error;

  lfit lf;
  design des;

void calllf(x,y,n)
double *x, *y;
int n;
{ lfdata *lfd;
  fitpt *fp;
  int i;

  lfit_init(&lf);

  lfd = &(lf.lfd);
  lfd->x[0] = x;
  lfd->y = y;
  lfd->n = n;
  lfd->d = 1;

  startlf(&des,&lf,procv,1);
//  startmodule(&lf,&des,"std",0);

  fp = &lf.fp;
  for (i=0; i<fp->nv; i++)
    printf("%8.5f %8.5f\n",evptx(fp,i,0),fp->coef[i]);
}

int main()
{ double x[10], y[10];

  x[0] = 0; x[1] = 1; x[2] = 2; x[3] = 3; x[4] = 4;
  x[5] = 5; x[6] = 6; x[7] = 7; x[8] = 8; x[9] = 9;

  y[0] = 0.3692449; y[1] = 0.8194270;
  y[2] = 1.6363139; y[3] =-0.9969944;
  y[4] = 0.5359200; y[5] = 1.8642622;
  y[6] = 0.3568127; y[7] = 0.4746753;
  y[8] =-2.0038246; y[9] = 1.6636109;

  calllf(x,y,10);
}
