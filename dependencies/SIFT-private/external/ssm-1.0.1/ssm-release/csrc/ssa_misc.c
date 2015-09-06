/**************************************\
    Filename: ssa_misc.c
    Author: Jyh-Ying Peng (´^´¼·­)
    Year: 2007
    Version: 1.2
\**************************************/

#include "ssa.h"

#define MAX_RAND_HALF 2147483648.0

/*====================*\
    Global variables
\*====================*/
extern true;
extern false;
extern g_zero;
extern g_one;
extern g_minusone;

/*===========================*\
    Missing data truncation
\*===========================*/
void truncate(int p, int m, int n, int N, const double *y, const int *miss,
    const int *pknown, int pd, const double *H, int Hd, const double *Z, int Zd,
    int *pmiss, int *imiss, double *ymiss, double *Hmiss, double *Zmiss)
{
    /*** Temporary variables ***/
    int i, j, t, pcur;
    
    for(t=0; t<n; t++) {
        pcur        = pd ? pknown[t] : p;
        pmiss[t]    = pcur;
        j           = 0;
        for(i=0; i<pcur; i++) if(miss[i]) pmiss[t]--; else imiss[j++] = i;
        
        if(pmiss[t] < pcur && pmiss[t] > 0) {
            for(i=0; i<pmiss[t]; i++) {
                if(y) for(j=0; j<N; j++) ymiss[j*pcur+i] = y[j*pcur+imiss[i]];
                if(H) for(j=0; j<pmiss[t]; j++) Hmiss[j*pcur+i] = H[imiss[j]*pcur+imiss[i]];
                if(Z) for(j=0; j<m; j++) Zmiss[j*pcur+i] = Z[j*pcur+imiss[i]];
            }
        }
        else {
            if(y) for(i=0; i<pcur*N; i++) ymiss[i] = y[i];
            if(H) for(i=0; i<pcur*pcur; i++) Hmiss[i] = H[i];
            if(Z) for(i=0; i<pcur*m; i++) Zmiss[i] = Z[i];
        }
        
        miss    += p;
        imiss   += p;
        if(y) { y += p*N; ymiss += p*N; }
        if(H) { if(Hd) H += p*p; Hmiss += p*p; }
        if(Z) { if(Zd) Z += p*m; Zmiss += p*m; }
    }
}

/*================*\
    Untruncation
\*================*/
void untruncate(int p, int n, int N, const int *pknown, int pd,
    const int *pmiss, const int *imiss, const double *vmiss, const double *Fmiss,
    double *v, double *F)
{
    int i, j, t, pcur;
    
    for(t=0; t<n; t++) {
        pcur    = pd ? pknown[t] : p;
        
        if(pmiss[t] > 0 && pmiss[t] < pcur) {
            for(i=0; i<pcur*N; i++) v[i] = 0;
            for(i=0; i<pcur*pcur; i++) F[i] = 0;
            for(i=0; i<pmiss[t]; i++) {
                for(j=0; j<N; j++) v[imiss[i]+j*pcur] = vmiss[i+j*pcur];
                for(j=0; j<pmiss[t]; j++) F[imiss[i]+imiss[j]*pcur] = Fmiss[i+j*pcur];
            }
        }
        else {
            for(i=0; i<pcur*N; i++) v[i] = vmiss[i];
            for(i=0; i<pcur*pcur; i++) F[i] = Fmiss[i];
        }
        
        vmiss   += p*N;
        Fmiss   += p*p;
        v       += p*N;
        F       += p*p;
    }
}

/*==================================*\
    Permute 2nd and 3rd dimensions
\*==================================*/
void permute(int p, int n, int N, const double *x, double *y)
{
    int i, j, t;
    
    /*** Reshape p*n*N to p*N*n ***/
    for(t=0; t<n; t++) {
        for(i=0; i<p; i++) for(j=0; j<N; j++) y[i+j*p] = x[i+(t+j*n)*p];
        y   += p*N;
    }
}

/*======================================================*\
    Generate random variables from normal distribution
\*======================================================*/
void genrandn(int N, double *x)
{
    int i;
    
    for(i=0; i<N/2; i++, x+=2) {
        double a, b, R;
        do {
            a   = (genrand_int32() - MAX_RAND_HALF + 0.5)/MAX_RAND_HALF;
            b   = (genrand_int32() - MAX_RAND_HALF + 0.5)/MAX_RAND_HALF;
            R   = a*a + b*b;
        } while(R == 0 || R > 1);
        R       = sqrt(-2*log(R)/R);
        *x      = a*R;
        *(x+1)  = b*R;
    }

    if(N%2) {
        double a, b, R;
        do {
            a   = (genrand_int32() - MAX_RAND_HALF + 0.5)/MAX_RAND_HALF;
            b   = (genrand_int32() - MAX_RAND_HALF + 0.5)/MAX_RAND_HALF;
            R   = a*a + b*b;
        } while(R == 0 || R > 1);
        R       = sqrt(-2*log(R)/R);
        *x      = a*R;
    }
}

/*========================*\
    Transform covariance
\*========================*/
void sigma(int p, int N, const double *C, const double *u, double *x)
{
    int i, j;
    int diag    = true;
    
    /*** Check for diagonal matrix ***/
    for(i=0; i<p; i++)
        for(j=0; j<i; j++) /* Only check the upper triangular part */
            if(C[i*p+j] != 0) { diag = false; break; }
    
    /*** Transform ***/
    if(diag) {
        /*** x += sqrt(diag(C))u ***/
        for(i=0; i<p; i++) {
            double c    = sqrt(C[i*p+i]);
            for(j=0; j<N; j++) x[i+j*p] += c*u[i+j*p];
        }
    }
    else {
        int info;
        double *U       = (double *)calloc(p*p, sizeof(double));
        double *lambda  = (double *)calloc(p, sizeof(double));
        int sym         = true;
        int lwork       = -1;
        double dlwork, *work;
        
        /*** Check for symmetric matrix ***/
        for(i=0; i<p; i++)
            for(j=0; j<i; j++)
                if(C[i*p+j] != C[j*p+i]) {
                    sym = false;
                    break;
                }
        
        if(sym) {
            for(i=0; i<p*p; i++) U[i] = C[i];
            
            dsyev("V", "U", &p, U, &p, lambda, &dlwork, &lwork, &info);
            lwork   = (int)dlwork;
            work    = (double *)calloc(lwork, sizeof(double));
            dsyev("V", "U", &p, U, &p, lambda, work, &lwork, &info);
            free(work);
        }
        else {
            /* Assume there's no complex eigenvalues */
            double *C2      = (double *)calloc(p*p, sizeof(double));
            double *lambda2 = (double *)calloc(p, sizeof(double));

            for(i=0; i<p*p; i++) C2[i] = C[i];

            dgeev("N", "V", &p, C2, &p, lambda, lambda2, (double *)0, &p, U, &p, &dlwork, &lwork, &info);
            lwork   = (int)dlwork;
            work    = (double *)calloc(lwork, sizeof(double));
            dgeev("N", "V", &p, C2, &p, lambda, lambda2, (double *)0, &p, U, &p, work, &lwork, &info);
            free(work);
            
            free(C2);
            free(lambda2);
        }
        
        /*** x += U(sqrt(lambda))u ***/
        for(i=0; i<p; i++) {
            double c    = sqrt(lambda[i]);
            for(j=0; j<p; j++) U[j+i*p] *= c;
        }
        dgemm("N", "N", &p, &N, &p, &g_one, U, &p, u, &p, &g_one, x, &p);
        
        free(U);
        free(lambda);
    }
}

/*============================*\
    Sample state space model
\*============================*/
void sample(int p, int m, int r, int n, int N, const double *H, int Hd, const double *Z, int Zd,
    const double *T, int Td, const double *R, int Rd, const double *Q, int Qd, const double *c, int cd,
    const double *a1, const double *P1, const double *samples, double *y, double *alpha, double *eps, double *eta)
{
    int i, t;
    double *alphaprev;
    const double *alpha10   = samples;
    const double *eps0      = alpha10 + m*N;
    const double *eta0      = eps0 + p*N*n;
    double *y_p0 = y, *alpha_p0 = alpha;

    /*** Initialize state sequence ***/
    for(i=0; i<m*N; i++) alpha[i] = a1[i%m];
    sigma(m, N, P1, alpha10, alpha);
    alphaprev = alpha; alpha += m*N;
    if(cd) {
        for(t=0; t<n-1; t++) {
            for(i=0; i<m*N; i++) alpha[i] = c[i%m];
            alpha   += m*N;
            c       += m;
        }
        alpha       = alpha_p0;
    }
    else for(i=0; i<m*N*(n-1); i++) alpha[i] = c[i%m];
    /* Here alpha(2:n) = c(1:n-1) and c is no longer needed */

    /*** Transform samples to disturbance ***/
    for(i=0; i<p*N*n; i++) eps[i] = 0;
    if(Hd) {
        sigma(p, N, H, eps0, eps);
        for(i=0; i<p*N; i++) y[i] = eps[i];
    }
    else {
        sigma(p, n*N, H, eps0, eps);
        for(i=0; i<p*N*n; i++) y[i] = eps[i];
        /* Here H, eps and eps0 are no longer needed */
    }
    for(i=0; i<r*N*n; i++) eta[i] = 0;
    if(Qd) {
        sigma(r, N, Q, eta0, eta);
        dgemm("N", "N", &m, &N, &r, &g_one, R, &m, eta, &r, &g_one, alpha, &m);
    }
    else {
        int N2  = (n-1)*N;
        sigma(r, n*N, Q, eta0, eta);
        dgemm("N", "N", &m, &N2, &r, &g_one, R, &m, eta, &r, &g_one, alpha, &m);
        /* Here Q, eta and eta0 are no longer needed */
    }
    
    /*** Initialize observation sequence ***/
    if(Zd) dgemm("N", "N", &p, &N, &m, &g_one, Z, &p, alphaprev, &m, &g_one, y, &p);
    
    /*** State space recursion ***/
    for(t=1; t<n; t++) {
        /*** alpha(t) = c(t-1) + T(t-1)alpha(t-1) + R(t-1)eta(t-1) ***/
        dgemm("N", "N", &m, &N, &m, &g_one, T, &m, alphaprev, &m, &g_one, alpha, &m);
        
        alphaprev   = alpha;
        alpha      += m*N;
        if(Td) T += m*m;
        if(Qd) {
            if(Rd) R += m*r;
            Q       += r*r;
            eta0    += r*N;
            eta     += r*N;
            sigma(r, N, Q, eta0, eta);
            /*** alpha(t+1) = c(t) + R(t)eta(t) ***/
            if(t < n-1) dgemm("N", "N", &m, &N, &r, &g_one, R, &m, eta, &r, &g_one, alpha, &m);
        }
        
        y   += p*N;
        if(Hd) {
            H       += p*p;
            eps0    += p*N;
            eps     += p*N;
            sigma(p, N, H, eps0, eps);
            for(i=0; i<p*N; i++) y[i] = eps[i];
        }
        if(Zd) {
            Z   += p*m;
            /*** y(t) = Z(t)alpha(t) + eps(t) ***/
            dgemm("N", "N", &p, &N, &m, &g_one, Z, &p, alphaprev, &m, &g_one, y, &p);
        }
    }
    if(!Zd) {
        int N2  = N*n;
        y       = y_p0;
        alpha   = alpha_p0;
        /*** y = Zalpha + eps ***/
        dgemm("N", "N", &p, &N2, &m, &g_one, Z, &p, alpha, &m, &g_one, y, &p);
    }
}

/*===============================================*\
    Generate signal from state space components
\*===============================================*/
void signal(int p, int m, int n, int ncom, const int *mcom,
    const double *Z, int Zd, const double *alpha, double *ycom)
{
    int i;
    
    if(Zd) {
        int t, j = 1;
        
        for(t=0; t<n; t++) {
            double *ycomc   = ycom + t*p;
            
            for(i=0; i<ncom; i++) {
                dgemv("N", &p, mcom + i, &g_one, Z, &p, alpha, &j, &g_zero, ycomc, &j);
                
                Z       += mcom[i]*p;
                alpha   += mcom[i];
                ycomc   += p*n;
            }
        }
    }
    else {
        for(i=0; i<ncom; i++) {
            dgemm("N", "N", &p, &n, mcom + i, &g_one, Z, &p, alpha, &m, &g_zero, ycom, &p);
            
            Z       += mcom[i]*p;
            alpha   += mcom[i];
            ycom    += p*n;
        }
    }
}
