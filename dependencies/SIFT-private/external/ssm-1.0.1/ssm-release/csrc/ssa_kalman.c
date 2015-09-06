/**************************************\
    Filename: ssa_kalman.c
    Author: Jyh-Ying Peng (´^´¼·­)
    Year: 2006-2007
    Version: 1.4
\**************************************/

#include "ssa.h"

/*====================*\
    Global variables
\*====================*/
const int true          = 1;
const int false         = 0;
const double g_zero     = 0;
const double g_one      = 1;
const double g_minusone = -1;

/*==================================*\
    Matrix inverse and determinant
\*==================================*/
void invdet(int p, double *F, double *detF, int method)
{
    if(p == 1) {
        if(detF) *detF  = F[0];
        F[0]    = 1/F[0];
    }
    else {
        int i, j, info;
        
        switch(method) {
            case 0: {
                /*** General matrix ***/
                int *ipiv   = (int *)calloc(p, sizeof(int));
                int lwork   = -1;
                double dlwork, *work;
                
                dgetrf(&p, &p, F, &p, ipiv, &info);
                if(detF) {
                    *detF   = 1;
                    for(i=0; i<p; i++) {
                        *detF   *= F[i*p+i];
                        if(ipiv[i] != i+1) *detF = -(*detF);
                    }
                }
                dgetri(&p, F, &p, ipiv, &dlwork, &lwork, &info);
                lwork   = (int)dlwork;
                work    = (double *)calloc(lwork, sizeof(double));
                dgetri(&p, F, &p, ipiv, work, &lwork, &info);
                
                free(work);
                free(ipiv);
                break;
            }
            case 1: {
                /*** Symmetric matrix ***/
                int *ipiv   = (int *)calloc(p, sizeof(int));
                int lwork   = -1;
                double dlwork, *work;
                
                dsytrf("U", &p, F, &p, ipiv, &dlwork, &lwork, &info);
                lwork   = (int)dlwork;
                work    = (double *)calloc(lwork, sizeof(double));
                dsytrf("U", &p, F, &p, ipiv, work, &lwork, &info);
                if(detF) {
                    *detF   = 1;
                    for(i=0; i<p; i++) {
                        if(ipiv[i] > 0) *detF *= F[i*p+i];
                        else {
                            *detF *= F[i*p+i]*F[(i+1)*p+(i+1)] - F[(i+1)*p+i]*F[(i+1)*p+i];
                            i++;
                        }
                    }
                }
                dsytri("U", &p, F, &p, ipiv, work, &info);
                
                free(work);
                free(ipiv);
                
                for(i=1; i<p; i++) for(j=0; j<i; j++) F[i+j*p] = F[j+i*p];
                break;
            }
            case 2: {
                /*** Positive matrix ***/
                dpotrf("U", &p, F, &p, &info);
                if(detF) {
                    *detF   = 1;
                    for(i=0; i<p; i++) *detF *= F[i*p+i];
                    *detF *= *detF;
                }
                dpotri("U", &p, F, &p, &info);
                for(i=1; i<p; i++) for(j=0; j<i; j++) F[i+j*p] = F[j+i*p];
                break;
            }
        }
    }
}

/*==================================*\
    Normal Kalman filter iteration
\*==================================*/
void kalman_iter_normal(int p, int m,
        const double *H, const double *Z, const double *T, const double *RQRt, const double *Pprev,
        double *detF, double *F, double *invF, double *K, double *L, double *P,
        int inv_method)
{
    int i;
    double *M   = (double *)calloc(m*p, sizeof(double));
    double *TM  = (double *)calloc(m*p, sizeof(double));
    double *TP  = (double *)calloc(m*m, sizeof(double));
    
    /*** M = PZ' ***/
    dgemm("N", "T", &m, &p, &m, &g_one, Pprev, &m, Z, &p, &g_zero, M, &m);
    
    /*** invF = inv(ZM + H) ***/
#ifdef DONT_USE_LAPACK_DGEMM_ADD
    dgemm("N", "N", &p, &p, &m, &g_one, Z, &p, M, &m, &g_zero, invF, &p);
    if(F) for(i=0; i<p*p; i++) F[i] = invF[i] = invF[i] + H[i]; else for(i=0; i<p*p; i++) invF[i] += H[i];
#else
    for(i=0; i<p*p; i++) invF[i] = H[i];
    dgemm("N", "N", &p, &p, &m, &g_one, Z, &p, M, &m, &g_one, invF, &p);
    if(F) for(i=0; i<p*p; i++) F[i] = invF[i];
#endif
    invdet(p, invF, detF, inv_method);
    
    /*** K = TMinvF ***/
    dgemm("N", "N", &m, &p, &m, &g_one, T, &m, M, &m, &g_zero, TM, &m);
    dgemm("N", "N", &m, &p, &p, &g_one, TM, &m, invF, &p, &g_zero, K, &m);
    
    /*** L = T - K*Z ***/
#ifdef DONT_USE_LAPACK_DGEMM_ADD
    dgemm("N", "N", &m, &m, &p, &g_minusone, K, &m, Z, &p, &g_zero, L, &m);
    for(i=0; i<m*m; i++) L[i] += T[i];
#else
    for(i=0; i<m*m; i++) L[i] = T[i];
    dgemm("N", "N", &m, &m, &p, &g_minusone, K, &m, Z, &p, &g_one, L, &m);
#endif
    
    /*** P = TPL' + RQRt ***/
#ifdef DONT_USE_LAPACK_DGEMM_ADD
    dgemm("N", "N", &m, &m, &m, &g_one, T, &m, Pprev, &m, &g_zero, TP, &m);
    dgemm("N", "T", &m, &m, &m, &g_one, TP, &m, L, &m, &g_zero, P, &m);
    for(i=0; i<m*m; i++) P[i] += RQRt[i];
#else
    for(i=0; i<m*m; i++) P[i] = RQRt[i];
    dgemm("N", "N", &m, &m, &m, &g_one, T, &m, Pprev, &m, &g_zero, TP, &m);
    dgemm("N", "T", &m, &m, &m, &g_one, TP, &m, L, &m, &g_one, P, &m);
#endif

    free(M);
    free(TM);
    free(TP);
}

/*=========================================*\
    Exact diffuse Kalman filter iteration
\*=========================================*/
void kalman_iter_init(int p, int m,
        const double *H, const double *Z, const double *T, const double *RQRt,
        const double *Pprev, const double *P_infprev,
        int *Fns, double *detF, double *F, double *invF, double *F2,
        double *K, double *K1, double *L, double *L1, double *P, double *P_inf,
        double tol0, int inv_method)
{
    int i;
    double *M_inf   = (double *)calloc(m*p, sizeof(double));
    double *TP      = (double *)calloc(m*m, sizeof(double));
    
    /*** M_inf = P_infZ' ***/
    dgemm("N", "T", &m, &p, &m, &g_one, P_infprev, &m, Z, &p, &g_zero, M_inf, &m);
    
    /*** Is F_inf zero? ***/
    *Fns    = false;
    for(i=0; i<m*p; i++)
        if(fabs(M_inf[i]) >= tol0) {
            *Fns    = true;
            break;
        }
    
    /*** F_inf nonsingular branch ***/
    if(*Fns) {
        /*** Temporary variables ***/
        double *M   = (double *)calloc(m*p, sizeof(double));
        double *F1F = (double *)calloc(p*p, sizeof(double));
        double *TM  = (double *)calloc(m*p, sizeof(double));
        double *MF  = TM;
#ifdef DONT_USE_LAPACK_DGEMM_ADD
        double *TPL = (double *)calloc(m*m, sizeof(double));
#endif
        /*** invF = inv(ZM_inf) (actually F1) ***/
        dgemm("N", "N", &p, &p, &m, &g_one, Z, &p, M_inf, &m, &g_zero, invF, &p);
        invdet(p, invF, detF, inv_method);
        
        /*** M = PZ' ***/
        dgemm("N", "T", &m, &p, &m, &g_one, Pprev, &m, Z, &p, &g_zero, M, &m);

        /*** F2 = -invF(ZM + H)invF ***/
        if(!F) F = F2; /* Borrow memory space from output */
#ifdef DONT_USE_LAPACK_DGEMM_ADD
        dgemm("N", "N", &p, &p, &m, &g_one, Z, &p, M, &m, &g_zero, F, &p);
        for(i=0; i<p*p; i++) F[i] += H[i];
#else
        for(i=0; i<p*p; i++) F[i] = H[i];
        dgemm("N", "N", &p, &p, &m, &g_one, Z, &p, M, &m, &g_one, F, &p);
#endif
        dgemm("N", "N", &p, &p, &p, &g_one, invF, &p, F, &p, &g_zero, F1F, &p);
        dgemm("N", "N", &p, &p, &p, &g_minusone, F1F, &p, invF, &p, &g_zero, F2, &p);
        
        /*** K = TM_infinvF ***/
        dgemm("N", "N", &m, &p, &m, &g_one, T, &m, M_inf, &m, &g_zero, TM, &m);
        dgemm("N", "N", &m, &p, &p, &g_one, TM, &m, invF, &p, &g_zero, K, &m);
        
        /*** K1 = T(MinvF + M_infF2) ***/
#ifdef DONT_USE_LAPACK_DGEMM_ADD
        dgemm("N", "N", &m, &p, &p, &g_one, M, &m, invF, &p, &g_zero, MF, &m);
        dgemm("N", "N", &m, &p, &p, &g_one, M_inf, &m, F2, &p, &g_zero, M, &m);
        for(i=0; i<m*p; i++) MF[i] += M[i];
#else
        dgemm("N", "N", &m, &p, &p, &g_one, M, &m, invF, &p, &g_zero, MF, &m);
        dgemm("N", "N", &m, &p, &p, &g_one, M_inf, &m, F2, &p, &g_one, MF, &m);
#endif
        dgemm("N", "N", &m, &p, &m, &g_one, T, &m, MF, &m, &g_zero, K1, &m);
        
        /*** L = T - KZ ***/
#ifdef DONT_USE_LAPACK_DGEMM_ADD
        dgemm("N", "N", &m, &m, &p, &g_minusone, K, &m, Z, &p, &g_zero, L, &m);
        for(i=0; i<m*m; i++) L[i] += T[i];
#else
        for(i=0; i<m*m; i++) L[i] = T[i];
        dgemm("N", "N", &m, &m, &p, &g_minusone, K, &m, Z, &p, &g_one, L, &m);
#endif
        /*** L1 = -K1Z ***/
        dgemm("N", "N", &m, &m, &p, &g_minusone, K1, &m, Z, &p, &g_zero, L1, &m);
        
        /*** P = TPL' + TP_infL1' + RQRt ***/
#ifdef DONT_USE_LAPACK_DGEMM_ADD
        dgemm("N", "N", &m, &m, &m, &g_one, T, &m, Pprev, &m, &g_zero, TP, &m);
        dgemm("N", "T", &m, &m, &m, &g_one, TP, &m, L, &m, &g_zero, P, &m);
        dgemm("N", "N", &m, &m, &m, &g_one, T, &m, P_infprev, &m, &g_zero, TP, &m);
        dgemm("N", "T", &m, &m, &m, &g_one, TP, &m, L1, &m, &g_zero, TPL, &m);
        for(i=0; i<m*m; i++) P[i] = P[i] + TPL[i] + RQRt[i];
#else
        for(i=0; i<m*m; i++) P[i] = RQRt[i];
        dgemm("N", "N", &m, &m, &m, &g_one, T, &m, Pprev, &m, &g_zero, TP, &m);
        dgemm("N", "T", &m, &m, &m, &g_one, TP, &m, L, &m, &g_one, P, &m);
        dgemm("N", "N", &m, &m, &m, &g_one, T, &m, P_infprev, &m, &g_zero, TP, &m);
        dgemm("N", "T", &m, &m, &m, &g_one, TP, &m, L1, &m, &g_one, P, &m);
#endif
        /*** P_inf = TP_infL' ***/
        dgemm("N", "T", &m, &m, &m, &g_one, TP, &m, L, &m, &g_zero, P_inf, &m);

        free(M);
        free(F1F);
        free(TM);
#ifdef DONT_USE_LAPACK_DGEMM_ADD
        free(TPL);
#endif
    }
    else {
        kalman_iter_normal(p, m, H, Z, T, RQRt, Pprev, detF, F, invF, K, L, P, inv_method);
        
        /*** F2 = 0 ***/
        for(i=0; i<p*p; i++) F2[i] = 0;
        
        /*** K1 = 0 ***/
        for(i=0; i<m*p; i++) K1[i] = 0;
        
        /*** L1 = 0 ***/
        for(i=0; i<m*m; i++) L1[i] = 0;
        
        /*** P_inf = TP_infT' ***/
        dgemm("N", "N", &m, &m, &m, &g_one, T, &m, P_infprev, &m, &g_zero, TP, &m);
        dgemm("N", "T", &m, &m, &m, &g_one, TP, &m, T, &m, &g_zero, P_inf, &m);
    }
    
    free(M_inf);
    free(TP);
}

/*========================================*\
    Missing data Kalman filter iteration
\*========================================*/
void kalman_iter_miss(int p, int m,
        const double *H, const double *Z, const double *T, const double *RQRt,
        const double *Pprev, const double *P_infprev,
        double *detF, double *F, double *invF, double *K, double *L, double *P, double *P_inf,
        int inv_method)
{
    /******* TODO: Does K, L have any meaning within this context? *******/
    
    /*** Temporary variables ***/
    int i;
    double *M   = (double *)calloc(m*p, sizeof(double));
    double *TP  = (double *)calloc(m*m, sizeof(double));
    
    /*** M = PZ' ***/
    dgemm("N", "T", &m, &p, &m, &g_one, Pprev, &m, Z, &p, &g_zero, M, &m);
    
    /*** invF = inv(ZM + H) ***/
#ifdef DONT_USE_LAPACK_DGEMM_ADD
    dgemm("N", "N", &p, &p, &m, &g_one, Z, &p, M, &m, &g_zero, invF, &p);
    if(F) for(i=0; i<p*p; i++) F[i] = invF[i] = invF[i] + H[i]; else for(i=0; i<p*p; i++) invF[i] += H[i];
#else
    for(i=0; i<p*p; i++) invF[i] = H[i];
    dgemm("N", "N", &p, &p, &m, &g_one, Z, &p, M, &m, &g_one, invF, &p);
    if(F) for(i=0; i<p*p; i++) F[i] = invF[i];
#endif
    invdet(p, invF, detF, inv_method);
    
    /*** K = 0 ***/
    for(i=0; i<m*p; i++) K[i] = 0;
    
    /*** L = T ***/
    for(i=0; i<m*m; i++) L[i] = T[i];
    
    /*** P = TPT' + RQRt ***/
#ifdef DONT_USE_LAPACK_DGEMM_ADD
    dgemm("N", "N", &m, &m, &m, &g_one, T, &m, Pprev, &m, &g_zero, TP, &m);
    dgemm("N", "T", &m, &m, &m, &g_one, TP, &m, T, &m, &g_zero, P, &m);
    for(i=0; i<m*m; i++) P[i] += RQRt[i];
#else
    for(i=0; i<m*m; i++) P[i] = RQRt[i];
    dgemm("N", "N", &m, &m, &m, &g_one, T, &m, Pprev, &m, &g_zero, TP, &m);
    dgemm("N", "T", &m, &m, &m, &g_one, TP, &m, T, &m, &g_one, P, &m);
#endif
    /*** Exact diffuse initialization branch ***/
    if(P_infprev) {
        /*** P_inf = TP_infT' ***/
        dgemm("N", "N", &m, &m, &m, &g_one, T, &m, P_infprev, &m, &g_zero, TP, &m);
        dgemm("N", "T", &m, &m, &m, &g_one, TP, &m, T, &m, &g_zero, P_inf, &m);
    }

    free(M);
    free(TP);
}

/*================================*\
    Kalman filter data iteration
\*================================*/
void kalman_iter_data(int p, int m, int N, const double *y,
        const double *Z, const double *T, const double *c, const double *K, const double *aprev,
        double *v, double *a)
{
    /*** Temporary variables ***/
    int i;
    
    if(y) {
#ifndef USE_LAPACK_DGEMM_ADD
        double *Kv  = (double *)calloc(m*N, sizeof(double));
#endif
        /*** v = y - Za ***/
        dgemm("N", "N", &p, &N, &m, &g_minusone, Z, &p, aprev, &m, &g_zero, v, &p);
        for(i=0; i<p*N; i++) v[i] += y[i];

        /*** a = c + Ta + Kv ***/
#ifdef USE_LAPACK_DGEMM_ADD
        for(i=0; i<m*N; i++) a[i] = c[i%m];
        dgemm("N", "N", &m, &N, &m, &g_one, T, &m, aprev, &m, &g_one, a, &m);
        dgemm("N", "N", &m, &N, &p, &g_one, K, &m, v, &p, &g_one, a, &m);
#else
        dgemm("N", "N", &m, &N, &p, &g_one, K, &m, v, &p, &g_zero, Kv, &m);
        dgemm("N", "N", &m, &N, &m, &g_one, T, &m, aprev, &m, &g_zero, a, &m);
        for(i=0; i<m*N; i++) a[i] = c[i%m] + a[i] + Kv[i];
        free(Kv);
#endif
    }
    else {
        /*** Missing data ***/
        /*** v = 0 ***/
        for(i=0; i<p*N; i++) v[i] = 0;
        
        /*** a = c + Ta ***/
        dgemm("N", "N", &m, &N, &m, &g_one, T, &m, aprev, &m, &g_zero, a, &m);
        for(i=0; i<m*N; i++) a[i] += c[i%m];
    }
}

/*=================*\
    Kalman filter
\*=================*/
void kalman(int p, int m, int r, int n, int N, const double *y, const int *pknown, int pd, const int *pmiss,
        int sta, const double *H, int Hd, const double *Z, int Zd, const double *T, int Td, const double *R, int Rd,
        const double *Q, int Qd, const double *c, int cd, const double *a1, const double *P1, const double *P1_inf,
        double *a, double *P, double *P_inf, int *d, int *Fns, double *v, double *F, double *invF, double *F2,
        double *K, double *K1, double *L, double *L1, double *snlogL, double *fvar, double *RQ, double *QRt, double *RQRt,
        double tol0, double tolP, int inv_method)
{
    /*** Kalman filter information ***/
    int miss    = pmiss != 0;
    int init    = P1_inf != 0;
    int conv    = false;
    int conv2   = false; /* for when P is not outputed */
    int RQd     = Rd || Qd;

    /*** Temporary variables ***/
    int i, j, t, pcur;
    int d_, Fns_;
    double detF_;
    double *detF    = snlogL ? &detF_ : (double *)0;
    double *aprev, *Pprev, *P_infprev;
    double *invFv;
    
    /*** "Allocate" buffer variables ***/
    int aout        = (a != 0);
    int Pout        = (P != 0);
    int P_infout    = (P_inf != 0);
    int vout        = (v != 0);
    int Fout        = (F != 0);
    int invFout     = (invF != 0);
    int F2out       = (F2 != 0);
    int Kout        = (K != 0);
    int K1out       = (K1 != 0);
    int Lout        = (L != 0);
    int L1out       = (L1 != 0);
    int RQout       = (RQ != 0);
    int QRtout      = (QRt != 0);
    int RQRtout     = (RQRt != 0);
    
    if(!aout) { a = (double *)calloc(m*N, sizeof(double)); aprev = (double *)calloc(m*N, sizeof(double)); } else { aprev = a; a += m*N; }
    if(!Pout) { P = (double *)calloc(m*m, sizeof(double)); Pprev = (double *)calloc(m*m, sizeof(double)); } else { Pprev = P; P += m*m; }
    if(!vout) v = (double *)calloc(p*N, sizeof(double));
    if(!invFout) invF = (double *)calloc(p*p, sizeof(double));
    if(!Kout) K = (double *)calloc(m*p, sizeof(double));
    if(!Lout) L = (double *)calloc(m*m, sizeof(double));
    if(!RQRtout) RQRt = (double *)calloc(m*m, sizeof(double));
    if(snlogL || fvar) {
        invFv = (double *)calloc(p*N, sizeof(double));
        if(snlogL) for(i=0; i<N; i++) snlogL[i] = 0;
        if(fvar) for(i=0; i<N; i++) fvar[i] = 0;
    }
    if(!QRtout) QRt = (double *)calloc(r*m, sizeof(double));
    
    /*** Initial state ***/
    if(d == 0) d = &d_; d[0] = init ? n+1 : 0;
    for(i=0; i<m*N; i++) aprev[i] = a1[i%m];
    for(i=0; i<m*m; i++) Pprev[i] = P1[i];
    if(init) {
        if(!P_infout) { P_inf = (double *)calloc(m*m, sizeof(double)); P_infprev = (double *)calloc(m*m, sizeof(double)); } else { P_infprev = P_inf; P_inf += m*m; }
        if(!F2out) F2 = (double *)calloc(p*p, sizeof(double));
        if(!K1out) K1 = (double *)calloc(m*p, sizeof(double));
        if(!L1out) L1 = (double *)calloc(m*m, sizeof(double));
        for(i=0; i<m*m; i++) P_infprev[i] = P1_inf[i];
    }

    /*** Main loop ***/
    for(t=0; t<n; t++) {
        /*** Calculate RQRt ***/
        if(RQd || t==0) {
            dgemm("N", "T", &r, &m, &r, &g_one, Q, &r, R, &m, &g_zero, QRt, &r);
            dgemm("N", "N", &m, &m, &r, &g_one, R, &m, QRt, &r, &g_zero, RQRt, &m);
            if(RQout) dgemm("N", "N", &m, &r, &r, &g_one, R, &m, Q, &r, &g_zero, RQ, &m);
        }

        /*** Get data dimension ***/
        pcur    = pd ? pknown[t] : p;
        
        /*** Missing data throws off convergence ***/
        if(conv && miss && pmiss[t] < p) conv = false;
        
        /*** Kalman filter ***/
        if(miss && pmiss[t] == 0) {
            /*** Missing data case ***/
            if(init) {
                kalman_iter_miss(pcur, m, H, Z, T, RQRt, Pprev, P_infprev, (double *)0, F, invF, K, L, P, P_inf, inv_method);
                Fns_    = false;
                if(F2out) for(i=0; i<p*p; i++) F2[i] = 0;
                if(K1out) for(i=0; i<p*m; i++) K1[i] = 0;
                if(L1out) for(i=0; i<m*m; i++) L1[i] = 0;
            }
            else
                kalman_iter_miss(pcur, m, H, Z, T, RQRt, Pprev, (double *)0, (double *)0, F, invF, K, L, P, (double *)0, inv_method);
            kalman_iter_data(pcur, m, N, (double *)0, Z, T, c, K, aprev, v, a);
        }
        else if(conv) {
            /*** Converged case ***/
            if(Fout) for(i=0; i<p*p; i++) F[i] = F[i-p*p];
            if(invFout) for(i=0; i<p*p; i++) invF[i] = invF[i-p*p];
            if(Kout) for(i=0; i<p*m; i++) K[i] = K[i-p*m];
            if(Lout) for(i=0; i<m*m; i++) L[i] = L[i-m*m];
            if(Pout) for(i=0; i<m*m; i++) P[i] = Pprev[i];
            conv2   = true; /* Pprev are now equal to the converged value, do not swap P, Pprev */
            kalman_iter_data(p, m, N, y, Z, T, c, K, aprev, v, a);
        }
        else {
            if(miss && pmiss[t] < pcur) pcur = pmiss[t];
            if(init) {
                /*** Exact diffuse case ***/
                kalman_iter_init(pcur, m, H, Z, T, RQRt, Pprev, P_infprev, &Fns_, detF, F, invF, F2, K, K1, L, L1, P, P_inf,
                    tol0, inv_method);
                /*** P_inf = 0? ***/
                init    = false;
                for(i=0; i<m*m; i++)
                    if(fabs(P_inf[i]) >= tol0) {
                        init    = true;
                        break;
                    }
                if(!init) d[0] = t + 1;
            }
            else {
                /*** Normal case ***/
                kalman_iter_normal(pcur, m, H, Z, T, RQRt, Pprev, detF, F, invF, K, L, P, inv_method);
                /*** Filter converged? ***/
                if(sta) {
                    conv    = true;
                    for(i=0; i<m*m; i++)
                        if(fabs(P[i] - Pprev[i]) >= tolP) {
                            conv = false;
                            break;
                        }
                }
            }
            kalman_iter_data(pcur, m, N, y, Z, T, c, K, aprev, v, a);
        }
        
        /*** Loglikelihood Calculations ***/
        if((snlogL || fvar) && (!miss || pmiss[t] > 0)) {
            if(t >= d[0] || !Fns_) {
                dgemm("N", "N", &pcur, &N, &pcur, &g_one, invF, &pcur, v, &pcur, &g_zero, invFv, &pcur);
                for(i=0; i<N; i++)
                    for(j=0; j<pcur; j++) {
                        double vinvFv   = v[i*pcur+j]*invFv[i*pcur+j];
                        if(fvar) fvar[i] += vinvFv;
                        if(snlogL) snlogL[i] += vinvFv;
                    }
            }
            if(snlogL) {
                double logdetF  = log(detF[0]);
                for(i=0; i<N; i++) snlogL[i] += logdetF;
            }
        }

        /*** Goto next observation ***/
        y   += p*N;
        /*** Goto next state matrices ***/
        if(Hd) H += p*p;
        if(Zd) Z += p*m;
        if(Td) T += m*m;
        if(Rd) R += m*r;
        if(Qd) Q += r*r;
        if(cd) c += m;
        /*** Goto next output matrices ***/
        if(aout) { aprev = a; a += m*N; } else { double *atemp = aprev; aprev = a; a = atemp; }
        if(Pout) { Pprev = P; P += m*m; } else if(!conv2) { double *Ptemp = Pprev; Pprev = P; P = Ptemp; }
        if(Fns) { if(t < d[0]) Fns[t] = Fns_; else Fns[t] = false; }
        if(vout) v += p*N;
        if(Fout) F += p*p;
        if(invFout) invF += p*p;
        if(Kout) K += m*p;
        if(Lout) L += m*m;
        if(RQout && RQd) RQ += m*r;
        if(QRtout && RQd) QRt += r*m;
        if(RQRtout && RQd) RQRt += m*m;
        if(init) {
            if(P_infout) { P_infprev = P_inf; P_inf += m*m; } else { double *Ptemp = P_infprev; P_infprev = P_inf; P_inf = Ptemp; }
            if(F2out) F2 += p*p;
            if(K1out) K1 += m*p;
            if(L1out) L1 += m*m;
        }
    }

    /*** Cleanup ***/
    if(!aout) { free(a); free(aprev); }
    if(!Pout) { free(P); free(Pprev); }
    if(!vout) free(v);
    if(!invFout) free(invF);
    if(!Kout) free(K);
    if(!Lout) free(L);
    if(!RQRtout) free(RQRt);
    if(snlogL || fvar) free(invFv);
    if(!QRtout) free(QRt);
    if(P1_inf) {
        if(!P_infout) { free(P_inf); free(P_infprev); }
        if(!F2out) free(F2);
        if(!K1out) free(K1);
        if(!L1out) free(L1);
    }
}

/*===================================*\
    Filtering and smoothing weights
\*===================================*/
void weights(int p, int m, int r, int n, int t0, const int *pknown, int pd,
        int sta, const double *H, int Hd, const double *Z, int Zd, const double *T, int Td,
        const double *R, int Rd, const double *Q, int Qd, const double *P1, const double *P1_inf,
        double *wt_a, double *wt_alpha, double tol0, double tolP, int inv_method)
{
    /******* TODO: Look into exact diffuse initialization for wt_alpha?! Requires O(d^2) LLasc and LLdesc?! *******/
    
    /*** Kalman filter information ***/
    int init    = P1_inf != 0;
    int conv    = false;
    int conv2   = false;
    int RQd     = Rd || Qd;

    /*** Temporary variables ***/
    int i, j, t, pcur;
    int d, Fns_, *Fns;
    double *P, *Pprev, *P_inf, *P_infprev, *Pt0;
    double *invF, *invF_p0, *F2, *K, *K_p0, *K1, *L, *L_p0, *L1, *QRt, *RQRt;
    double *LLasc, *LLasc_p0, *LLdsc, *Ltemp;
    
    /*** Allocate buffer variables ***/
    int alphaout    = (wt_alpha != 0);
    
    t0          = t0 - 1; /* Adjust to 0-based index */
    if(t0 == 0)
        /*** a(1) = a1 does not depend on y(t) at all!!! ***/
        for(i=0; i<m*p*n; i++) wt_a[i] = 0;
    if(alphaout) {
        Pt0     = (double *)calloc(m*m, sizeof(double));
        if(t0 == 0) for(i=0; i<m*m; i++) Pt0[i] = P1[i];
        invF    = invF_p0 = (double *)calloc((n > t0) ? p*p*(n-t0) : p*p, sizeof(double)); /* Only invF(t0:n-1) is needed */
        K       = K_p0 = (double *)calloc(m*p*n, sizeof(double));
        L       = L_p0 = (double *)calloc(m*m*n, sizeof(double));
        if(n > t0+1) LLasc = LLasc_p0 = (double *)calloc(m*m*(n-t0-1), sizeof(double));
        else LLasc = LLasc_p0 = (double *)0;
        Fns     = (int *)calloc(n, sizeof(int));
    }
    else if(t0 > 0) {
        invF    = invF_p0 = (double *)calloc(p*p, sizeof(double));
        K       = K_p0 = (double *)calloc(m*p*t0, sizeof(double)); /* Only K(0:t0) is needed */
        L       = L_p0 = (double *)calloc(m*m*t0, sizeof(double)); /* Only L(0:t0) is needed */
        Fns     = &Fns_;
    }
    else return;
    P           = (double *)calloc(m*m, sizeof(double));
    Pprev       = (double *)calloc(m*m, sizeof(double));
    RQRt        = (double *)calloc(m*m, sizeof(double));
    QRt         = (double *)calloc(r*m, sizeof(double));
    
    /*** Initial state ***/
    d = init ? n+1 : 0;
    for(i=0; i<m*m; i++) Pprev[i] = P1[i];
    if(init) {
        P_inf       = (double *)calloc(m*m, sizeof(double));
        P_infprev   = (double *)calloc(m*m, sizeof(double));
        F2          = (double *)calloc(p*p, sizeof(double));
        K1          = (double *)calloc(m*p, sizeof(double));
        L1          = (double *)calloc(m*m, sizeof(double));
        for(i=0; i<m*m; i++) P_infprev[i] = P1_inf[i];
    }

    /*** Main loop ***/
    for(t=0; t<n; t++) {
        /*** Calculate RQRt ***/
        if(RQd || t==0) {
            dgemm("N", "T", &r, &m, &r, &g_one, Q, &r, R, &m, &g_zero, QRt, &r);
            dgemm("N", "N", &m, &m, &r, &g_one, R, &m, QRt, &r, &g_zero, RQRt, &m);
        }

        /*** Get data dimension ***/
        pcur    = pd ? pknown[t] : p;
        
        /*** Kalman filter ***/
        if(conv) {
            /*** Converged case ***/
            if(t > t0) for(i=0; i<p*p; i++) invF[i] = invF[i-p*p];
            for(i=0; i<p*m; i++) K[i] = K[i-p*m];
            for(i=0; i<m*m; i++) L[i] = L[i-m*m];
            conv2   = true;
        }
        else {
            if(init) {
                /*** Exact diffuse case ***/
                kalman_iter_init(pcur, m, H, Z, T, RQRt, Pprev, P_infprev, alphaout ? Fns + t : Fns, (double *)0, (double *)0, invF, F2, K, K1, L, L1, P, P_inf,
                    tol0, inv_method);
                /*** P_inf = 0? ***/
                init    = false;
                for(i=0; i<m*m; i++)
                    if(fabs(P_inf[i]) >= tol0) {
                        init    = true;
                        break;
                    }
                if(!init) d = t + 1;
            }
            else {
                /*** Normal case ***/
                kalman_iter_normal(pcur, m, H, Z, T, RQRt, Pprev, (double *)0, (double *)0, invF, K, L, P, inv_method);
                /*** Filter converged? ***/
                if(sta) {
                    conv    = true;
                    for(i=0; i<m*m; i++)
                        if(fabs(P[i] - Pprev[i]) >= tolP) {
                            conv = false;
                            break;
                        }
                }
            }
        }
        
        /*** Leave loop when all data is obtained ***/
        if((!alphaout && t == t0-1) || (t == n-1)) break;
        /*** Goto next state matrices ***/
        if(Hd) H += p*p;
        if(Zd) Z += p*m;
        if(Td) T += m*m;
        if(Rd) R += m*r;
        if(Qd) Q += r*r;
        /*** Goto next output matrices ***/
        if(!conv2) { double *Ptemp = Pprev; Pprev = P; P = Ptemp; }
        if(t == t0-1) for(i=0; i<m*m; i++) Pt0[i] = Pprev[i]; /* Here Pprev = P(t0) */
        if(t >= t0) {
            invF    += p*p;
            if(t == t0) for(i=0; i<m*m; i++) LLasc[i] = L[i];
            else dgemm("N", "N", &m, &m, &m, &g_one, L, &m, LLasc - m*m, &m, &g_zero, LLasc, &m);
            if(t < n-2) LLasc += m*m;
        }
        K += m*p;
        L += m*m;
        if(init) { double *Ptemp = P_infprev; P_infprev = P_inf; P_inf = Ptemp; }
    }
    
    free(P);
    free(Pprev);
    free(RQRt);
    free(QRt);
    if(P1_inf) {
        free(P_inf);
        free(P_infprev);
        free(F2);
        free(K1);
        free(L1);
    }

    LLdsc   = (double *)calloc(m*m, sizeof(double));
    Ltemp   = (double *)calloc(m*m, sizeof(double));
    
    if(alphaout) {
        double *N       = (double *)calloc(m*m, sizeof(double));
        double *Nprev   = (double *)calloc(m*m, sizeof(double));
        double *IPN     = Nprev;
        double *KtN     = (double *)calloc(p*m, sizeof(double));
        double *W       = (double *)calloc(p*m, sizeof(double));
        double *LtWt    = KtN;
        
        wt_a        += m*p*(t0-1);
        wt_alpha    += m*p*(n-1);
        
        if(t0 == n) for(i=0; i<m; i++) for(j=0; j<m; j++) IPN[i*m+j] = LLdsc[i*m+j] = (i==j) ? 1 : 0;
        else for(i=0; i<m*m; i++) Nprev[i] = 0;

        for(t=n-1; t>=0; t--) {
            /*** Get data dimension ***/
            pcur    = pd ? pknown[t] : p;
            
            if(t >= t0) {
                if(t < d && Fns[t]) {
                    /*** W(t) = -K(t)'N(t)L(t) ***/
                    dgemm("T", "N", &pcur, &m, &m, &g_one, K, &m, Nprev, &m, &g_zero, KtN, &pcur);
                    dgemm("N", "N", &pcur, &m, &m, &g_minusone, KtN, &pcur, L, &m, &g_zero, W, &pcur);
                    
                    backward_iter_init(pcur, m, 1, Z, (double *)0, true, (double *)0, invF, (double *)0, L, (double *)0,
                        (double *)0, (double *)0, Nprev, (double *)0, (double *)0, (double *)0, (double *)0, N, (double *)0, (double *)0);
                }
                else {
                    /*** W(t) = invF(t)Z(t) - K(t)'N(t)L(t) ***/
                    dgemm("N", "N", &pcur, &m, &pcur, &g_one, invF, &pcur, Z, &pcur, &g_zero, W, &pcur);
                    dgemm("T", "N", &pcur, &m, &m, &g_one, K, &m, Nprev, &m, &g_zero, KtN, &pcur);
                    dgemm("N", "N", &pcur, &m, &m, &g_minusone, KtN, &pcur, L, &m, &g_one, W, &pcur);
                    
                    backward_iter_normal(pcur, m, 1, Z, (double *)0, invF, L, (double *)0, Nprev, (double *)0, N);
                }
                if(t == t0) {
                    /*** wt_alpha(t) = P(t0)W(t0)' ***/
                    dgemm("N", "T", &m, &pcur, &m, &g_one, Pt0, &m, W, &pcur, &g_zero, wt_alpha, &m);
                    
                    /*** IPN = I - P(t0)N(t0-1) ***/
                    IPN     = Nprev; /* Borrow memory space (Nprev will not be used after) */
                    for(i=0; i<m; i++) for(j=0; j<m; j++) IPN[i*m+j] = LLdsc[i*m+j] = (i==j) ? 1 : 0;
                    dgemm("N", "N", &m, &m, &m, &g_minusone, Pt0, &m, N, &m, &g_one, IPN, &m);
                }
                else {
                    double *Ntemp;

                    /*** wt_alpha(t) = P(t0)L(t0)'L(t0+1)' ... L(t-1)'W(t)' ***/
                    dgemm("T", "T", &m, &pcur, &m, &g_one, LLasc, &m, W, &pcur, &g_zero, LtWt, &m);
                    dgemm("N", "N", &m, &pcur, &m, &g_one, Pt0, &m, LtWt, &m, &g_zero, wt_alpha, &m);
                    Ntemp = Nprev; Nprev = N; N = Ntemp;
                    if(Zd) Z -= p*m;
                    invF    -= p*p;
                    LLasc   -= m*m;
                }
            }
            else {
                /*** wt_a(t) = L(t0-1) ... L(t+1)K(t) ***/
                dgemm("N", "N", &m, &pcur, &m, &g_one, LLdsc, &m, K, &m, &g_zero, wt_a, &m);
                dgemm("N", "N", &m, &m, &m, &g_one, LLdsc, &m, L, &m, &g_zero, Ltemp, &m);
                for(i=0; i<m*m; i++) LLdsc[i] = Ltemp[i];
                
                /*** wt_alpha(t) = (I - P(t0)N(t0-1))L(t0-1) ... L(t+1)K(t) ***/
                dgemm("N", "N", &m, &pcur, &m, &g_one, IPN, &m, wt_a, &m, &g_zero, wt_alpha, &m);
                
                wt_a    -= m*p;
            }
            wt_alpha    -= m*p;
            K           -= m*p;
            L           -= m*m;
        }
        
        free(N);
        free(Nprev);
        free(KtN);
        free(W);
    }
    else {
        wt_a    += m*p*(t0-1);
        
        for(i=0; i<m; i++) for(j=0; j<m; j++) LLdsc[i*m+j] = (i==j) ? 1 : 0;
        
        for(t=t0-1; t>=0; t--) {
            /*** Get data dimension ***/
            pcur    = pd ? pknown[t] : p;
            
            /*** wt_a(t) = L(t0-1) ... L(t+1)K(t) ***/
            dgemm("N", "N", &m, &pcur, &m, &g_one, LLdsc, &m, K, &m, &g_zero, wt_a, &m);
            dgemm("N", "N", &m, &m, &m, &g_one, LLdsc, &m, L, &m, &g_zero, Ltemp, &m);
            for(i=0; i<m*m; i++) LLdsc[i] = Ltemp[i];
            
            /*** Next ***/
            wt_a    -= m*p;
            K       -= m*p;
            L       -= m*m;
        }
    }
    
    /*** Cleanup ***/
    free(LLdsc);
    free(Ltemp);
    free(invF_p0);
    free(K_p0);
    free(L_p0);
    if(alphaout) {
        free(Pt0);
        if(n > t0+1) free(LLasc_p0);
        free(Fns);
    }
}
