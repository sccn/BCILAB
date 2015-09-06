/**************************************\
    Filename: ssa_smooth.c
    Author: Jyh-Ying Peng (´^´¼·­)
    Year: 2007
    Version: 1.2
\**************************************/

#include "ssa.h"

/*====================*\
    Global variables
\*====================*/
extern true;
extern false;
extern g_zero;
extern g_one;
extern g_minusone;

/*========================================*\
    Normal backward recursion iteration
\*========================================*/
void backward_iter_normal(int p, int m, int NN, const double *Z,
        const double *v, const double *invF, const double *L,
        const double *rprev, const double *Nprev, double *r, double *N)
{
    /*** Temporary variables ***/
    double *ZtF     = (double *)calloc(m*p, sizeof(double));
    
    dgemm("T", "N", &m, &p, &p, &g_one, Z, &p, invF, &p, &g_zero, ZtF, &m);

    if(rprev) {
        /*** r = Z'invFv + L'r ***/
#ifdef DONT_USE_LAPACK_DGEMM_ADD
        int i;
        double *ZtFv    = (double *)calloc(m*NN, sizeof(double));
        dgemm("N", "N", &m, &NN, &p, &g_one, ZtF, &m, v, &p, &g_zero, ZtFv, &m);
        dgemm("T", "N", &m, &NN, &m, &g_one, L, &m, rprev, &m, &g_zero, r, &m);
        for(i=0; i<m*NN; i++) r[i] += ZtFv[i];
        free(ZtFv);
#else
        dgemm("N", "N", &m, &NN, &p, &g_one, ZtF, &m, v, &p, &g_zero, r, &m);
        dgemm("T", "N", &m, &NN, &m, &g_one, L, &m, rprev, &m, &g_one, r, &m);
#endif
    }
    
    if(Nprev) {
        double *LtN     = (double *)calloc(m*m, sizeof(double));
#ifdef DONT_USE_LAPACK_DGEMM_ADD
        int i;
        double *LtNL    = (double *)calloc(m*m, sizeof(double));
#endif
        /*** N = Z'invFZ + L'NL ***/
        dgemm("N", "N", &m, &m, &p, &g_one, ZtF, &m, Z, &p, &g_zero, N, &m);
        dgemm("T", "N", &m, &m, &m, &g_one, L, &m, Nprev, &m, &g_zero, LtN, &m);
#ifdef DONT_USE_LAPACK_DGEMM_ADD
        dgemm("N", "N", &m, &m, &m, &g_one, LtN, &m, L, &m, &g_zero, LtNL, &m);
        for(i=0; i<m*m; i++) N[i] += LtNL[i];
        free(LtNL);
#else
        dgemm("N", "N", &m, &m, &m, &g_one, LtN, &m, L, &m, &g_one, N, &m);
#endif
        free(LtN);
    }
    
    free(ZtF);
}

/*==============================================*\
    Exact diffuse backward recursion iteration
\*==============================================*/
void backward_iter_init(int p, int m, int NN, const double *Z, const double *T,
        int Fns, const double *v, const double *invF, const double *F2,
        const double *L, const double *L1, const double *rprev, const double *r1prev,
        const double *Nprev, const double *N1prev, const double *N2prev,
        double *r, double *r1, double *N, double *N1, double *N2)
{
    if(Fns) {
        /*** Temporary variables ***/
        double *ZtF     = (double *)calloc(m*p, sizeof(double));
        double *LtN     = (double *)calloc(m*m, sizeof(double));
        
        if(r1prev) {
            /*** r1 = Z'invFv + L'r1 + L1'r ***/
#ifdef DONT_USE_LAPACK_DGEMM_ADD
            int i;
            double *ZtFv    = (double *)calloc(m*NN, sizeof(double));
            double *Ltr     = (double *)calloc(m*NN, sizeof(double));
            dgemm("T", "N", &m, &p, &p, &g_one, Z, &p, invF, &p, &g_zero, ZtF, &m);
            dgemm("N", "N", &m, &NN, &p, &g_one, ZtF, &m, v, &p, &g_zero, ZtFv, &m);
            dgemm("T", "N", &m, &NN, &m, &g_one, L, &m, r1prev, &m, &g_zero, Ltr, &m);
            dgemm("T", "N", &m, &NN, &m, &g_one, L1, &m, rprev, &m, &g_zero, r1, &m);
            for(i=0; i<m*NN; i++) r1[i] = ZtFv[i] + Ltr[i] + r1[i];
            free(ZtFv);
            free(Ltr);
#else
            dgemm("T", "N", &m, &p, &p, &g_one, Z, &p, invF, &p, &g_zero, ZtF, &m);
            dgemm("N", "N", &m, &NN, &p, &g_one, ZtF, &m, v, &p, &g_zero, r1, &m);
            dgemm("T", "N", &m, &NN, &m, &g_one, L, &m, r1prev, &m, &g_one, r1, &m);
            dgemm("T", "N", &m, &NN, &m, &g_one, L1, &m, rprev, &m, &g_one, r1, &m);
#endif
        }

        if(rprev) {
            /*** r = L'r ***/
            dgemm("T", "N", &m, &NN, &m, &g_one, L, &m, rprev, &m, &g_zero, r, &m);
        }
        
        if(N1prev) {
#ifdef DONT_USE_LAPACK_DGEMM_ADD
            int i;
            double *LtN2    = (double *)calloc(m*m, sizeof(double));
            double *LtNL    = LtN2;
#endif
            /*** N1 = Z'invFZ + (L'N1 + L1'N)L ***/
            if(!r1prev) dgemm("T", "N", &m, &p, &p, &g_one, Z, &p, invF, &p, &g_zero, ZtF, &m);
            dgemm("N", "N", &m, &m, &p, &g_one, ZtF, &m, Z, &p, &g_zero, N1, &m);
            dgemm("T", "N", &m, &m, &m, &g_one, L, &m, N1prev, &m, &g_zero, LtN, &m);
#ifdef DONT_USE_LAPACK_DGEMM_ADD
            dgemm("T", "N", &m, &m, &m, &g_one, L1, &m, Nprev, &m, &g_zero, LtN2, &m);
            for(i=0; i<m*m; i++) LtN[i] += LtN2[i];
            dgemm("N", "N", &m, &m, &m, &g_one, LtN, &m, L, &m, &g_zero, LtNL, &m);
            for(i=0; i<m*m; i++) N1[i] += LtNL[i];
            free(LtNL);
#else
            dgemm("T", "N", &m, &m, &m, &g_one, L1, &m, Nprev, &m, &g_one, LtN, &m);
            dgemm("N", "N", &m, &m, &m, &g_one, LtN, &m, L, &m, &g_one, N1, &m);
#endif
        }
        
        if(N2prev) {
#ifdef DONT_USE_LAPACK_DGEMM_ADD
            int i;
            double *ZtFZ    = (double *)calloc(m*m, sizeof(double));
            double *LtN2    = (double *)calloc(m*m, sizeof(double));
            double *LtNL    = LtN2;
#endif
            if(!N1prev) {
#ifdef DONT_USE_LAPACK_DGEMM_ADD
                dgemm("T", "N", &m, &m, &m, &g_one, L, &m, N1prev, &m, &g_zero, LtN, &m);
                dgemm("T", "N", &m, &m, &m, &g_one, L1, &m, Nprev, &m, &g_zero, LtN2, &m);
                for(i=0; i<m*m; i++) LtN[i] += LtN2[i];
#else
                dgemm("T", "N", &m, &m, &m, &g_one, L, &m, N1prev, &m, &g_zero, LtN, &m);
                dgemm("T", "N", &m, &m, &m, &g_one, L1, &m, Nprev, &m, &g_one, LtN, &m);
#endif
            }
            /*** N2 = Z'F2Z + (L'N1 + L1'N)L1 + (L'N2 + L1'N1)L ***/
            dgemm("N", "N", &m, &m, &m, &g_one, LtN, &m, L1, &m, &g_zero, N2, &m);
            dgemm("T", "N", &m, &p, &p, &g_one, Z, &p, F2, &p, &g_zero, ZtF, &m);
#ifdef DONT_USE_LAPACK_DGEMM_ADD
            dgemm("N", "N", &m, &m, &p, &g_one, ZtF, &m, Z, &p, &g_zero, ZtFZ, &m);
            dgemm("T", "N", &m, &m, &m, &g_one, L, &m, N2prev, &m, &g_zero, LtN, &m);
            dgemm("T", "N", &m, &m, &m, &g_one, L1, &m, N1prev, &m, &g_zero, LtN2, &m);
            for(i=0; i<m*m; i++) LtN[i] += LtN2[i];
            dgemm("N", "N", &m, &m, &m, &g_one, LtN, &m, L, &m, &g_zero, LtNL, &m);
            for(i=0; i<m*m; i++) N2[i] = ZtFZ[i] + N2[i] + LtNL[i];
            free(ZtFZ);
            free(LtNL);
#else
            dgemm("N", "N", &m, &m, &p, &g_one, ZtF, &m, Z, &p, &g_one, N2, &m);
            dgemm("T", "N", &m, &m, &m, &g_one, L, &m, N2prev, &m, &g_zero, LtN, &m);
            dgemm("T", "N", &m, &m, &m, &g_one, L1, &m, N1prev, &m, &g_one, LtN, &m);
            dgemm("N", "N", &m, &m, &m, &g_one, LtN, &m, L, &m, &g_one, N2, &m);
#endif
        }

        if(Nprev) {
            /*** N = L'NL ***/
            dgemm("T", "N", &m, &m, &m, &g_one, L, &m, Nprev, &m, &g_zero, LtN, &m);
            dgemm("N", "N", &m, &m, &m, &g_one, LtN, &m, L, &m, &g_zero, N, &m);
        }
        
        free(ZtF);
        free(LtN);
    }
    else {
        double *TtN = (double *)calloc(m*m, sizeof(double));
        
        if(r1prev) {
            /*** r1 = T'r1 ***/
            dgemm("T", "N", &m, &NN, &m, &g_one, T, &m, r1prev, &m, &g_zero, r1, &m);
        }
        
        if(N1prev) {
            /*** N1 = T'N1L ***/
            dgemm("T", "N", &m, &m, &m, &g_one, T, &m, N1prev, &m, &g_zero, TtN, &m);
            dgemm("N", "N", &m, &m, &m, &g_one, TtN, &m, L, &m, &g_zero, N1, &m);
        }
        
        if(N2prev) {
            /*** N2 = T'N2T ***/
            dgemm("T", "N", &m, &m, &m, &g_one, T, &m, N2prev, &m, &g_zero, TtN, &m);
            dgemm("N", "N", &m, &m, &m, &g_one, TtN, &m, T, &m, &g_zero, N2, &m);
        }
        
        free(TtN);
        
        backward_iter_normal(p, m, NN, Z, v, invF, L, rprev, Nprev, r, N);
    }
}

/*=============================================*\
    Missing data backward recursion iteration
\*=============================================*/
void backward_iter_miss(int m, int NN, const double *T,
        const double *rprev, const double *r1prev,
        const double *Nprev, const double *N1prev, const double *N2prev,
        double *r, double *r1, double *N, double *N1, double *N2)
{
    double *TtN = (double *)calloc(m*m, sizeof(double));
    
    if(rprev) {
        /*** r = T'r ***/
        dgemm("T", "N", &m, &NN, &m, &g_one, T, &m, rprev, &m, &g_zero, r, &m);
    }
    
    if(r1prev) {
        /*** r1 = T'r1 ***/
        dgemm("T", "N", &m, &NN, &m, &g_one, T, &m, r1prev, &m, &g_zero, r1, &m);
    }

    if(Nprev) {
        /*** N = T'NT ***/
        dgemm("T", "N", &m, &m, &m, &g_one, T, &m, Nprev, &m, &g_zero, TtN, &m);
        dgemm("N", "N", &m, &m, &m, &g_one, TtN, &m, T, &m, &g_zero, N, &m);
    }

    if(N1prev) {
        /*** N1 = T'N1T ***/
        dgemm("T", "N", &m, &m, &m, &g_one, T, &m, N1prev, &m, &g_zero, TtN, &m);
        dgemm("N", "N", &m, &m, &m, &g_one, TtN, &m, T, &m, &g_zero, N1, &m);
    }
    
    if(N2prev) {
        /*** N2 = T'N2T ***/
        dgemm("T", "N", &m, &m, &m, &g_one, T, &m, N2prev, &m, &g_zero, TtN, &m);
        dgemm("N", "N", &m, &m, &m, &g_one, TtN, &m, T, &m, &g_zero, N2, &m);
    }
    
    free(TtN);
}

/*=================*\
    Fast smoother
\*=================*/
void fastsmo(int p, int m, int rr, int n, int N, const int *pknown, int pd, const int *pmiss,
        const double *H, int Hd, const double *Z, int Zd, const double *T, int Td, const double *R, int Rd,
        const double *c, int cd, const double *a1, const double *P1, const double *P1_inf,
        int d, const int *Fns, const double *v, const double *invF, const double *F2,
        const double *K, const double *L, const double *L1, const double *QRt, const double *RQRt, int RQd,
        double *r, double *r1, double *alphahat, double *epshat, double *etahat)
{
    /*** Fast smoother information ***/
    int miss    = pmiss != 0;
    
    /*** Temporary variables ***/
    int i, t, pcur;
    double *rprev, *r1prev;
    const double *T_p0;
    double *etahat_p0, *r_p0;
    double *invFv, *Ktr, *HKt;
    
    /*** Initialize ***/
    int rout        = (r != 0);
    int r1out       = (r1 != 0);
    int alphaout    = (alphahat != 0);
    int epsout      = (epshat != 0);
    int etaout      = (etahat != 0);
    int rout_0      = rout;
    int etaout_0    = etaout;
    
    /*** Allocate buffer variables ***/
    if(rout) { r_p0 = r; rprev = r + m*N*n; r = rprev - m*N; }
    if(r1out) { r1prev = r1 + m*N*d; r1 = r1prev - m*N; }
    if(alphaout) {
        if(!r1out) {
            r1      = (double *)calloc(m*N, sizeof(double));
            r1prev  = (double *)calloc(m*N, sizeof(double));
        }
        if(m > 2*rr) {
            /*** Favor calculating alphahat from R*etahat where etahat = QRt*r ***/
            if(!etaout) {
                etaout  = true;
                etahat  = (double *)calloc(rr*N*n, sizeof(double));
            }
            if(!rout) {
                r       = (double *)calloc(m*N, sizeof(double));
                rprev   = (double *)calloc(m*N, sizeof(double));
            }
        }
        else {
            /*** Favor calculation alphahat from RQRt*r ***/
            if(!rout) {
                rout    = true;
                r_p0    = (double *)calloc(m*N*(n+1), sizeof(double));
                rprev   = r_p0 + m*N*n;
                r       = rprev - m*N;
            }
        }
        /*** r(n) = r1(n) = 0 ***/
        for(i=0; i<m*N; i++) rprev[i] = r1prev[i] = 0;
    }
    else {
        if(!rout) {
            r       = (double *)calloc(m*N, sizeof(double));
            rprev   = (double *)calloc(m*N, sizeof(double));
        }
        /*** r(n) = 0 ***/
        for(i=0; i<m*N; i++) rprev[i] = 0;
        r1prev  = (double *)0;
    }
    if(epsout) {
        if(Hd) H += p*p*(n-1);
        K       += m*p*(n-1);
        epshat  += p*N*(n-1);
        invFv    = (double *)calloc(p*N, sizeof(double));
        if(m > N) Ktr = invFv;
        else HKt = (double *)calloc(p*m, sizeof(double));
    }
    if(etaout) {
        if(RQd) QRt += rr*m*(n-1);
        etahat_p0    = etahat;
        etahat      += rr*N*(n-1);
    }

    /*** Set model matrices ***/
    if(Zd) Z += p*m*(n-1);
    if(Td) { T_p0 = T; T += m*m*(n-1); }
    
    /*** Set Kalman results ***/
    v       += p*N*(n-1);
    invF    += p*p*(n-1);
    F2      += p*p*(d-1);
    L       += m*m*(n-1);
    L1      += m*m*(d-1);
    
    /*** Backward recursion ***/
    for(t=n-1; t>=0; t--) {
        /*** etahat(t) = QRtr(t) ***/
        if(etaout) dgemm("N", "N", &rr, &N, &m, &g_one, QRt, &rr, rprev, &m, &g_zero, etahat, &rr);
        
        /*** Get data dimension ***/
        pcur    = pd ? pknown[t] : p;
        
        if(miss && pmiss[t] == 0) {
            /*** Missing data case ***/
            if(epsout)
                /*** epshat(t) = 0 ***/
                for(i=0; i<pcur*N; i++) epshat[i] = 0;
            
            if(t < d) /* t+1 <= d */
                backward_iter_miss(m, N, T, rprev, r1prev, (double *)0, (double *)0, (double *)0,
                    r, r1, (double *)0, (double *)0, (double *)0);
            else
                backward_iter_miss(m, N, T, rprev, (double *)0, (double *)0, (double *)0, (double *)0,
                    r, (double *)0, (double *)0, (double *)0, (double *)0);
        }
        else {
            if(miss && pmiss[t] < pcur) pcur = pmiss[t];
            if(t < d && Fns[t]) {
                /*** Exact diffuse case ***/
                if(epsout) {
                    /*** epshat(t) = -HK'r(t) ***/
                    if(m > N) {
                        dgemm("T", "N", &pcur, &N, &m, &g_one, K, &m, rprev, &m, &g_zero, Ktr, &pcur);
                        dgemm("N", "N", &pcur, &N, &pcur, &g_minusone, H, &pcur, Ktr, &pcur, &g_zero, epshat, &pcur);
                    }
                    else {
                        dgemm("N", "T", &pcur, &m, &pcur, &g_one, H, &pcur, K, &m, &g_zero, HKt, &pcur);
                        dgemm("N", "N", &pcur, &N, &m, &g_minusone, HKt, &pcur, rprev, &m, &g_zero, epshat, &pcur);
                    }
                }
                
                backward_iter_init(pcur, m, N, Z, T, true, v, invF, F2, L, L1, rprev, r1prev, (double *)0,
                    (double *)0, (double *)0, r, r1, (double *)0, (double *)0, (double *)0);
            }
            else {
                /*** Normal case ***/
                if(epsout) {
                    /*** epshat(t) = H(invFv - K'r(t)) ***/
                    dgemm("N", "N", &pcur, &N, &pcur, &g_one, invF, &pcur, v, &pcur, &g_zero, invFv, &pcur);
                    dgemm("T", "N", &pcur, &N, &m, &g_minusone, K, &m, rprev, &m, &g_one, invFv, &pcur);
                    dgemm("N", "N", &pcur, &N, &pcur, &g_one, H, &pcur, invFv, &pcur, &g_zero, epshat, &pcur);
                }

                if(t < d)
                    backward_iter_init(pcur, m, N, Z, T, false, v, invF, F2, L, L1, rprev, r1prev, (double *)0,
                        (double *)0, (double *)0, r, r1, (double *)0, (double *)0, (double *)0);
                else
                    backward_iter_normal(pcur, m, N, Z, v, invF, L, rprev, (double *)0, r, (double *)0);
            }
        }
        
        /*** Get next model matrices ***/
        if(Zd) Z -= p*m;
        if(Td) T -= m*m;
        
        /*** Get next Kalman results ***/
        v       -= p*N;
        invF    -= p*p;
        L       -= m*m;
        if(t < d) {
            F2      -= p*p;
            L1      -= m*m;
            if(r1out) { r1prev = r1; r1 -= m*N; }
            else if(r1prev) { double *rtemp = r1prev; r1prev = r1; r1 = rtemp; }
        }
        
        /*** Get next output ***/
        if(rout) { rprev = r; r -= m*N; }
        else { double *rtemp = rprev; rprev = r; r = rtemp; }
        if(epsout) {
            if(Hd) H -= p*p;
            K       -= m*p;
            epshat  -= p*N;
        }
        if(etaout) {
            if(RQd) QRt -= rr*m;
            etahat  -= rr*N;
        }
    }
    
    if(alphaout) {
        if(Td) T = T_p0;
        if(m > 2*rr) faststatesmo2(m, rr, n, N, T, Td, R, Rd, c, cd, a1, P1, P1_inf, rprev, r1prev, etahat_p0, alphahat);
        else faststatesmo1(m, n, N, T, Td, c, cd, a1, P1, P1_inf, RQRt, RQd, r_p0, r1prev, alphahat);
    }

    /*** Cleanup ***/
    if(alphaout) {
        if(!r1out) { free(r1); free(r1prev); }
        if(m > 2*rr) {
            if(!etaout_0) free(etahat_p0);
            if(!rout) { free(r); free(rprev); }
        }
        else if(!rout_0) free(r_p0);
    }
    else if(!rout) { free(r); free(rprev); }
    if(epsout) {
        free(invFv);
        if(m <= N) free(HKt);
    }
}

/*=======================*\
    Fast state smoother
\*=======================*/
void faststatesmo1(int m, int n, int N, const double *T, int Td, const double *c, int cd,
        const double *a1, const double *P1, const double *P1_inf, const double *RQRt, int RQd,
        const double *r, const double *r1, double *alphahat)
{
    int i, t;
    double *alphahatprev = alphahat;
    
    /*** alphahat(1) = a1 + P1r(0) + P1_infr1(0) ***/
    for(i=0; i<m*N; i++) alphahat[i] = a1[i%m];
    dgemm("N", "N", &m, &N, &m, &g_one, P1, &m, r, &m, &g_one, alphahat, &m);
    if(P1_inf) dgemm("N", "N", &m, &N, &m, &g_one, P1_inf, &m, r1, &m, &g_one, alphahat, &m);
    r           += m*N;
    alphahat    += m*N;
    
    for(t=0; t<n-1; t++) {
        /*** alphahat(t+1) = c(t) + T(t)alphahat(t) + R(t)Q(t)R(t)'r(t) ***/
        for(i=0; i<m*N; i++) alphahat[i] = c[i%m];
        dgemm("N", "N", &m, &N, &m, &g_one, T, &m, alphahatprev, &m, &g_one, alphahat, &m);
        dgemm("N", "N", &m, &N, &m, &g_one, RQRt, &m, r, &m, &g_one, alphahat, &m);
        
        alphahatprev     = alphahat;
        alphahat        += m*N;
        r               += m*N;
        if(Td) T += m*m;
        if(cd) c += m;
        if(RQd) RQRt += m*m;
    }
}

void faststatesmo2(int m, int rr, int n, int N, const double *T, int Td, const double *R, int Rd,
        const double *c, int cd, const double *a1, const double *P1, const double *P1_inf,
        const double *r, const double *r1, const double *etahat, double *alphahat)
{
    int i, t;
    double *alphahatprev = alphahat;
    
    /*** alphahat(1) = a1 + P1r(0) + P1_infr1(0) ***/
    for(i=0; i<m*N; i++) alphahat[i] = a1[i%m];
    dgemm("N", "N", &m, &N, &m, &g_one, P1, &m, r, &m, &g_one, alphahat, &m);
    if(P1_inf) dgemm("N", "N", &m, &N, &m, &g_one, P1_inf, &m, r1, &m, &g_one, alphahat, &m);
    alphahat    += m*N;
    
    for(t=0; t<n-1; t++) {
        /*** alphahat(t+1) = c(t) + T(t)alphahat(t) + R(t)etahat(t) ***/
        for(i=0; i<m*N; i++) alphahat[i] = c[i%m];
        dgemm("N", "N", &m, &N, &m, &g_one, T, &m, alphahatprev, &m, &g_one, alphahat, &m);
        dgemm("N", "N", &m, &N, &rr, &g_one, R, &m, etahat, &rr, &g_one, alphahat, &m);

        alphahatprev     = alphahat;
        alphahat        += m*N;
        etahat          += rr*N;
        if(Td) T += m*m;
        if(Rd) R += m*rr;
        if(cd) c += m;
    }
}

/*==================*\
    State smoother
\*==================*/
void statesmo(int p, int m, int n, int NN, const int *pknown, int pd,
        const int *pmiss, const double *Z, int Zd, const double *T, int Td,
        const double *a, const double *P, const double *P_inf,
        int d, const int *Fns, const double *v, const double *invF, const double *F2,
        const double *L, const double *L1, double *r, double *r1,
        double *N, double *N1, double *N2, double *alphahat, double *V)
{
    /*** State smoother information ***/
    int miss    = pmiss != 0;
    
    /*** Temporary variables ***/
    int i, j, t, pcur;
    double *rprev   = (double *)0;
    double *r1prev  = (double *)0;
    double *Nprev   = (double *)0;
    double *N1prev  = (double *)0;
    double *N2prev  = (double *)0;
    double *PN, *PNP;
    
    /*** Initialize ***/
    int rout        = (r != 0);
    int r1out       = (r1 != 0);
    int Nout        = (N != 0);
    int N1out       = (N1 != 0);
    int N2out       = (N2 != 0);
    int alphaout    = (alphahat != 0);
    int Vout        = (V != 0);
    if(rout) {
        rprev       = r + m*NN*n;
        r           = rprev - m*NN;
        /*** r(n) = 0 ***/
        for(i=0; i<m*NN; i++) rprev[i] = 0;
    }
    if(r1out) {
        r1prev      = r1 + m*NN*d;
        r1          = r1prev - m*NN;
        /*** r1(d) = 0 ***/
        for(i=0; i<m*NN; i++) r1prev[i] = 0;
    }
    if(Nout) {
        Nprev       = N + m*m*n;
        N           = Nprev - m*m;
        /*** N(n) = 0 ***/
        for(i=0; i<m*m; i++) Nprev[i] = 0;
    }
    if(N1out) {
        N1prev      = N1 + m*m*d;
        N1          = N1prev - m*m;
        /*** N1(d) = 0 ***/
        for(i=0; i<m*m; i++) N1prev[i] = 0;
    }
    if(N2out) {
        N2prev      = N2 + m*m*d;
        N2          = N2prev - m*m;
        /*** N2(d) = 0 ***/
        for(i=0; i<m*m; i++) N2prev[i] = 0;
    }
    if(alphaout) {
        for(i=0; i<m*NN*n; i++) alphahat[i] = a[i];
        alphahat    += m*NN*(n-1);
        if(!rout) {
            r       = (double *)calloc(m*NN, sizeof(double));
            rprev   = (double *)calloc(m*NN, sizeof(double));
            for(i=0; i<m*NN; i++) rprev[i] = 0;
        }
        if(!r1out) {
            r1      = (double *)calloc(m*NN, sizeof(double));
            r1prev  = (double *)calloc(m*NN, sizeof(double));
            for(i=0; i<m*NN; i++) r1prev[i] = 0;
        }
    }
    if(Vout) {
        for(i=0; i<m*m*n; i++) V[i] = P[i];
        V   += m*m*(n-1);
        if(!Nout) {
            N       = (double *)calloc(m*m, sizeof(double));
            Nprev   = (double *)calloc(m*m, sizeof(double));
            for(i=0; i<m*m; i++) Nprev[i] = 0;
        }
        if(!N1out) {
            N1      = (double *)calloc(m*m, sizeof(double));
            N1prev  = (double *)calloc(m*m, sizeof(double));
            for(i=0; i<m*m; i++) N1prev[i] = 0;
        }
        if(!N2out) {
            N2      = (double *)calloc(m*m, sizeof(double));
            N2prev  = (double *)calloc(m*m, sizeof(double));
            for(i=0; i<m*m; i++) N2prev[i] = 0;
        }
        PN  = (double *)calloc(m*m, sizeof(double));
        PNP = (double *)calloc(m*m, sizeof(double));
    }
    
    /*** Set model matrices ***/
    if(Zd) Z += p*m*(n-1);
    if(Td) T += m*m*(n-1);
    
    /*** Set Kalman results ***/
    a       += m*NN*(n-1);
    P       += m*m*(n-1);
    v       += p*NN*(n-1);
    invF    += p*p*(n-1);
    L       += m*m*(n-1);
    if(d > 0) {
        P_inf   += m*m*(d-1);
        F2      += p*p*(d-1);
        L1      += m*m*(d-1);
    }
    
    /*** Backward recursion ***/
    for(t=n-1; t>=0; t--) {
        /*** Get data dimension ***/
        pcur    = pd ? pknown[t] : p;
        
        if(miss && pmiss[t] == 0) {
            /*** Missing data case ***/
            if(t < d) /* t+1 <= d */
                backward_iter_miss(m, NN, T, rprev, r1prev, Nprev, N1prev, N2prev, r, r1, N, N1, N2);
            else
                backward_iter_miss(m, NN, T, rprev, (double *)0, Nprev, (double *)0, (double *)0,
                    r, (double *)0, N, (double *)0, (double *)0);
        }
        else {
            if(miss && pmiss[t] < pcur) pcur = pmiss[t];
            if(t < d)
                /*** Exact diffuse case ***/
                backward_iter_init(pcur, m, NN, Z, T, Fns[t], v, invF, F2, L, L1, rprev, r1prev, Nprev, N1prev, N2prev,
                    r, r1, N, N1, N2);
            else
                /*** Normal case ***/
                backward_iter_normal(pcur, m, NN, Z, v, invF, L, rprev, Nprev, r, N);
        }

        if(alphaout) {
#ifdef DONT_USE_LAPACK_DGEMM_ADD
            double *Pr  = (double *)calloc(m*NN, sizeof(double));

            /*** alphahat(t) = a(t) + P(t)r(t-1) ***/
            dgemm("N", "N", &m, &NN, &m, &g_one, P, &m, r, &m, &g_zero, Pr, &m);
            for(i=0; i<m*NN; i++) alphahat[i] = a[i] + Pr[i];
            
            if(t < d) {
                /*** alphahat(t) = ... + P_inf(t)r1(t-1) ***/
                for(i=0; i<m*NN; i++) alphahat[i] = Pr[i];
                dgemm("N", "N", &m, &NN, &m, &g_one, P_inf, &m, r1, &m, &g_zero, Pr, &m);
                for(i=0; i<m*NN; i++) alphahat[i] = a[i] + (alphahat[i] + Pr[i]);
            }
            
            free(Pr);
#else
            /*** alphahat(t) = a(t) + P(t)r(t-1) ***/
            dgemm("N", "N", &m, &NN, &m, &g_one, P, &m, r, &m, &g_one, alphahat, &m);
            
            if(t < d)
                /*** alphahat(t) = ... + P_inf(t)r1(t-1) ***/
                dgemm("N", "N", &m, &NN, &m, &g_one, P_inf, &m, r1, &m, &g_one, alphahat, &m);
#endif
        }

        if(Vout) {
            /*** V(t) = P(t) - P(t)N(t-1)P(t) ***/
            dgemm("N", "N", &m, &m, &m, &g_one, P, &m, N, &m, &g_zero, PN, &m);
            dgemm("N", "N", &m, &m, &m, &g_minusone, PN, &m, P, &m, &g_one, V, &m);
            
            if(t < d) {
                /*** V(t) = ... - (P_inf(t)N1(t-1)P(t))' - P_inf(t)N1(t-1)P(t) - P_inf(t)N2(t-1)P_inf(t) ***/
                dgemm("N", "N", &m, &m, &m, &g_one, P_inf, &m, N1, &m, &g_zero, PN, &m);
                dgemm("N", "N", &m, &m, &m, &g_one, PN, &m, P, &m, &g_zero, PNP, &m);
                for(i=0; i<m; i++) for(j=0; j<m; j++) V[i*m+j] -= PNP[i*m+j] + PNP[j*m+i];
                dgemm("N", "N", &m, &m, &m, &g_one, P_inf, &m, N2, &m, &g_zero, PN, &m);
                dgemm("N", "N", &m, &m, &m, &g_minusone, PN, &m, P_inf, &m, &g_one, V, &m);
            }
        }

        /*** Get next model matrices ***/
        if(Zd) Z -= p*m;
        if(Td) T -= m*m;
        
        /*** Get next Kalman results ***/
        a       -= m*NN;
        P       -= m*m;
        v       -= p*NN;
        invF    -= p*p;
        L       -= m*m;
        if(t < d) {
            P_inf   -= m*m;
            F2      -= p*p;
            L1      -= m*m;
            if(r1out) { r1prev = r1; r1 -= m*NN; }
            else if(r1prev) { double *rtemp = r1prev; r1prev = r1; r1 = rtemp; }
            if(N1out) { N1prev = N1; N1 -= m*m; }
            else if(N1prev) { double *Ntemp = N1prev; N1prev = N1; N1 = Ntemp; }
            if(N2out) { N2prev = N2; N2 -= m*m; }
            else if(N2prev) { double *Ntemp = N2prev; N2prev = N2; N2 = Ntemp; }
        }
        
        /*** Get next output ***/
        if(rout) { rprev = r; r -= m*NN; }
        else if(rprev) { double *rtemp = rprev; rprev = r; r = rtemp; }
        if(Nout) { Nprev = N; N -= m*m; }
        else if(Nprev) { double *Ntemp = Nprev; Nprev = N; N = Ntemp; }
        if(alphaout) alphahat -= m*NN;
        if(Vout) V -= m*m;
    }
    
    /*** Cleanup ***/
    if(alphaout) {
        if(!rout) { free(r); free(rprev); }
        if(!r1out) { free(r1); free(r1prev); }
    }
    if(Vout) {
        if(!Nout) { free(N); free(Nprev); }
        if(!N1out) { free(N1); free(N1prev); }
        if(!N2out) { free(N2); free(N2prev); }
        free(PN);
        free(PNP);
    }
}

/*========================*\
    Disturbance smoother
\*========================*/
void disturbsmo(int p, int m, int rr, int n, int NN, const int *pknown, int pd, const int *pmiss,
        const double *H, int Hd, const double *Z, int Zd, const double *T, int Td, const double *Q, int Qd,
        int d, const int *Fns, const double *v, const double *invF, const double *K,
        const double *L, const double *RQ, const double *QRt, int RQd,
        double *r, double *N, double *epshat, double *etahat, double *epsvarhat, double *etavarhat)
{
    /*** Disturbance smoother information ***/
    int miss    = pmiss != 0;
    
    /*** Temporary variables ***/
    int i, t, pcur;
    double *rprev   = (double *)0;
    double *Nprev   = (double *)0;
    double *QRtN, *invFv, *Ktr, *HKt, *KtN, *KtNK, *FKtNK, *HKtNK, *HFKtNK;
    
    /*** Initialize ***/
    int rout        = (r != 0);
    int Nout        = (N != 0);
    int epsout      = (epshat != 0);
    int etaout      = (etahat != 0);
    int epsvarout   = (epsvarhat != 0);
    int etavarout   = (etavarhat != 0);
    if(rout) {
        rprev       = r + m*NN*n;
        r           = rprev - m*NN;
        /*** r(n) = 0 ***/
        for(i=0; i<m*NN; i++) rprev[i] = 0;
    }
    if(Nout) {
        Nprev       = N + m*m*n;
        N           = Nprev - m*m;
        /*** N(n) = 0 ***/
        for(i=0; i<m*m; i++) Nprev[i] = 0;
    }
    if(epsout || etaout) {
        if(epsout) {
            epshat += p*NN*(n-1);
            invFv    = (double *)calloc(p*NN, sizeof(double));
            if(m > NN) Ktr = invFv;
            else HKt = (double *)calloc(p*m, sizeof(double));
        }
        if(etaout) etahat += rr*NN*(n-1);
        if(!rout) {
            r       = (double *)calloc(m*NN, sizeof(double));
            rprev   = (double *)calloc(m*NN, sizeof(double));
            for(i=0; i<m*NN; i++) rprev[i] = 0;
        }
    }
    if(epsvarout || etavarout) {
        if(epsvarout) {
            if(Hd) for(i=0; i<p*p*n; i++) epsvarhat[i] = H[i]; /* May waste computation due to changing p. */
            else for(i=0; i<p*p*n; i++) epsvarhat[i] = H[i%(p*p)];
            epsvarhat += p*p*(n-1);
            KtN     = (double *)calloc(p*m, sizeof(double));
            KtNK    = FKtNK = (double *)calloc(p*p, sizeof(double));
            HKtNK   = HFKtNK = (double *)calloc(p*p, sizeof(double));
        }
        if(etavarout) {
            if(Qd) for(i=0; i<rr*rr*n; i++) etavarhat[i] = Q[i];
            else for(i=0; i<rr*rr*n; i++) etavarhat[i] = Q[i%(rr*rr)];
            etavarhat += rr*rr*(n-1);
            QRtN    = (double *)calloc(rr*m, sizeof(double));
        }
        if(!Nout) {
            N       = (double *)calloc(m*m, sizeof(double));
            Nprev   = (double *)calloc(m*m, sizeof(double));
            for(i=0; i<m*m; i++) Nprev[i] = 0;
        }
    }
    
    /*** Set model matrices ***/
    if(Hd) H += p*p*(n-1);
    if(Zd) Z += p*m*(n-1);
    if(Td) T += m*m*(n-1);
    
    /*** Set Kalman results ***/
    v       += p*NN*(n-1);
    invF    += p*p*(n-1);
    K       += m*p*(n-1);
    L       += m*m*(n-1);
    if(RQd) {
        RQ  += m*rr*(n-1);
        QRt += rr*m*(n-1);
    }
    
    /*** Backward recursion ***/
    for(t=n-1; t>=0; t--) {
        /*** etahat(t) = Q(t)R(t)'r(t) ***/
        if(etaout) dgemm("N", "N", &rr, &NN, &m, &g_one, QRt, &rr, rprev, &m, &g_zero, etahat, &rr);
        
        /*** etavarhat(t) = Q(t) - Q(t)R(t)'N(t)R(t)Q(t) ***/
        if(etavarout) {
            dgemm("N", "N", &rr, &m, &m, &g_one, QRt, &rr, Nprev, &m, &g_zero, QRtN, &rr);
            dgemm("N", "N", &rr, &rr, &m, &g_minusone, QRtN, &rr, RQ, &m, &g_one, etavarhat, &rr);
        }
        
        /*** Get data dimension ***/
        pcur    = pd ? pknown[t] : p;
        
        if(miss && pmiss[t] == 0) {
            /*** Missing data case ***/
            if(epsout)
                /*** epshat(t) = 0 ***/
                for(i=0; i<pcur*NN; i++) epshat[i] = 0;
            
            /******* TODO: Confirm epsvarhat(t) = H(t) when data missing. *******/

            backward_iter_miss(m, NN, T, rprev, (double *)0, Nprev, (double *)0, (double *)0,
                r, (double *)0, N, (double *)0, (double *)0);
        }
        else {
            if(miss && pmiss[t] < pcur) pcur = pmiss[t];
            if(t < d && Fns[t]) {
                /*** Exact diffuse case ***/
                if(epsout) {
                    /*** epshat(t) = -H(t)K(t)'r(t) ***/
                    if(m > NN) {
                        dgemm("T", "N", &pcur, &NN, &m, &g_one, K, &m, rprev, &m, &g_zero, Ktr, &pcur);
                        dgemm("N", "N", &pcur, &NN, &pcur, &g_minusone, H, &pcur, Ktr, &pcur, &g_zero, epshat, &pcur);
                    }
                    else {
                        dgemm("N", "T", &pcur, &m, &pcur, &g_one, H, &pcur, K, &m, &g_zero, HKt, &pcur);
                        dgemm("N", "N", &pcur, &NN, &m, &g_minusone, HKt, &pcur, rprev, &m, &g_zero, epshat, &pcur);
                    }
                }
                
                if(epsvarout) {
                    /*** epsvarhat(t) = H(t) - H(t)K(t)'N(t)K(t)H(t) ***/
                    dgemm("T", "N", &pcur, &m, &m, &g_one, K, &m, Nprev, &m, &g_zero, KtN, &pcur);
                    dgemm("N", "N", &pcur, &pcur, &m, &g_one, KtN, &pcur, K, &m, &g_zero, KtNK, &pcur);
                    dgemm("N", "N", &pcur, &pcur, &pcur, &g_one, H, &pcur, KtNK, &pcur, &g_zero, HKtNK, &pcur);
                    dgemm("N", "N", &pcur, &pcur, &pcur, &g_minusone, HKtNK, &pcur, H, &pcur, &g_one, epsvarhat, &pcur);
                }

                backward_iter_init(pcur, m, NN, Z, T, true, v, invF, (double *)0, L, (double *)0, rprev, (double *)0, Nprev,
                    (double *)0, (double *)0, r, (double *)0, N, (double *)0, (double *)0);
            }
            else {
                /*** Normal case ***/
                if(epsout) {
                    /*** epshat(t) = H(invFv - K'r(t)) ***/
                    dgemm("N", "N", &pcur, &NN, &pcur, &g_one, invF, &pcur, v, &pcur, &g_zero, invFv, &pcur);
                    dgemm("T", "N", &pcur, &NN, &m, &g_minusone, K, &m, rprev, &m, &g_one, invFv, &pcur);
                    dgemm("N", "N", &pcur, &NN, &pcur, &g_one, H, &pcur, invFv, &pcur, &g_zero, epshat, &pcur);
                }
                
                if(epsvarout) {
                    /*** epsvarhat(t) = H(t) - H(t)(invF(t) + K(t)'N(t)K(t))H(t) ***/
                    dgemm("T", "N", &pcur, &m, &m, &g_one, K, &m, Nprev, &m, &g_zero, KtN, &pcur);
                    for(i=0; i<pcur*pcur; i++) FKtNK[i] = invF[i];
                    dgemm("N", "N", &pcur, &pcur, &m, &g_one, KtN, &pcur, K, &m, &g_one, FKtNK, &pcur);
                    dgemm("N", "N", &pcur, &pcur, &pcur, &g_one, H, &pcur, FKtNK, &pcur, &g_zero, HFKtNK, &pcur);
                    dgemm("N", "N", &pcur, &pcur, &pcur, &g_minusone, HFKtNK, &pcur, H, &pcur, &g_one, epsvarhat, &pcur);
                }
                
                backward_iter_normal(pcur, m, NN, Z, v, invF, L, rprev, Nprev, r, N);
            }
        }
        
        /*** Get next model matrices ***/
        if(Hd) H -= p*p;
        if(Zd) Z -= p*m;
        if(Td) T -= m*m;
        
        /*** Get next Kalman results ***/
        v       -= p*NN;
        invF    -= p*p;
        K       -= m*p;
        L       -= m*m;
        if(RQd) {
            RQ  -= m*rr;
            QRt -= rr*m;
        }
        
        /*** Get next output ***/
        if(rout) { rprev = r; r -= m*NN; }
        else if(rprev) { double *rtemp = rprev; rprev = r; r = rtemp; }
        if(Nout) { Nprev = N; N -= m*m; }
        else if(Nprev) { double *Ntemp = Nprev; Nprev = N; N = Ntemp; }
        if(epsout) epshat -= p*NN;
        if(etaout) etahat -= rr*NN;
        if(epsvarout) epsvarhat -= p*p;
        if(etavarout) etavarhat -= rr*rr;
    }
    
    /*** Cleanup ***/
    if(epsout || etaout) {
        if(epsout) { free(invFv); if(m <= NN) free(HKt); }
        if(!rout) { free(r); free(rprev); }
    }
    if(epsvarout || etavarout) {
        if(epsvarout) { free(KtN); free(KtNK); free(HKtNK); }
        if(etavarout) free(QRtN);
        if(!Nout) { free(N); free(Nprev); }
    }
}

/*======================*\
    Smoothing cumulant
\*======================*/
void smocum(int p, int m, int n, const int *pknown, int pd, const int *pmiss,
        const double *Z, int Zd, const double *T, int Td, int d, const int *Fns,
        const double *v, const double *invF, const double *K, const double *L,
        double *uuD, double *rrN)
{
    /*** Disturbance smoother information ***/
    int miss    = pmiss != 0;
    
    /*** Temporary variables ***/
    int i, j, t, pcur;
    double *r, *rprev, *rtemp;
    double *N, *Nprev, *Ntemp;
    double *u, *D, *KtN;
    
    /*** Initialize ***/
    r       = (double *)calloc(m, sizeof(double));
    rprev   = (double *)calloc(m, sizeof(double));
    N       = (double *)calloc(m*m, sizeof(double));
    Nprev   = (double *)calloc(m*m, sizeof(double));
    u       = (double *)calloc(p, sizeof(double));
    D       = (double *)calloc(p*p, sizeof(double));
    KtN     = (double *)calloc(p*m, sizeof(double));
    for(i=0; i<m; i++) rprev[i] = 0;
    for(i=0; i<m*m; i++) Nprev[i] = 0;
    
    /*** Set model matrices ***/
    if(Zd) Z += p*m*(n-1);
    if(Td) T += m*m*(n-1);
    
    /*** Set Kalman results ***/
    v       += p*(n-1);
    invF    += p*p*(n-1);
    K       += m*p*(n-1);
    L       += m*m*(n-1);
    
    /*** Set output ***/
    uuD    += p*p*(n-1);
    rrN    += m*m*(n-1);
    
    /*** Backward recursion ***/
    for(t=n-1; t>=0; t--) {
        /*** Compute rrN(t) = r(t)*r(t)' - N(t) ***/
        for(i=0; i<m; i++) for(j=0; j<m; j++) rrN[i*m+j] = rprev[i]*rprev[j] - Nprev[i*m+j];
        
        /*** Get data dimension ***/
        pcur    = pd ? pknown[t] : p;
        
        if(miss && pmiss[t] == 0) {
            /*** Missing data case ***/
            /*** uuD(t) = 0 ***/
            for(i=0; i<pcur*pcur; i++) uuD[i] = 0;

            backward_iter_miss(m, 1, T, rprev, (double *)0, Nprev, (double *)0, (double *)0,
                r, (double *)0, N, (double *)0, (double *)0);
        }
        else {
            if(miss && pmiss[t] < pcur) pcur = pmiss[t];
            if(t < d && Fns[t]) {
                /*** Exact diffuse case ***/
                /*** u(t) = -K(t)'r(t) ***/
                i   = 1;
                dgemv("T", &m, &pcur, &g_minusone, K, &m, rprev, &i, &g_zero, u, &i);
                
                /*** D(t) = K(t)'N(t)K(t) ***/
                dgemm("T", "N", &pcur, &m, &m, &g_one, K, &m, Nprev, &m, &g_zero, KtN, &pcur);
                dgemm("N", "N", &pcur, &pcur, &m, &g_one, KtN, &pcur, K, &m, &g_zero, D, &pcur);

                /*** Compute uuD(t) = u(t)*u(t)' - D(t) ***/
                for(i=0; i<pcur; i++) for(j=0; j<pcur; j++) uuD[i*pcur+j] = u[i]*u[j] - D[i*pcur+j];

                backward_iter_init(pcur, m, 1, Z, T, true, v, invF, (double *)0, L, (double *)0, rprev, (double *)0, Nprev,
                    (double *)0, (double *)0, r, (double *)0, N, (double *)0, (double *)0);
            }
            else {
                /*** Normal case ***/
                /*** u(t) = invFv - K'r(t) ***/
                i   = 1;
                dgemv("N", &pcur, &pcur, &g_one, invF, &pcur, v, &i, &g_zero, u, &i);
                dgemv("T", &m, &pcur, &g_minusone, K, &m, rprev, &i, &g_one, u, &i);

                /*** D(t) = invF(t) + K(t)'N(t)K(t) ***/
                for(i=0; i<pcur*pcur; i++) D[i] = invF[i];
                dgemm("T", "N", &pcur, &m, &m, &g_one, K, &m, Nprev, &m, &g_zero, KtN, &pcur);
                dgemm("N", "N", &pcur, &pcur, &m, &g_one, KtN, &pcur, K, &m, &g_one, D, &pcur);
                
                /*** Compute uuD(t) = u(t)*u(t)' - D(t) ***/
                for(i=0; i<pcur; i++) for(j=0; j<pcur; j++) uuD[i*pcur+j] = u[i]*u[j] - D[i*pcur+j];

                backward_iter_normal(pcur, m, 1, Z, v, invF, L, rprev, Nprev, r, N);
            }
        }
        
        /*** Get next model matrices ***/
        if(Zd) Z -= p*m;
        if(Td) T -= m*m;
        
        /*** Get next Kalman results ***/
        v       -= p;
        invF    -= p*p;
        K       -= m*p;
        L       -= m*m;
        
        /*** Get next output ***/
        rtemp   = rprev; rprev = r; r = rtemp;
        Ntemp   = Nprev; Nprev = N; N = Ntemp;
        uuD     -= p*p;
        rrN     -= m*m;
    }
    
    /*** Cleanup ***/
    free(r);
    free(rprev);
    free(N);
    free(Nprev);
    free(u);
    free(D);
    free(KtN);
}
