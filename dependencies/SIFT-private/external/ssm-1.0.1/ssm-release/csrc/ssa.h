/**************************************\
    Filename: ssa.h
    Author: Jyh-Ying Peng (´^´¼·­)
    Year: 2006-2007
    Version: 1.4
\**************************************/

#ifndef SSA_H_INCLUDED_20070119
#define SSA_H_INCLUDED_20070119

#include <stdlib.h>
#include <math.h>
#include "mt19937ar.h"

/*=======================================================================*\
                          Kalman filter iteration
  -----------------------------------------------------------------------
    Input:
        p           - data dimension.
        m           - state dimension.
        N           - number of data sets.
        y           - p*N matrix, observation data.
        H           - p*p matrix, observation disturbance variance.
        Z           - p*m matrix, state to observation transform.
        T           - m*m matrix, state update transform.
        c           - m*1 vector, state update constant.
        RQRt        - m*m matrix, state disturbance variance.
        aprev       - m*N matrix, current filtered state vector.
        Pprev       - m*m matrix, current filtered state variance.
        P_infprev   - m*m matrix, current diffuse state variance.
    Output:
        Fns     - true if F_inf = Z*P_infprev*Z' is non-singular.
        v       - p*N matrix, innovation.
        detF    - determinant of F, can be NULL.
        F       - p*p matrix, innovation variance, can be NULL.
        invF    - p*p matrix, inverse of F.
        F2      - p*p matrix, diffuse component of F.
        K       - m*p matrix.
        K1      - m*p matrix, diffuse component of K.
        L       - m*m matrix.
        L1      - m*m matrix, diffuse component of L.
        a       - m*N matrix, next filtered state vector.
        P       - m*m matrix, next filtered state variance.
        P_inf   - m*m matrix, next diffuse state variance.
    Options:
        tol0        - tolerance from 0, any number less than tol0 will be
                      treated as 0.
        inv_method  - 0 is general; 1 is symmetric; 2 is positive.
\*=======================================================================*/
void kalman_iter_normal(int p, int m, const double *H, const double *Z, const double *T,
        const double *RQRt, const double *Pprev,
        double *detF, double *F, double *invF, double *K, double *L, double *P,
        int inv_method);
void kalman_iter_init(int p, int m, const double *H, const double *Z, const double *T,
        const double *RQRt, const double *Pprev, const double *P_infprev,
        int *Fns, double *detF, double *F, double *invF, double *F2,
        double *K, double *K1, double *L, double *L1, double *P, double *P_inf,
        double tol0, int inv_method);
void kalman_iter_miss(int p, int m, const double *H, const double *Z, const double *T,
        const double *RQRt, const double *Pprev, const double *P_infprev,
        double *detF, double *F, double *invF, double *K, double *L, double *P, double *P_inf,
        int inv_method);
void kalman_iter_data(int p, int m, int N, const double *y, const double *Z, const double *T,
        const double *c, const double *K, const double *aprev, double *v, double *a);

/*=======================================================================*\
                              Kalman filter
  -----------------------------------------------------------------------
    Input:
        p           - data dimension.
        m           - state dimension.
        r           - state disturbance dimension.
        n           - data size/length.
        N           - number of data sets.
        y           - p*N*n matrix, observation data.
        pknown      - 1*n vector, dynamic data dimensions, can be NULL if
                      pd = false.
        pd          - true if data dimensions are dynamic.
        pmiss       - 1*n vector, data dimension after missing data is
                      removed, can be NULL if there's no missing data.
        sta         - true if H, Z, T, R and Q are all stationary.
        H           - p*p or p*p*n matrix, observation disturbance variance.
        Hd          - true if H is dynamic (p*p*n).
        Z           - p*m or p*m*n matrix, state to observation transform.
        Zd          - true if Z is dynamic (p*m*n).
        T           - m*m or m*m*n matrix, state update transform.
        Td          - true if T is dynamic (m*m*n).
        R           - m*r or m*r*n matrix, state disturbance transform.
        Rd          - true if R is dynamic (m*r*n).
        Q           - r*r or r*r*n matrix, state disturbance variance.
        Qd          - true if Q is dynamic (r*r*n).
        c           - m*1 or m*n matrix, state update constant.
        cd          - true if c is dynamic (m*n).
        a1          - m*1 vector, initial state vector.
        P1          - m*m matrix, initial state variance.
        P1_inf      - m*m matrix, initial diffuse state variance.
    Output: (all output pointers can be NULL)
        a       - m*N*n matrix, filtered state vector sequence.
        P       - m*m*n matrix, filtered state variance sequence.
        P_inf   - m*m*min(d, n) matrix, diffuse state variance sequence.
        d       - length of non-zero diffuse variances, 0 if exact diffuse
                  initialization is not used, n+1 if the exact diffuse
                  initialization failed to converge.
        Fns     - 1*min(d, n) vector, Fns[t] is true if F_inf(t)
                  = Z(t)*P_inf(t)*Z(t)' is non-singular.
        v       - p*N*n matrix, innovation sequence.
        F       - p*p*n matrix, innovation variance sequence.
        invF    - p*p*n matrix, inverse of F.
        F2      - p*p*min(d, n) matrix, diffuse component of F.
        K       - m*p*n matrix.
        K1      - m*p*min(d, n) matrix, diffuse component of K.
        L       - m*m*n matrix.
        L1      - m*m*min(d, n) matrix, diffuse component of L.
        snlogL  - 1*N vector, semi-negative loglikelihood.
        fvar    - 1*N vector, forecast variance.
        RQ      - m*r or m*r*n matrix, R*Q.
        QRt     - r*m or r*m*n matrix, Q*R'.
        RQRt    - m*m or m*m*n matrix, R*Q*R'.
    Options:
        tol0        - tolerance from 0.
        tolP        - tolerance for Kalman filter convergence.
        inv_method  - 0 is general; 1 is symmetric; 2 is positive.
\*=======================================================================*/
void kalman(int p, int m, int r, int n, int N, const double *y,
        const int *pknown, int pd, const int *pmiss, int sta,
        const double *H, int Hd, const double *Z, int Zd, const double *T, int Td,
        const double *R, int Rd, const double *Q, int Qd, const double *c, int cd,
        const double *a1, const double *P1, const double *P1_inf,
        double *a, double *P, double *P_inf, int *d, int *Fns, double *v,
        double *F, double *invF, double *F2, double *K, double *K1, double *L, double *L1,
        double *snlogL, double *fvar, double *RQ, double *QRt, double *RQRt,
        double tol0, double tolP, int inv_method);

/*=======================================================================*\
                      Filtering and smoothing weights
  -----------------------------------------------------------------------
    Input:
        p           - data dimension.
        m           - state dimension.
        r           - state disturbance dimension.
        n           - time duration.
        t0          - time point of the weighting (1 <= t0 <= n).
        pknown      - 1*n vector, dynamic data dimensions, can be NULL if
                      pd = false.
        pd          - true if data dimensions are dynamic.
        sta         - true if H, Z, T, R and Q are all stationary.
        H           - p*p or p*p*n matrix, observation disturbance variance.
        Hd          - true if H is dynamic (p*p*n).
        Z           - p*m or p*m*n matrix, state to observation transform.
        Zd          - true if Z is dynamic (p*m*n).
        T           - m*m or m*m*n matrix, state update transform.
        Td          - true if T is dynamic (m*m*n).
        R           - m*r or m*r*n matrix, state disturbance transform.
        Rd          - true if R is dynamic (m*r*n).
        Q           - r*r or r*r*n matrix, state disturbance variance.
        Qd          - true if Q is dynamic (r*r*n).
        P1          - m*m matrix, initial state variance.
        P1_inf      - m*m matrix, initial diffuse state variance.
    Output:
        wt_a        - m*p*n matrix, filter weights a(t0) = sum_t wt_a(t)*y(t).
        wt_alpha    - m*p*n matrix, smoother weights alphahat(t0) = sum_t
                      wt_alpha(t)*y(t), can be NULL.
    Options:
        tol0        - tolerance from 0.
        tolP        - tolerance for Kalman filter convergence.
        inv_method  - 0 is general; 1 is symmetric; 2 is positive.
\*=======================================================================*/
void weights(int p, int m, int r, int n, int t0, const int *pknown, int pd,
        int sta, const double *H, int Hd, const double *Z, int Zd, const double *T, int Td,
        const double *R, int Rd, const double *Q, int Qd, const double *P1, const double *P1_inf,
        double *wt_a, double *wt_alpha, double tol0, double tolP, int inv_method);

/*=======================================================================*\
                       Backward recursion iteration
  -----------------------------------------------------------------------
    Input:
        p       - data dimension.
        m       - state dimension.
        NN      - number of data sets.
        Z       - p*m matrix, state to observation transform.
        T       - m*m matrix, state update transform.
        Fns     - true if F_inf = Z*P_infprev*Z' is non-singular.
        v       - p*NN matrix, innovation.
        invF    - p*p matrix, inverse of F.
        F2      - p*p matrix, diffuse component of F.
        L       - m*m matrix.
        L1      - m*m matrix, diffuse component of L.
        (the following inputs can be NULL)
        rprev   - m*NN matrix, current smoothing cumulant.
        r1prev  - m*NN matrix, current diffuse smoothing cumulant.
        Nprev   - m*m matrix, current smoothing cumulant variance.
        N1prev  - m*m matrix, current diffuse smoothing cumulant variance.
        N2prev  - m*m matrix, current diffuse smoothing cumulant variance.
    Output:
        (output pointers can be NULL if their correponding input is NULL)
        r   - m*NN matrix, next smoothing cumulant.
        r1  - m*NN matrix, next diffuse smoothing cumulant.
        N   - m*m matrix, next smoothing cumulant variance.
        N1  - m*m matrix, next diffuse smoothing cumulant variance.
        N2  - m*m matrix, next diffuse smoothing cumulant variance.
\*=======================================================================*/
void backward_iter_normal(int p, int m, int NN, const double *Z,
        const double *v, const double *invF, const double *L,
        const double *rprev, const double *Nprev, double *r, double *N);
void backward_iter_init(int p, int m, int NN, const double *Z, const double *T,
        int Fns, const double *v, const double *invF, const double *F2,
        const double *L, const double *L1, const double *rprev, const double *r1prev,
        const double *Nprev, const double *N1prev, const double *N2prev,
        double *r, double *r1, double *N, double *N1, double *N2);
void backward_iter_miss(int m, int NN, const double *T,
        const double *rprev, const double *r1prev,
        const double *Nprev, const double *N1prev, const double *N2prev,
        double *r, double *r1, double *N, double *N1, double *N2);

/*=======================================================================*\
                              Fast smoother
  -----------------------------------------------------------------------
    Input:
        p       - data dimension.
        m       - state dimension.
        rr      - state disturbance dimension.
        n       - data size/length.
        N       - number of data sets.
        pknown  - 1*n vector, dynamic data dimensions, can be NULL if
                  pd = false.
        pd      - true if data dimensions are dynamic.
        pmiss   - 1*n vector, data dimension after missing data is
                  removed, can be NULL if there's no missing data.
        H       - p*p or p*p*n matrix, observation disturbance variance.
        Hd      - true if H is dynamic (p*p*n).
        Z       - p*m or p*m*n matrix, state to observation transform.
        Zd      - true if Z is dynamic (p*m*n).
        T       - m*m or m*m*n matrix, state update transform.
        Td      - true if T is dynamic (m*m*n).
        R       - m*rr or m*rr*n matrix, state disturbance transform.
        Rd      - true if R is dynamic (m*rr*n).
        c       - m*1 or m*n matrix, state update constant.
        cd      - true if c is dynamic (m*n).
        a1      - m*1 vector, initial state vector.
        P1      - m*m matrix, initial state variance.
        P1_inf  - m*m matrix, initial diffuse state variance.
        d       - length of non-zero diffuse variances.
        Fns     - 1*d vector, Fns[t] is true if F_inf(t)
                  = Z(t)*P_inf(t)*Z(t)' is non-singular.
        v       - p*N*n matrix, innovation sequence.
        invF    - p*p*n matrix, inverse of F.
        F2      - p*p*d matrix, diffuse component of F.
        K       - m*p*n matrix.
        L       - m*m*n matrix.
        L1      - m*m*d matrix, diffuse component of L.
        QRt     - rr*m or rr*m*n matrix, Q*R'.
        RQRt    - m*m or m*m*n matrix, R*Q*R'.
        RQd     - true if one of R or Q is dynamic.
    Output:
        r           - m*N*(n+1) matrix, smoothing cumulant sequence.
        r1          - m*N*(d+1) matrix, diffuse smoothing cumulant sequence.
        alphahat    - m*N*n matrix, smoothed state sequence.
        epshat      - p*N*n matrix, smoothed observation disturbance sequence.
        etahat      - rr*N*n matrix, smoothed state disturbance sequence.
\*=======================================================================*/
void fastsmo(int p, int m, int rr, int n, int N, const int *pknown, int pd,
        const int *pmiss, const double *H, int Hd, const double *Z, int Zd,
        const double *T, int Td, const double *R, int Rd, const double *c, int cd,
        const double *a1, const double *P1, const double *P1_inf,
        int d, const int *Fns, const double *v, const double *invF, const double *F2,
        const double *K, const double *L, const double *L1,
        const double *QRt, const double *RQRt, int RQd,
        double *r, double *r1, double *alphahat, double *epshat, double *etahat);

/*=======================================================================*\
                           Fast state smoother
  -----------------------------------------------------------------------
    Input:
        m       - state dimension.
        rr      - state disturbance dimension.
        n       - data size/length.
        N       - number of data sets.
        T       - m*m or m*m*n matrix, state update transform.
        Td      - true if T is dynamic (m*m*n).
        R       - m*rr or m*rr*n matrix, state disturbance transform.
        Rd      - true if R is dynamic (m*rr*n).
        c       - m*1 or m*n matrix, state update constant.
        cd      - true if c is dynamic (m*n).
        a1      - m*1 vector, initial state vector.
        P1      - m*m matrix, initial state variance.
        P1_inf  - m*m matrix, initial diffuse state variance.
        RQRt    - m*m or m*m*n matrix, R*Q*R'.
        RQd     - true if one of R or Q is dynamic.
        r       - m*N or m*N*(n+1) matrix, smoothing cumulant, only needed
                  at time 0 if etahat is provided.
        r1      - m*N matrix, diffuse smoothing cumulant at time 0.
        etahat  - rr*N*n matrix, smoothed state disturbance sequence.
    Output:
        alphahat    - m*N*n matrix, smoothed state sequence.
\*=======================================================================*/
void faststatesmo1(int m, int n, int N, const double *T, int Td, const double *c, int cd,
        const double *a1, const double *P1, const double *P1_inf, const double *RQRt, int RQd,
        const double *r, const double *r1, double *alphahat);
void faststatesmo2(int m, int rr, int n, int N, const double *T, int Td, const double *R, int Rd,
        const double *c, int cd, const double *a1, const double *P1, const double *P1_inf,
        const double *r, const double *r1, const double *etahat, double *alphahat);

/*=======================================================================*\
                              State smoother
  -----------------------------------------------------------------------
    Input:
        p       - data dimension.
        m       - state dimension.
        n       - data size/length.
        NN      - number of data sets.
        pknown  - 1*n vector, dynamic data dimensions, can be NULL if
                  pd = false.
        pd      - true if data dimensions are dynamic.
        pmiss   - 1*n vector, data dimension after missing data is
                  removed, can be NULL if there's no missing data.
        Z       - p*m or p*m*n matrix, state to observation transform.
        Zd      - true if Z is dynamic (p*m*n).
        T       - m*m or m*m*n matrix, state update transform.
        Td      - true if T is dynamic (m*m*n).
        a       - m*N*n matrix, filtered state vector sequence.
        P       - m*m*n matrix, filtered state variance sequence.
        P_inf   - m*m*d matrix, diffuse state variance sequence.
        d       - length of non-zero diffuse variances.
        Fns     - 1*d vector, Fns[t] is true if F_inf(t)
                  = Z(t)*P_inf(t)*Z(t)' is non-singular.
        v       - p*N*n matrix, innovation sequence.
        invF    - p*p*n matrix, inverse of F.
        F2      - p*p*d matrix, diffuse component of F.
        L       - m*m*n matrix.
        L1      - m*m*d matrix, diffuse component of L.
    Output:
        r           - m*N*(n+1) matrix, smoothing cumulant sequence.
        r1          - m*N*(d+1) matrix, diffuse smoothing cumulant sequence.
        N           - m*m*(n+1) matrix, smoothing cumulant variance sequence.
        N1          - m*m*(d+1) matrix, diffuse smoothing cumulant variance
                      sequence.
        N2          - m*m*(d+1) matrix, diffuse smoothing cumulant variance
                      sequence.
        alphahat    - m*N*n matrix, smoothed state sequence.
        V           - m*m*n matrix, smoothed state variance sequence.
\*=======================================================================*/
void statesmo(int p, int m, int n, int NN, const int *pknown, int pd,
        const int *pmiss, const double *Z, int Zd, const double *T, int Td,
        const double *a, const double *P, const double *P_inf,
        int d, const int *Fns, const double *v, const double *invF, const double *F2,
        const double *L, const double *L1, double *r, double *r1,
        double *N, double *N1, double *N2, double *alphahat, double *V);

/*=======================================================================*\
                           Disturbance smoother
  -----------------------------------------------------------------------
    Input:
        p       - data dimension.
        m       - state dimension.
        rr      - state disturbance dimension.
        n       - data size/length.
        NN      - number of data sets.
        pknown  - 1*n vector, dynamic data dimensions, can be NULL if
                  pd = false.
        pd      - true if data dimensions are dynamic.
        pmiss   - 1*n vector, data dimension after missing data is
                  removed, can be NULL if there's no missing data.
        H       - p*p or p*p*n matrix, observation disturbance variance.
        Hd      - true if H is dynamic (p*p*n).
        Z       - p*m or p*m*n matrix, state to observation transform.
        Zd      - true if Z is dynamic (p*m*n).
        T       - m*m or m*m*n matrix, state update transform.
        Td      - true if T is dynamic (m*m*n).
        Q       - r*r or r*r*n matrix, state disturbance variance.
        Qd      - true if Q is dynamic (r*r*n).
        d       - length of non-zero diffuse variances.
        Fns     - 1*d vector, Fns[t] is true if F_inf(t)
                  = Z(t)*P_inf(t)*Z(t)' is non-singular.
        v       - p*NN*n matrix, innovation sequence.
        invF    - p*p*n matrix, inverse of F.
        K       - m*p*n matrix.
        L       - m*m*n matrix.
        RQ      - m*rr or m*rr*n matrix, R*Q.
        QRt     - rr*m or rr*m*n matrix, Q*R'.
        RQd     - true if one of R or Q is dynamic.
    Output:
        r           - m*N*(n+1) matrix, smoothing cumulant sequence.
        N           - m*m*(n+1) matrix, smoothing cumulant variance sequence.
        epshat      - p*N*n matrix, smoothed observation disturbance sequence.
        etahat      - rr*N*n matrix, smoothed state disturbance sequence.
        epsvarhat   - p*p*n matrix, smoothed observation disturbance variance.
        etavarhat   - rr*rr*n matrix, smoothed state disturbace variance.
\*=======================================================================*/
void disturbsmo(int p, int m, int rr, int n, int NN, const int *pknown, int pd, const int *pmiss,
        const double *H, int Hd, const double *Z, int Zd, const double *T, int Td, const double *Q, int Qd,
        int d, const int *Fns, const double *v, const double *invF, const double *K,
        const double *L, const double *RQ, const double *QRt, int RQd,
        double *r, double *N, double *epshat, double *etahat, double *epsvarhat, double *etavarhat);

/*=======================================================================*\
                            Smoothing cumulant
  -----------------------------------------------------------------------
    Input:
        p       - data dimension.
        m       - state dimension.
        n       - data size/length.
        pknown  - 1*n vector, dynamic data dimensions, can be NULL if
                  pd = false.
        pd      - true if data dimensions are dynamic.
        pmiss   - 1*n vector, data dimension after missing data is
                  removed, can be NULL if there's no missing data.
        Z       - p*m or p*m*n matrix, state to observation transform.
        Zd      - true if Z is dynamic (p*m*n).
        T       - m*m or m*m*n matrix, state update transform.
        Td      - true if T is dynamic (m*m*n).
        d       - length of non-zero diffuse variances.
        Fns     - 1*d vector, Fns[t] is true if F_inf(t)
                  = Z(t)*P_inf(t)*Z(t)' is non-singular.
        v       - p*NN*n matrix, innovation sequence.
        invF    - p*p*n matrix, inverse of F.
        K       - m*p*n matrix.
        L       - m*m*n matrix.
    Output:
        uuD     - p*p*n matrix, uuD(t) = u(t)*u(t)' - D(t).
        rrN     - m*m*n matrix, rrN(t) = r(t)*r(t)' - N(t).
\*=======================================================================*/
void smocum(int p, int m, int n, const int *pknown, int pd, const int *pmiss,
        const double *Z, int Zd, const double *T, int Td, int d, const int *Fns,
        const double *v, const double *invF, const double *K, const double *L,
        double *uuD, double *rrN);

/*=======================================================================*\
                          Missing data truncation
  -----------------------------------------------------------------------
    Input:
        p       - data dimension.
        m       - state dimension.
        n       - data size/length.
        N       - number of data sets.
        y       - p*N*n matrix, observation data.
        miss    - p*n matrix, true if corresponding element of y is missing.
        pknown  - 1*n vector, dynamic data dimensions, can be NULL if
                  pd = false.
        pd      - true if data dimensions are dynamic.
        H       - p*p or p*p*n matrix, observation disturbance variance.
        Hd      - true if H is dynamic (p*p*n).
        Z       - p*m or p*m*n matrix, state to observation transform.
        Zd      - true if Z is dynamic (p*m*n).
    Output:
        pmiss   - 1*n vector, data dimension after missing data is removed.
        imiss   - p*n matrix, index of non-missing data, there are pmiss[t]
                  indices at time point t.
        ymiss   - p*N*n matrix, observation data with missing data truncated.
        Hmiss   - p*p*n matrix, observation disturbance variance with
                  missing data locations truncated.
        Zmiss   - p*m*n matrix, state to observation transform with missing
                  data locations truncated.
\*=======================================================================*/
void truncate(int p, int m, int n, int N, const double *y, const int *miss,
    const int *pknown, int pd, const double *H, int Hd, const double *Z, int Zd,
    int *pmiss, int *imiss, double *ymiss, double *Hmiss, double *Zmiss);

/*=======================================================================*\
                               Untruncation
  -----------------------------------------------------------------------
    Input:
        p       - data dimension.
        n       - data size/length.
        N       - number of data sets.
        pknown  - 1*n vector, dynamic data dimensions, can be NULL if
                  pd = false.
        pd      - true if data dimensions are dynamic.
        pmiss   - 1*n vector, data dimension after missing data is removed.
        imiss   - p*n matrix, index of non-missing data, there are pmiss[t]
                  indices at time point t.
        vmiss   - p*N*n matrix, innovation with missing data location
                  truncated.
        Fmiss   - p*p*n matrix, innovation variance with missing data
                  location truncated.
    Output:
        v       - p*N*n matrix, innovation.
        F       - p*p*n matrix, innovation variance.
\*=======================================================================*/
void untruncate(int p, int n, int N, const int *pknown, int pd,
    const int *pmiss, const int *imiss, const double *vmiss, const double *Fmiss,
    double *v, double *F);

/*==================================*\
    Permute 2nd and 3rd dimensions
  ----------------------------------
    p   - 1st dimension.
    n   - 2nd dimension.
    N   - 3rd dimension.
    x   - p*n*N matrix.
    y   - p*N*n matrix.
\*==================================*/
void permute(int p, int n, int N, const double *x, double *y);

/*======================================================*\
    Generate random variables from normal distribution
  ------------------------------------------------------
    N   - number of samples to generate.
    x   - buffer to hold samples.
\*======================================================*/
void genrandn(int N, double *x);

/*======================================================*\
                   Transform covariance
  ------------------------------------------------------
    p   - number of variables.
    N   - number of samples.
    C   - p*p covariance matrix.
    u   - independent Gaussian random variables.
    x   - Gaussian random variables with covariance C.
\*======================================================*/
void sigma(int p, int N, const double *C, const double *u, double *x);

/*=======================================================================*\
                         Sample state space model
  -----------------------------------------------------------------------
    Input:
        p       - data dimension.
        m       - state dimension.
        r       - state disturbance dimension.
        n       - data size/length.
        N       - number of data sets.
        H       - p*p or p*p*n matrix, observation disturbance variance.
        Hd      - true if H is dynamic (p*p*n).
        Z       - p*m or p*m*n matrix, state to observation transform.
        Zd      - true if Z is dynamic (p*m*n).
        T       - m*m or m*m*n matrix, state update transform.
        Td      - true if T is dynamic (m*m*n).
        R       - m*r or m*r*n matrix, state disturbance transform.
        Rd      - true if R is dynamic (m*r*n).
        Q       - r*r or r*r*n matrix, state disturbance variance.
        Qd      - true if Q is dynamic (r*r*n).
        c       - m*1 or m*n matrix, state update constant.
        cd      - true if c is dynamic (m*n).
        a1      - m*1 vector, initial state vector.
        P1      - m*m matrix, initial state variance.
        samples - (m + (p + r)*n)*N independent Gaussian samples.
    Output:
        y       - p*N*n matrix, sampled observation.
        alpha   - m*N*n matrix, sampled state sequence.
        eps     - p*N*n matrix, sampled observation disturbance sequence.
        eta     - r*N*n matrix, sampled state disturbance sequence.
\*=======================================================================*/
void sample(int p, int m, int r, int n, int N, const double *H, int Hd,
    const double *Z, int Zd, const double *T, int Td, const double *R, int Rd,
    const double *Q, int Qd, const double *c, int cd, const double *a1,
    const double *P1, const double *samples, double *y, double *alpha,
    double *eps, double *eta);

/*=======================================================================*\
               Generate signal from state space components
  -----------------------------------------------------------------------
    Input:
        p       - data dimension.
        m       - state dimension.
        n       - data size/length.
        ncom    - number of signal components.
        mcom    - state dimension of each component, must sum to m.
        Z       - p*m or p*m*n matrix, state to observation transform.
        Zd      - true if Z is dynamic (p*m*n).
        alpha   - m*n matrix, state vector sequence.
    Output:
        ycom    - p*n*ncom matrix, signal for each component.
\*=======================================================================*/
void signal(int p, int m, int n, int ncom, const int *mcom,
    const double *Z, int Zd, const double *alpha, double *ycom);

#endif
