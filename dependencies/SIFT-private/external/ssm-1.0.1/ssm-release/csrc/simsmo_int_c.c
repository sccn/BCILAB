#include "ssa_matlab.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /***************************************************************************\
     MATLAB
     [alphatilde epstilde etatilde alphaplus]
        = simsmo_int(y, N, H, Hd, Z, Zd, T, Td, R, Rd, Q, Qd, c, cd, a1, P1,
                        antithetic, tol0, tolP, inv_method, postproc)
     simsmo_int('seed', S)
    \***************************************************************************/
    
    if(mxIsChar(prhs[0])) {
        int length          = mxGetNumberOfElements(prhs[0]) + 1;
        char *str           = (char *)mxCalloc(length, sizeof(char));
        
        mxGetString(prhs[0], str, length);
        if(strcmp(str, "seed") == 0) {
            int n               = mxGetNumberOfElements(prhs[1]);
            unsigned long *seed = (unsigned long *)mxGetData(prhs[1]);
            if(n == 1) init_genrand(*seed);
            else init_by_array(seed, n);
        }
        
        mxFree(str);
        return;
    }
    else {
        /*** Input variables ***/
        int p, m, r, n, N, Nsamp;
        double *y;
        int *miss, *pmiss, *imiss;
        int Hd, Zd, Td, Rd, Qd, cd, sta;
        double *H, *Z, *T, *R, *Q, *c, *a1, *P1, *P1_inf;
        double tol0, tolP;
        int antithetic, inv_method, postproc;
        
        /*** Output variables ***/
        int dims[3];
        double *alphatilde, *epstilde, *etatilde, *alphaplus;
        
        /*** Buffer variables ***/
        double *yplus, *epsplus, *etaplus, *samples;
        double *alphaplushat, *epsplushat, *etaplushat;
        double *alphahat, *epshat, *etahat;
        int d, *Fns;
        double *v, *invF, *F2, *K, *L, *L1, *QRt, *RQRt;
        
        /*** Temporary variables ***/
        int i, j, t;
        
        /*** Retrieve data ***/
        get_data(prhs[0], &p, &n, (int *)0, &y, &miss);
        
        /*** Retrieve model ***/
        get_M(prhs + 2, &p, &H, &Hd);
        get_M(prhs + 4, (int *)0, &Z, &Zd);
        get_M(prhs + 6, &m, &T, &Td);
        get_M(prhs + 8, (int *)0, &R, &Rd);
        get_M(prhs + 10, &r, &Q, &Qd);
        get_M(prhs + 12, (int *)0, &c, &cd);
        get_M(prhs + 14, (int *)0, &a1, (int *)0);
        get_P1(prhs[15], (int *)0, &P1, &P1_inf);
        sta     = !Hd && !Zd && !Td && !Rd && !Qd;
        
        /*** Retrieve analysis settings ***/
        N           = (int)mxGetScalar(prhs[1]);
        antithetic  = (int)mxGetScalar(prhs[16]);
        Nsamp       = (antithetic > 0) ? (N+1)/2 : N;
        tol0        = mxGetScalar(prhs[17]);
        tolP        = mxGetScalar(prhs[18]);
        inv_method  = (int)mxGetScalar(prhs[19]);
        postproc    = (int)mxGetScalar(prhs[20]);
        
        /*** Allocate buffers ***/
        yplus           = (double *)mxCalloc(p*Nsamp*n, sizeof(double));
        epsplus         = (double *)mxCalloc(p*Nsamp*n, sizeof(double));
        etaplus         = (double *)mxCalloc(r*Nsamp*n, sizeof(double));
        alphaplushat    = (double *)mxCalloc(m*Nsamp*n, sizeof(double));
        epsplushat      = (double *)mxCalloc(p*Nsamp*n, sizeof(double));
        etaplushat      = (double *)mxCalloc(r*Nsamp*n, sizeof(double));
        alphahat        = (double *)mxCalloc(m*n, sizeof(double));
        epshat          = (double *)mxCalloc(p*n, sizeof(double));
        etahat          = (double *)mxCalloc(r*n, sizeof(double));
        Fns             = (int *)mxCalloc(n, sizeof(int));
        v               = (double *)mxCalloc(p*Nsamp*n, sizeof(double));
        invF            = (double *)mxCalloc(p*p*n, sizeof(double));
        F2              = (double *)mxCalloc(p*p*n, sizeof(double));
        K               = (double *)mxCalloc(m*p*n, sizeof(double));
        L               = (double *)mxCalloc(m*m*n, sizeof(double));
        L1              = (double *)mxCalloc(m*m*n, sizeof(double));
        QRt             = (double *)((Rd || Qd) ? mxCalloc(r*m*n, sizeof(double)) : mxCalloc(r*m, sizeof(double)));
        RQRt            = (double *)((Rd || Qd) ? mxCalloc(m*m*n, sizeof(double)) : mxCalloc(m*m, sizeof(double)));
        
        /*** Allocate output variables ***/
        if(postproc) {
            dims[0]     = m; dims[1] = n; dims[2] = N;
            alphatilde  = mxGetPr(plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL));
            dims[0]     = p; dims[1] = n; dims[2] = N;
            epstilde    = mxGetPr(plhs[1] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL));
            dims[0]     = r; dims[1] = n; dims[2] = N;
            etatilde    = mxGetPr(plhs[2] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL));
        }
        else {
            dims[0]     = m; dims[1] = N; dims[2] = n;
            alphatilde  = mxGetPr(plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL));
            dims[0]     = p; dims[1] = N; dims[2] = n;
            epstilde    = mxGetPr(plhs[1] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL));
            dims[0]     = r; dims[1] = N; dims[2] = n;
            etatilde    = mxGetPr(plhs[2] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL));
        }
        dims[0]     = m; dims[1] = N; dims[2] = n;
        alphaplus   = mxGetPr(plhs[3] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL));
        
        /*** Unconditional sampling ***/
        samples = (double *)mxCalloc((m+(p+r)*n)*Nsamp, sizeof(double));
        genrandn((m+(p+r)*n)*Nsamp, samples);
        sample(p, m, r, n, Nsamp, H, Hd, Z, Zd, T, Td, R, Rd, Q, Qd, c, cd, a1, P1, samples, yplus, alphaplus, epsplus, etaplus);
        mxFree(samples);
        
        /*** Batch fast smoothing ***/
        kalman(p, m, r, n, Nsamp, yplus, (int *)0, false, (int *)0, sta, H, Hd, Z, Zd, T, Td, R, Rd, Q, Qd, c, cd,
            a1, P1, P1_inf, (double *)0, (double *)0, (double *)0, &d, Fns, v, (double *)0, invF, F2, K, (double *)0,
            L, L1, (double *)0, (double *)0, (double *)0, QRt, RQRt, tol0, tolP, inv_method);
        fastsmo(p, m, r, n, Nsamp, (int *)0, false, (int *)0, H, Hd, Z, Zd, T, Td, R, Rd, c, cd, a1, P1, P1_inf,
            d, Fns, v, invF, F2, K, L, L1, QRt, RQRt, Rd || Qd, (double *)0, (double *)0, alphaplushat, epsplushat,
            etaplushat);

        /*** Truncate Z and H w.r.t. missing data ***/
        if(miss) truncate_yHZ(p, m, n, 1, &y, miss, &pmiss, &imiss, &H, &Hd, &Z, &Zd);
        else { pmiss = (int *)0; imiss = (int *)0; }
        
        /*** Fast smoothing ***/
        kalman(p, m, r, n, 1, y, (int *)0, false, pmiss, sta, H, Hd, Z, Zd, T, Td, R, Rd, Q, Qd, c, cd, a1, P1, P1_inf,
            (double *)0, (double *)0, (double *)0, &d, Fns, v, (double *)0, invF, F2, K, (double *)0,
            L, L1, (double *)0, (double *)0, (double *)0, QRt, RQRt, tol0, tolP, inv_method);
        fastsmo(p, m, r, n, 1, (int *)0, false, pmiss, H, Hd, Z, Zd, T, Td, R, Rd, c, cd, a1, P1, P1_inf,
            d, Fns, v, invF, F2, K, L, L1, QRt, RQRt, Rd || Qd, (double *)0, (double *)0, alphahat, epshat, etahat);

        /*** Conditional samples ***/
        if(antithetic >= 1) {
            if(postproc) {
                if(N > 1) {
                    for(t=0; t<n; t++) {
                        for(i=0; i<m; i++)
                            for(j=0; j<N/2; j++) {
                                alphatilde[i+(t+2*j*n)*m]       = alphahat[i+t*m] - alphaplushat[i+(j+t*Nsamp)*m] + alphaplus[i+(j+t*Nsamp)*m];
                                alphatilde[i+(t+(2*j+1)*n)*m]   = alphahat[i+t*m] + alphaplushat[i+(j+t*Nsamp)*m] - alphaplus[i+(j+t*Nsamp)*m];
                            }
                        for(i=0; i<p; i++)
                            for(j=0; j<N/2; j++) {
                                epstilde[i+(t+2*j*n)*p]         = epshat[i+t*p] - epsplushat[i+(j+t*Nsamp)*p] + epsplus[i+(j+t*Nsamp)*p];
                                epstilde[i+(t+(2*j+1)*n)*p]     = epshat[i+t*p] + epsplushat[i+(j+t*Nsamp)*p] - epsplus[i+(j+t*Nsamp)*p];
                            }
                        for(i=0; i<r; i++)
                            for(j=0; j<N/2; j++) {
                                etatilde[i+(t+2*j*n)*r]         = etahat[i+t*r] - etaplushat[i+(j+t*Nsamp)*r] + etaplus[i+(j+t*Nsamp)*r];
                                etatilde[i+(t+(2*j+1)*n)*r]     = etahat[i+t*r] + etaplushat[i+(j+t*Nsamp)*r] - etaplus[i+(j+t*Nsamp)*r];
                            }
                    }
                }
                if(N % 2) {
                    for(t=0; t<n; t++) {
                        for(i=0; i<m; i++)
                            alphatilde[i+(t+(N-1)*n)*m] = alphahat[i+t*m] - alphaplushat[i+(Nsamp-1+t*Nsamp)*m] + alphaplus[i+(Nsamp-1+t*Nsamp)*m];
                        for(i=0; i<p; i++)
                            epstilde[i+(t+(N-1)*n)*p]   = epshat[i+t*p] - epsplushat[i+(Nsamp-1+t*Nsamp)*p] + epsplus[i+(Nsamp-1+t*Nsamp)*p];
                        for(i=0; i<r; i++)
                            etatilde[i+(t+(N-1)*n)*r]   = etahat[i+t*r] - etaplushat[i+(Nsamp-1+t*Nsamp)*r] + etaplus[i+(Nsamp-1+t*Nsamp)*r];
                    }
                }
                permute_data(plhs[3], m, N, n);
            }
            else {
                if(N > 1) {
                    for(t=0; t<n; t++) {
                        for(i=0; i<m; i++)
                            for(j=0; j<N/2; j++) {
                                alphatilde[i+(2*j+t*N)*m]   = alphahat[i+t*m] - alphaplushat[i+(j+t*Nsamp)*m] + alphaplus[i+(j+t*Nsamp)*m];
                                alphatilde[i+(2*j+1+t*N)*m] = alphahat[i+t*m] + alphaplushat[i+(j+t*Nsamp)*m] - alphaplus[i+(j+t*Nsamp)*m];
                            }
                        for(i=0; i<p; i++)
                            for(j=0; j<N/2; j++) {
                                epstilde[i+(2*j+t*N)*p]     = epshat[i+t*p] - epsplushat[i+(j+t*Nsamp)*p] + epsplus[i+(j+t*Nsamp)*p];
                                epstilde[i+(2*j+1+t*N)*p]   = epshat[i+t*p] + epsplushat[i+(j+t*Nsamp)*p] - epsplus[i+(j+t*Nsamp)*p];
                            }
                        for(i=0; i<r; i++)
                            for(j=0; j<N/2; j++) {
                                etatilde[i+(2*j+t*N)*r]     = etahat[i+t*r] - etaplushat[i+(j+t*Nsamp)*r] + etaplus[i+(j+t*Nsamp)*r];
                                etatilde[i+(2*j+1+t*N)*r]   = etahat[i+t*r] + etaplushat[i+(j+t*Nsamp)*r] - etaplus[i+(j+t*Nsamp)*r];
                            }
                    }
                }
                if(N % 2) {
                    for(t=0; t<n; t++) {
                        for(i=0; i<m; i++)
                            alphatilde[i+(N-1+t*N)*m]   = alphahat[i+t*m] - alphaplushat[i+(Nsamp-1+t*Nsamp)*m] + alphaplus[i+(Nsamp-1+t*Nsamp)*m];
                        for(i=0; i<p; i++)
                            epstilde[i+(N-1+t*N)*p]     = epshat[i+t*p] - epsplushat[i+(Nsamp-1+t*Nsamp)*p] + epsplus[i+(Nsamp-1+t*Nsamp)*p];
                        for(i=0; i<r; i++)
                            etatilde[i+(N-1+t*N)*r]     = etahat[i+t*r] - etaplushat[i+(Nsamp-1+t*Nsamp)*r] + etaplus[i+(Nsamp-1+t*Nsamp)*r];
                    }
                }
            }
        }
        else {
            if(postproc) {
                for(t=0; t<n; t++) {
                    for(j=0; j<N; j++) {
                        for(i=0; i<m; i++)
                            alphatilde[i+(t+j*n)*m] = alphahat[i+t*m] - alphaplushat[i+(j+t*N)*m] + alphaplus[i+(j+t*N)*m];
                        for(i=0; i<p; i++)
                            epstilde[i+(t+j*n)*p]   = epshat[i+t*p] - epsplushat[i+(j+t*N)*p] + epsplus[i+(j+t*N)*p];
                        for(i=0; i<r; i++)
                            etatilde[i+(t+j*n)*r]   = etahat[i+t*r] - etaplushat[i+(j+t*N)*r] + etaplus[i+(j+t*N)*r];
                    }
                }
                permute_data(plhs[3], m, N, n);
            }
            else {
                for(t=0; t<n; t++) {
                    for(j=0; j<N; j++) {
                        for(i=0; i<m; i++)
                            alphatilde[i+(j+t*N)*m] = alphahat[i+t*m] - alphaplushat[i+(j+t*N)*m] + alphaplus[i+(j+t*N)*m];
                        for(i=0; i<p; i++)
                            epstilde[i+(j+t*N)*p]   = epshat[i+t*p] - epsplushat[i+(j+t*N)*p] + epsplus[i+(j+t*N)*p];
                        for(i=0; i<r; i++)
                            etatilde[i+(j+t*N)*r]   = etahat[i+t*r] - etaplushat[i+(j+t*N)*r] + etaplus[i+(j+t*N)*r];
                    }
                }
            }
        }
        
        /*** Cleanup ***/
        cleanup();
        mxFree(yplus);
        mxFree(epsplus);
        mxFree(etaplus);
        mxFree(alphaplushat);
        mxFree(epsplushat);
        mxFree(etaplushat);
        mxFree(alphahat);
        mxFree(epshat);
        mxFree(etahat);
        mxFree(Fns);
        mxFree(v);
        mxFree(invF);
        mxFree(F2);
        mxFree(K);
        mxFree(L);
        mxFree(L1);
        mxFree(QRt);
        mxFree(RQRt);
    }
}
