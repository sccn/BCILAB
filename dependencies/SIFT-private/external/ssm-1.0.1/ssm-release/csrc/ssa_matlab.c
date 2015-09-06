/**************************************\
 * Filename: ssa_matlab.c
 * Author: Jyh-Ying Peng (´^´¼·­)
 * Year: 2007
 * Version: 1.2
 * \**************************************/

#include "ssa_matlab.h"

double *g_y         = (double *)0;
int *g_miss         = (int *)0;
double *g_P1        = (double *)0;
double *g_P1_inf    = (double *)0;
double *g_ymiss     = (double *)0;
int *g_pmiss        = (int *)0;
int *g_imiss        = (int *)0;
double *g_Hmiss     = (double *)0;
double *g_Zmiss     = (double *)0;

void get_data(const mxArray *data, int *pp, int *pn, int *pN, double **y, int **miss)
{
    int i;
    
    double *ydata       = mxGetPr(data);
    
    if(pN) {
        int yndims          = mxGetNumberOfDimensions(data);
        const int *ydims    = mxGetDimensions(data);
        int p   = *pp       = ydims[0];
        int n   = *pn       = ydims[1];
        int N   = *pN       = (yndims > 2) ? ydims[2] : 1;
        
        *miss   = (int *)0;
        for(i=0; i<p*N*n; i++)
            if(mxIsNaN(ydata[i])) {
                *miss   = g_miss = (int *)mxCalloc(p*n, sizeof(int));
                break;
            }
        
        if(N > 1) {
            int j, t;
            double *yiter   = *y = g_y = (double *)mxCalloc(p*N*n, sizeof(double));
            for(t=0; t<n; t++, yiter += p*N)
                for(i=0; i<p; i++) {
                    if(*miss) (*miss)[i+t*p] = false;
                    for(j=0; j<N; j++) {
                        yiter[i+j*p]    = ydata[i+t*p+j*p*n];
                        if(*miss && !(*miss)[i+t*p] && mxIsNaN(yiter[i+j*p])) (*miss)[i+t*p] = true;
                    }
                }
        }
        else {
            *y  = ydata;
            if(*miss) for(i=0; i<p*n; i++) (*miss)[i] = mxIsNaN((*y)[i]);
        }
    }
    else {
        int p   = *pp       = mxGetM(data);
        int n   = *pn       = mxGetN(data);
        
        *miss  = (int *)0;
        for(i=0; i<p*n; i++)
            if(mxIsNaN(ydata[i])) {
                *miss   = g_miss = (int *)mxCalloc(p*n, sizeof(int));
                break;
            }
        
        *y  = ydata;
        if(*miss) for(i=0; i<p*n; i++) (*miss)[i] = mxIsNaN((*y)[i]);
    }
}

void get_M(const mxArray *model[], int *pp, double **M, int *Md)
{
    if(pp) *pp = mxGetM(model[0]);
    *M  = mxGetPr(model[0]);
    if(Md) *Md = (int)mxGetScalar(model[1]);
}

void get_P1(const mxArray *model, int *pm, double **P1, double **P1_inf)
{
    int m, i, init;
    
    m       = mxGetM(model);
    if(pm) *pm = m;
    *P1     = mxGetPr(model);
    
    /*** Determine initialization ***/
    for(i=0, init=false; i<m*m; i++) if(mxIsInf((*P1)[i])) { init = true; break; }
    
    if(init) {
        /*** Make a copy first ***/
        double *P1new   = g_P1 = (double *)mxCalloc(m*m, sizeof(double));
        *P1_inf         = g_P1_inf = (double *)mxCalloc(m*m, sizeof(double));
        for(i=0; i<m*m; i++) {
            P1new[i]    = (*P1)[i];
            if(mxIsInf(P1new[i])) {
                P1new[i]        = 0;
                (*P1_inf)[i]    = 1;
            }
            else (*P1_inf)[i]   = 0;
        }
        *P1             = P1new;
    }
    else *P1_inf = (double *)0;
}

void truncate_yHZ(int p, int m, int n, int N, double **py, const int *miss,
int **ppmiss, int **pimiss, double **pH, int *pHd, double **pZ, int *pZd)
{
    double *y, *H, *Z, *ymiss, *Hmiss, *Zmiss;
    
    *ppmiss = g_pmiss = (int *)mxCalloc(n, sizeof(int));
    *pimiss = g_imiss = (int *)mxCalloc(p*n, sizeof(int));
    y       = py ? *py : (double *)0;
    H       = pH ? *pH : (double *)0;
    Z       = pZ ? *pZ : (double *)0;
    ymiss   = py ? (g_ymiss = (double *)mxCalloc(p*N*n, sizeof(double))) : (double *)0;
    Hmiss   = pH ? (g_Hmiss = (double *)mxCalloc(p*p*n, sizeof(double))) : (double *)0;
    Zmiss   = pZ ? (g_Zmiss = (double *)mxCalloc(p*m*n, sizeof(double))) : (double *)0;
    truncate(p, m, n, N, y, miss, (int *)0, false, H, *pHd, Z, *pZd, *ppmiss, *pimiss, ymiss, Hmiss, Zmiss);
    if(py) *py = ymiss;
    if(pH) { *pH = Hmiss; *pHd = true; }
    if(pZ) { *pZ = Zmiss; *pZd = true; }
}

void permute_data(mxArray *data, int p, int N, int n)
{
    int dims[3];
    
    if(N > 1) {
        int i;
        double *pdata   = mxGetPr(data);
        double *temp    = (double *)mxCalloc(p*n*N, sizeof(double));
        
        permute(p, N, n, pdata, temp);
        for(i=0; i<p*n*N; i++) pdata[i] = temp[i];
        mxFree(temp);
        
        dims[0] = p; dims[1] = n; dims[2] = N;
        mxSetDimensions(data, dims, 3);
    }
    else {
        dims[0] = p; dims[1] = n;
        mxSetDimensions(data, dims, 2);
    }
}

void cleanup()
{
    if(g_y) { mxFree(g_y); g_y = (double *)0; }
    if(g_miss) { mxFree(g_miss); g_miss = (int *)0; }
    if(g_P1) { mxFree(g_P1); g_P1 = (double *)0; }
    if(g_P1_inf) { mxFree(g_P1_inf); g_P1_inf = (double *)0; }
    if(g_ymiss) { mxFree(g_ymiss); g_ymiss = (double *)0; }
    if(g_pmiss) { mxFree(g_pmiss); g_pmiss = (int *)0; }
    if(g_imiss) { mxFree(g_imiss); g_imiss = (int *)0; }
    if(g_Hmiss) { mxFree(g_Hmiss); g_Hmiss = (double *)0; }
    if(g_Zmiss) { mxFree(g_Zmiss); g_Zmiss = (double *)0; }
}
