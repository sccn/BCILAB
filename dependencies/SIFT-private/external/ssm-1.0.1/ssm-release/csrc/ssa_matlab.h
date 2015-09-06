/**************************************\
    Filename: ssa_matlab.h
    Author: Jyh-Ying Peng (´^´¼·­)
    Year: 2007
    Version: 1.2
\**************************************/

#ifndef SSA_MATLAB_INCLUDED_20070215
#define SSA_MATLAB_INCLUDED_20070215

#include "mex.h"
#include "ssa.h"

void get_data(const mxArray *data, int *pp, int *pn, int *pN, double **y, int **miss);

void get_M(const mxArray *model[], int *pp, double **M, int *Md);

void get_P1(const mxArray *model, int *pm, double **P1, double **P1_inf);

void truncate_yHZ(int p, int m, int n, int N, double **py, const int *miss,
    int **ppmiss, int **pimiss, double **pH, int *pHd, double **pZ, int *pZd);

void permute_data(mxArray *data, int p, int N, int n);

void cleanup();

#endif
