/* FAST WAVELET TRANSFORM in D dimensions
 * X      - array of size N[0]xN[1]x..xN[D-1] containing the data
 * N      - D integers containing the size of X/W with N[i]=2^L*n[i]
 * L      - number of levels in the pyramid
 * f      - 1d low pass quadrature mirror filter 
 * Nf     - length of the 1d low pass quadrature mirror filter 
 * W      - the double array containing the wavelet coefficients same size as X
 * tmp    - temporary array of size 3*max(N[i]) for auxiliary calculations
 *
 *  Written by Hannes Nickisch, 17 Jan 2008
 * 
 * 
 *  The 4 auxiliary routines {up,down}{hi,lo} were inspired by code taken from 
 *  WaveLab Version 850.
 *  http://www-stat.stanford.edu/~wavelab/Wavelab_850/download.html
 *
 */

/* Wavelet Transform */
void  wtn(double *X, int *N, int D, int L, double *f, int Nf, double *W, double *tmp);

/* Inverse Wavelet Transform */
void iwtn(double *W, int *N, int D, int L, double *f, int Nf, double *X, double *tmp);
