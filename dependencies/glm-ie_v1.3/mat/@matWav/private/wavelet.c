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

#include "wavelet.h"
#include <stdlib.h>
#include <string.h>

/* 1d downsampling using high-pass filter 
 * IN:  x - vector of (even) length nx
 *      f - high-pass filter of length nf
 * OUT: y - vector of length nx/2
 */
void downhi( double *x, int nx,  double *f, int nf,  double *y ) {
    int nx2, nf2, h, i, j;  double s;

    nx2 = nx/2;
    nf2 = nf/2-1;  if(2*nf2+1-(nf-1) < 0) nf2++;
    
    for(i=nf2; i<nx2; i++) {           /* Calculate y[nf2] till y[nx2-1] */
        s = 0.0;
        for(h=0; h<nf; h++)  s += f[h]*x[2*i+1-h];
        y[i] = s;
    }
    if(nf2 > nx2) nf2 = nx2;
    
    for(i=0; i<nf2; i++) {             /* Calculate y[0] till y[nf2-1] */
        s = 0.0;
        j = 2*i+1;
        for( h=0; h<nf; h++) {
            if(j < 0) j += nx;
            s += f[h]*x[j];
            --j;
        }
        y[i] = s;
    }
}


/* 1d downsampling using low-pass filter
 * IN:  x - vector of (even) length nx
 *      f - low-pass filter of length nf
 * OUT: y - vector of length nx/2
 */
void downlo(double *x, int nx, double *f, int nf, double *y) {
    int nx2, nf2, dn2, h, i, j;  double s;

    nx2 = nx/2;
    nf2 = nf/2; 
    dn2 = nx2-nf2;
    if(2*dn2 + (nf2-1) >= nx) --dn2;
    if(dn2 < 0) dn2 = -1;
    
    for(i=0; i<=dn2; i++){             /* Calculate y[0] till y[nx2-nf2-1] */
        s = 0.0;
        for(h=0; h<nf; h++)  s += f[h]*x[2*i+h];
        y[i] = s;
    }

    for(i=dn2+1; i<nx2; i++){          /* Calculate y[nx2-nf2] till y[nf2-1] */
        s = 0.0;
        j = 2*i;
        for(h=0; h<nf; h++){
            if(j >= nx) j -= nx;
            s += f[h]*x[j];
            j++;
        }
        y[i] = s;
    }
}


/* 1d upsampling using high-pass filter
 * IN:  x - vector of length nx
 *      f - high-pass filter of length nf
 * OUT: y - vector of length 2*nx
 */
void uphi(double *x, int nx, double *f, int nf, double *y) {
    int even, odd, min, h, i, j;  double s;

    even = (nf+1)/2; odd  = nf/2; min = nx-even; if(min < 0) min = 0;

    for(i=0; i+even<nx; i++){         /* Calculate y[0] till y[2*nf-1] */
        s = 0.0;
        for(h=0; h<even; h++)  s += f[2*h  ]*x[i+h];
        y[2*i+1] = s;
        
        s = 0.0;
        for(h=0; h<odd;  h++)  s += f[2*h+1]*x[i+h];
        y[2*i  ] = s;
    }

    for(i=min; i<nx; i++){            /* Calculate y[2*nf] till y[2*nx-1] */
        s = 0.0;
        j = i;
        for(h=0; h<even; h++){
            if(j >= nx) j -= nx;
            s += f[2*h]*x[j];
            j++;
        }
        y[2*i+1] = s;
        
        s = 0.0;
        j = i;
        for(h=0; h<odd; h++){
            if(j >= nx) j -= nx;
            s += f[2*h+1]*x[j];
            j++;
        }
        y[2*i] = s;
    }
}


/* 1d upsampling using low-pass filter
 * IN:  x - vector of length nx
 *      f - low-pass filter of length nf
 * OUT: y - vector of length 2*nx
 */
void uplo(double *x, int nx, double *f, int nf, double *y) {
    int even, odd, max, h, i, j;  double s;

    even = (nf+1)/2; odd = nf/2; max = even; if(max > nx) max = nx;
    
    for(i=even; i<nx; i++){            /* Calculate y[0] till y[2*nf-1] */
        s = 0.0;
        for(h=0; h<even; h++) s += f[2*h  ]*x[i-h];
        y[2*i  ] = s;
        
        s = 0.0;
        for(h=0; h<odd; h++) s += f[2*h+1]*x[i-h];
        y[2*i+1] = s;				
    }
    
    for(i=0; i<max; i++){              /* Calculate y[2*nf] till y[2*nx-1] */
        s = 0.0;
        j = i;
        for(h=0; h<even; h++){
            if(j < 0) j += nx;
            s += f[2*h]*x[j];
            --j;
        }
        y[2*i] = s;
        
        s = 0.0;
        j = i;
        for(h=0; h<odd; h++){
            if(j < 0) j += nx;
            s += f[2*h+1]*x[j];
            --j;
        }
        y[2*i+1] = s;
    }
}




/* Wavelet Transform in D dimensions and L levels */
void  wtn(double *X, int *N, int D, int L, double *f, int Nf, double *W, double *tmp) { 

    double *tmplo, *tmphi, *fm;
    int *Nl, *Noff, *stride;
    int l, d, dd, k, s, length, offset, numel;
    /* Size of the subcube in the current level */
    Nl     = (int*) malloc(sizeof(int)*D);   for(d=0; d<D; d++) Nl[d]=N[d];
    /* Index in the subcube in the current level, cube of valid offsets */
    Noff   = (int*) malloc(sizeof(int)*D);
    /* Strides or 1d neighbor distances */
    stride = (int*) malloc(sizeof(int)*D);
    for(d=0,dd=1; d<D; d++) {stride[d]=dd; dd*=N[d];}
    
    for(numel=1,d=0; d<D; d++) numel*=N[d];       /* Total number of elements */   
    memcpy(W, X, numel*sizeof(double));               /* Copy input to output */    
    /* Construct mirrored filter fm from the filter f */
    fm = (double*) malloc( sizeof(double)*Nf );
    for(k=0,s=1; k<Nf; k++,s*=-1) fm[k] = s*f[k];
    
    for(l=0; l<L; l++) {                             /* Iterate over levels l */
        for(d=0; d<D; d++) {                     /* Iterate over dimensions d */
    		tmplo = &tmp[Nl[d]]; tmphi = &tmp[3*Nl[d]/2];   /* Temp. pointers */
			length = Nl[d]; Nl[d] = 1;              /* Length of the 1d slice */  
            
            offset=0; for(dd=0; dd<D; dd++) Noff[dd]=0;  /* Start at (0,..,0) */
            while(Noff[D-1]<Nl[D-1]) {   /* Go through cube N, stay inside Nl */

                /* 1d transform: only (offset, stride, length) are needed */
                for(k=0;k<length;k++) tmp[k]=W[offset+stride[d]*k];   /*Read  */
				downlo(tmp, length, f , Nf, tmplo);   /* 1d Transform along d */
                downhi(tmp, length, fm, Nf, tmphi);   /* 1d Transform along d */
				for(k=0;k<length;k++) W[offset+stride[d]*k]=tmplo[k]; /*Write */
                
                Noff[0]++;                        /* Go one step forward in N */
                dd=0;
                while(Noff[dd]>=Nl[dd] && dd<D-1){ /* Carry over if not in Nl */
                    Noff[dd+1]++;
                    Noff[dd]=0;
                    dd++;
                }
                if(dd==0) offset++; else                   /* Increase offset */
                    for(dd=0,offset=0; dd<D; dd++) offset+=Noff[dd]*stride[dd];
            }

        Nl[d] = length; /* Along the current dimension, the subcube has size 1*/
		}                                        /* Iterate over dimensions d */
        for(d=0; d<D; d++) Nl[d] /= 2;                   /* Halve the subcube */
    }                                                /* Iterate over levels l */   
    free(fm); free(Nl); free(Noff);                            /* Free memory */
}



/* Inverse Wavelet Transform in D dimensions and L levels */
void iwtn(double *W, int *N, int D, int L, double *f, int Nf, double *X, double *tmp){

    double *tmplo, *tmphi, *tmptop, *fm;
    int *Nl, *Noff, *stride;
    int l, d, dd, j, k, s, length, offset, numel;
    /* Size of the subcube in the current level */
    for(l=0,j=1; l<L-1; l++) j *= 2;      /* Scaling from N to coarsest level */
    Nl     = (int*) malloc(sizeof(int)*D); for(d=0; d<D; d++) Nl[d]=N[d]/j;
    /* Index in the subcube in the current level, cube of valid offsets */
    Noff   = (int*) malloc(sizeof(int)*D);
    /* Strides or 1d neighbor distances */
    stride = (int*) malloc(sizeof(int)*D);
    for(d=0,dd=1; d<D; d++) {stride[d]=dd; dd*=N[d];}
    
    for(numel=1,d=0; d<D; d++) numel*=N[d];       /* Total number of elements */   
    memcpy(X, W, numel*sizeof(double));               /* Copy input to output */
    /* Construct mirrored filter fm from the filter f */
    fm = (double*) malloc( sizeof(double)*Nf );
    for(k=0,s=1; k<Nf; k++,s*=-1) fm[k] = s*f[k];;     
        
    for(l=L-1; l>=0; l--) {                          /* Iterate over levels l */
        for(d=D-1; d>=0; d--) {                  /* Iterate over dimensions d */
			length = Nl[d]/2; Nl[d] = 1;            /* Length of the 1d slice */  
            tmplo = &tmp[2*length]; 
            tmphi = &tmp[3*length]; tmptop = &tmp[4*length]; /* Temp pointers */
            
            offset=0; for(dd=0; dd<D; dd++) Noff[dd]=0;  /* Start at (0,..,0) */
            while(Noff[D-1]<Nl[D-1]) {   /* Go through cube N, stay inside Nl */

                /* 1d transform: only (offset, stride, length) are needed */
                for(k=0;k<2*length;k++) tmplo[k]=X[offset+stride[d]*k]; /*Read*/
                uplo(tmplo, length, f,  Nf, tmp);
                uphi(tmphi, length, fm, Nf, tmptop);
                for(k=0;k<2*length;k++) tmp[k]=tmp[k]+tmptop[k];         /*Add*/
				for(k=0;k<2*length;k++) X[offset+stride[d]*k]=tmp[k];  /*Write*/
                
                Noff[0]++;                        /* Go one step forward in N */
                dd=0;
                while(Noff[dd]>=Nl[dd] && dd<D-1){ /* Carry over if not in Nl */
                    Noff[dd+1]++;
                    Noff[dd]=0;
                    dd++;
                }
                if(dd==0) offset++; else                   /* Increase offset */
                    for(dd=0,offset=0; dd<D; dd++) offset+=Noff[dd]*stride[dd];
            }

        Nl[d]=2*length; /* Along the current dimension, the subcube has size 1*/
		}                                        /* Iterate over dimensions d */
        for(d=0; d<D; d++) Nl[d] *= 2;                  /* Double the subcube */
    }                                                /* Iterate over levels l */
    free(fm); free(Nl); free(Noff);                            /* Free memory */
}
