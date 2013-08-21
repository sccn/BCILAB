/******************************************************************************
 *
 * MATLAB (R) is a trademark of The Mathworks (R) Corporation
 *
 * Function:    bsxarg
 * Filename:    bsxarg.c
 * Programmer:  James Tursa
 * Version:     1.0
 * Date:        February 09, 2008
 * Copyright:   (c) 2008 by James Tursa, All Rights Reserved
 * Permission:  Permission is granted to freely distribute and use this code
 *              as long as the header information is included.
 *
 * Returns the singleton expanded arrays that are virtually used in bsxfun.
 *
 * bsxfun performs binary operations on input arrays, where singleton
 * dimensions are virtually expanded to perform the operation. bsxarg will
 * actually do the singleton expansion and return the expanded arrays.
 *
 * Building:
 *
 * First review the bit length defines below and make any necessary changes
 * to match your particular system. You may need to change the bit16, bit32,
 * and/or bit64 definitions. Basically, you simply need to have an integer
 * type that is the specified bit length so that the copy method I use will
 * work.
 *
 * >> mex -setup
 *   (then follow instructions to select a C / C++ compiler of your choice)
 * >> mex bsxarg.c
 *
 * Syntax
 *
 * [C D] = bsxarg(A,B)
 *
 * Description
 *
 * C = the expanded version of A.
 * D = the expanded version of B.
 *
 * Each dimension of A and B must either be equal to each other, or equal to 1.
 * Whenever a dimension of A or B is singleton (equal to 1), the array is 
 * virtually replicated along the dimension to match the other array. The array
 * may be diminished if the corresponding dimension of the other array is 0.
 *
 * The size of the output arrays C and D are equal to:
 *
 * max(size(A),size(B)).*(size(A)>0 & size(B)>0).
 *
 * Examples
 * 
 *  >> A = 10*rand(2,1,2)
 *  A(:,:,1) =
 *    9.3181
 *    4.6599
 *  A(:,:,2) =
 *    4.1865
 *    8.4622
 *
 *  >> B = 10*rand(1,3)
 *  B =
 *    5.2515    2.0265    6.7214
 *
 *  >> [C D] = bsxarg(A,B)
 *  C(:,:,1) =
 *   9.3181    9.3181    9.3181
 *   4.6599    4.6599    4.6599
 *  C(:,:,2) =
 *   4.1865    4.1865    4.1865
 *   8.4622    8.4622    8.4622
 *
 *  D(:,:,1) =
 *   5.2515    2.0265    6.7214
 *   5.2515    2.0265    6.7214
 *  D(:,:,2) =
 *   5.2515    2.0265    6.7214
 *   5.2515    2.0265    6.7214
 *
 * User's with older versions of MATLAB that do not have the bsxfun intrinsic
 * available to them can use this simple m-file to get that capability:
 *
 *  function C = bsxfun(fun,A,B)
 *  if( nargin ~= 3 )
 *      error('Need 3 arguments for bsxfun')
 *  end
 *  [AX BX] = bsxarg(A,B);
 *  if( ischar(fun) )
 *      C = eval([fun '(AX,BX)']);
 *  else
 *      C = fun(AX,BX);
 *  end
 *  return
 *  end
 *
 ********************************************************************************/

// Bit-length definitions. These may need to be adjusted for your system.

#define bit8 char
#define bit16 short
#define bit32 long
#define bit64 long long

// Include headers & defines

#include <string.h>
#include "mex.h"
#include "matrix.h"

#ifndef mwSize
#define mwSize int
#endif

// Macro for extracting the singleton expansion element of an array

#define AX(x,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12) x[i1*A1+i2*A2+i3*A3+i4*A4+i5*A5+i6*A6+i7*A7+i8*A8+i9*A9+i10*A10+i11*A11+i12*A12]

// Macro for MAX

#define MAX(X,Y) (((X) >= (Y)) ? (X) : (Y))

// Prototypes

void expand_64bits(void *cpr, void *cpi, void *apr, void *api, mxComplexity complexity);
void expand_32bits(void *cpr, void *cpi, void *apr, void *api, mxComplexity complexity);
void expand_16bits(void *cpr, void *cpi, void *apr, void *api, mxComplexity complexity);
void expand_8bits(void *cpr, void *cpi, void *apr, void *api, mxComplexity complexity);

// global index variables

int i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12;
int k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12;
int A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12;
int B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12;

// Gateway routine

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    mxClassID classid;
    mxArray *Amx, *Bmx;
    mwSize Adims[12] = {1,1,1,1,1,1,1,1,1,1,1,1};
    mwSize Bdims[12] = {1,1,1,1,1,1,1,1,1,1,1,1};
    mwSize Cdims[12];
    int Andim, Bndim, Cndim;
    mxComplexity complexity;
    void *Apr, *Api, *Cpr, *Cpi;
    int i;
    
// Check for bit length definitions
    
    if( sizeof(bit16) != 2 )
        mexErrMsgTxt("Error using bsxarg: bit16 length is not 16-bits.");
    if( sizeof(bit32) != 4 )
        mexErrMsgTxt("Error using bsxarg: bit32 length is not 32-bits.");
    if( sizeof(bit64) != 8 )
        mexErrMsgTxt("Error using bsxarg: bit64 length is not 64-bits.");
    
// Check for proper inputs
    
    if( nrhs != 2 )
        mexErrMsgTxt("Error using bsxarg: Invalid number of inputs. Need two.");
    
    if( nlhs > 2 )
        mexErrMsgTxt("Error using bsxarg: Invalid number of outputs. Not more than two.");

// Set temporary operand pointers to the inputs.
        
    Amx = prhs[0];
    Bmx = prhs[1];
    
// Check for supported operand classes
    
    switch( mxGetClassID(Amx) )
    {
    case mxDOUBLE_CLASS:
    case mxINT64_CLASS:
    case mxUINT64_CLASS:
    case mxSINGLE_CLASS:
    case mxINT32_CLASS:
    case mxUINT32_CLASS:
    case mxCHAR_CLASS:
    case mxINT16_CLASS:
    case mxUINT16_CLASS:
    case mxLOGICAL_CLASS:
    case mxINT8_CLASS:
    case mxUINT8_CLASS:
        break;
    default:
        mexErrMsgTxt("Error using bsxarg: 1st argument class not supported");
    }
        
    switch( mxGetClassID(Bmx) )
    {
    case mxDOUBLE_CLASS:
    case mxINT64_CLASS:
    case mxUINT64_CLASS:
    case mxSINGLE_CLASS:
    case mxINT32_CLASS:
    case mxUINT32_CLASS:
    case mxCHAR_CLASS:
    case mxINT16_CLASS:
    case mxUINT16_CLASS:
    case mxLOGICAL_CLASS:
    case mxINT8_CLASS:
    case mxUINT8_CLASS:
        break;
    default:
        mexErrMsgTxt("Error using bsxarg: 2nd argument class not supported");
    }

// Generate the indexing multipliers    
    
    Andim = mxGetNumberOfDimensions(Amx);
    if( Andim > 12 )
        mexErrMsgTxt("Error using bsxarg: 1st array argument has too many dimensions.");
    memcpy(Adims, mxGetDimensions(Amx), Andim*sizeof(mwSize));
    
    A1 = 1;
    A2 = Adims[0];
    A3 = Adims[1]*A2;
    A4 = Adims[2]*A3;
    A5 = Adims[3]*A4;
    A6 = Adims[4]*A5;
    A7 = Adims[5]*A6;
    A8 = Adims[6]*A7;
    A9 = Adims[7]*A8;
    A10 = Adims[8]*A9;
    A11 = Adims[9]*A10;
    A12 = Adims[10]*A11;

    if( Adims[0] < 2 ) A1 = 0;
    if( Adims[1] < 2 ) A2 = 0;
    if( Adims[2] < 2 ) A3 = 0;
    if( Adims[3] < 2 ) A4 = 0;
    if( Adims[4] < 2 ) A5 = 0;
    if( Adims[5] < 2 ) A6 = 0;
    if( Adims[6] < 2 ) A7 = 0;
    if( Adims[7] < 2 ) A8 = 0;
    if( Adims[8] < 2 ) A9 = 0;
    if( Adims[9] < 2 ) A10 = 0;
    if( Adims[10] < 2 ) A11 = 0;
    if( Adims[11] < 2 ) A12 = 0;    
    
    Bndim = mxGetNumberOfDimensions(Bmx);
    if( Bndim > 12 )
        mexErrMsgTxt("Error using bsxarg: 2nd array argument has too many dimensions.");
    memcpy(Bdims, mxGetDimensions(Bmx), Bndim*sizeof(mwSize));

    B1 = 1;
    B2 = Bdims[0];
    B3 = Bdims[1]*B2;
    B4 = Bdims[2]*B3;
    B5 = Bdims[3]*B4;
    B6 = Bdims[4]*B5;
    B7 = Bdims[5]*B6;
    B8 = Bdims[6]*B7;
    B9 = Bdims[7]*B8;
    B10 = Bdims[8]*B9;
    B11 = Bdims[9]*B10;
    B12 = Bdims[10]*B11;

    if( Bdims[0] < 2 ) B1 = 0;
    if( Bdims[1] < 2 ) B2 = 0;
    if( Bdims[2] < 2 ) B3 = 0;
    if( Bdims[3] < 2 ) B4 = 0;
    if( Bdims[4] < 2 ) B5 = 0;
    if( Bdims[5] < 2 ) B6 = 0;
    if( Bdims[6] < 2 ) B7 = 0;
    if( Bdims[7] < 2 ) B8 = 0;
    if( Bdims[8] < 2 ) B9 = 0;
    if( Bdims[9] < 2 ) B10 = 0;
    if( Bdims[10] < 2 ) B11 = 0;
    if( Bdims[11] < 2 ) B12 = 0;
    
    k1 = MAX( Adims[0], Bdims[0] );
    k2 = MAX( Adims[1], Bdims[1] );
    k3 = MAX( Adims[2], Bdims[2] );
    k4 = MAX( Adims[3], Bdims[3] );
    k5 = MAX( Adims[4], Bdims[4] );
    k6 = MAX( Adims[5], Bdims[5] );
    k7 = MAX( Adims[6], Bdims[6] );
    k8 = MAX( Adims[7], Bdims[7] );
    k9 = MAX( Adims[8], Bdims[8] );
    k10 = MAX( Adims[9], Bdims[9] );
    k11 = MAX( Adims[10], Bdims[10] );
    k12 = MAX( Adims[11], Bdims[11] );
    
    for( i=0; i<12; i++ )
    {
        if( Adims[i] != Bdims[i] && Adims[i] != 1 && Bdims[i] != 1 )
            mexErrMsgTxt("Error using bsxarg: All non-singleton dimensions must match.");
        Cdims[i] = (Adims[i] && Bdims[i]) ? MAX( Adims[i], Bdims[i] ) : 0;
    }
    Cndim = MAX( Andim, Bndim );

// Fill in the values ---------------------------------------------------------

    classid = mxGetClassID(prhs[0]);
    complexity = mxIsComplex(prhs[0]) ? mxCOMPLEX : mxREAL;
    plhs[0] = mxCreateNumericArray(Cndim, Cdims, classid, complexity);
    
    if( mxGetNumberOfElements(plhs[0]) )
    {
        Cpr = mxGetData(plhs[0]);
        Cpi = mxGetImagData(plhs[0]);
        Apr = mxGetData(prhs[0]);
        Api = mxGetImagData(prhs[0]);
        
        switch( classid )
        {
        case mxDOUBLE_CLASS:
        case mxINT64_CLASS:
        case mxUINT64_CLASS:
            expand_64bits(Cpr, Cpi, Apr, Api, complexity);
            break;
            
        case mxSINGLE_CLASS:
        case mxINT32_CLASS:
        case mxUINT32_CLASS:
            expand_32bits(Cpr, Cpi, Apr, Api, complexity);
            break;
            
        case mxCHAR_CLASS:
        case mxINT16_CLASS:
        case mxUINT16_CLASS:
            expand_16bits(Cpr, Cpi, Apr, Api, complexity);
            break;
            
        case mxLOGICAL_CLASS:
        case mxINT8_CLASS:
        case mxUINT8_CLASS:
            expand_8bits(Cpr, Cpi, Apr, Api, complexity);
            break;
        }
    }

    if( nlhs < 2 ) return;
    
    classid = mxGetClassID(prhs[1]);
    complexity = mxIsComplex(prhs[1]) ? mxCOMPLEX : mxREAL;
    plhs[1] = mxCreateNumericArray(Cndim, Cdims, classid, complexity);
    
    if( mxGetNumberOfElements(plhs[1]) )
    {
        Cpr = mxGetData(plhs[1]);
        Cpi = mxGetImagData(plhs[1]);
        Apr = mxGetData(prhs[1]);
        Api = mxGetImagData(prhs[1]);
        
        A1 = B1;
        A2 = B2;
        A3 = B3;
        A4 = B4;
        A5 = B5;
        A6 = B6;
        A7 = B7;
        A8 = B8;
        A9 = B9;
        A10 = B10;
        A11 = B11;
        A12 = B12;
        
        switch( classid )
        {
        case mxDOUBLE_CLASS:
        case mxINT64_CLASS:
        case mxUINT64_CLASS:
            expand_64bits(Cpr, Cpi, Apr, Api, complexity);
            break;
            
        case mxSINGLE_CLASS:
        case mxINT32_CLASS:
        case mxUINT32_CLASS:
            expand_32bits(Cpr, Cpi, Apr, Api, complexity);
            break;
            
        case mxCHAR_CLASS:
        case mxINT16_CLASS:
        case mxUINT16_CLASS:
            expand_16bits(Cpr, Cpi, Apr, Api, complexity);
            break;
            
        case mxLOGICAL_CLASS:
        case mxINT8_CLASS:
        case mxUINT8_CLASS:
            expand_8bits(Cpr, Cpi, Apr, Api, complexity);
            break;
        }
    }

}

//-----------------------------------------------------------------------------
//
// Expansion function for 64-bit classes
//
//-----------------------------------------------------------------------------

void expand_64bits(void *cpr, void *cpi, void *apr, void *api, mxComplexity complexity)
{
    bit64 *Cpr = (bit64 *) cpr;
    bit64 *Cpi = (bit64 *) cpi;
    bit64 *Apr = (bit64 *) apr;
    bit64 *Api = (bit64 *) api;
    
    if( complexity == mxCOMPLEX )
    {
        for( i12=0; i12<k12; i12++ )
        for( i11=0; i11<k11; i11++ )
        for( i10=0; i10<k10; i10++ )
        for( i9=0; i9<k9; i9++ )
        for( i8=0; i8<k8; i8++ )
        for( i7=0; i7<k7; i7++ )
        for( i6=0; i6<k6; i6++ )
        for( i5=0; i5<k5; i5++ )
        for( i4=0; i4<k4; i4++ )
        for( i3=0; i3<k3; i3++ )
        for( i2=0; i2<k2; i2++ )
        for( i1=0; i1<k1; i1++ )
        {
            *Cpr++ = AX(Apr,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12);
            *Cpi++ = AX(Api,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12);
        }
    }
    else
    {
        for( i12=0; i12<k12; i12++ )
        for( i11=0; i11<k11; i11++ )
        for( i10=0; i10<k10; i10++ )
        for( i9=0; i9<k9; i9++ )
        for( i8=0; i8<k8; i8++ )
        for( i7=0; i7<k7; i7++ )
        for( i6=0; i6<k6; i6++ )
        for( i5=0; i5<k5; i5++ )
        for( i4=0; i4<k4; i4++ )
        for( i3=0; i3<k3; i3++ )
        for( i2=0; i2<k2; i2++ )
        for( i1=0; i1<k1; i1++ )
        {
            *Cpr++ = AX(Apr,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12);
        }
    }

}

//-----------------------------------------------------------------------------
//
// Expansion function for 32-bit classes
//
//-----------------------------------------------------------------------------

void expand_32bits(void *cpr, void *cpi, void *apr, void *api, mxComplexity complexity)
{
    bit32 *Cpr = (bit32 *) cpr;
    bit32 *Cpi = (bit32 *) cpi;
    bit32 *Apr = (bit32 *) apr;
    bit32 *Api = (bit32 *) api;
    
    if( complexity == mxCOMPLEX )
    {
        for( i12=0; i12<k12; i12++ )
        for( i11=0; i11<k11; i11++ )
        for( i10=0; i10<k10; i10++ )
        for( i9=0; i9<k9; i9++ )
        for( i8=0; i8<k8; i8++ )
        for( i7=0; i7<k7; i7++ )
        for( i6=0; i6<k6; i6++ )
        for( i5=0; i5<k5; i5++ )
        for( i4=0; i4<k4; i4++ )
        for( i3=0; i3<k3; i3++ )
        for( i2=0; i2<k2; i2++ )
        for( i1=0; i1<k1; i1++ )
        {
            *Cpr++ = AX(Apr,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12);
            *Cpi++ = AX(Api,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12);
        }
    }
    else
    {
        for( i12=0; i12<k12; i12++ )
        for( i11=0; i11<k11; i11++ )
        for( i10=0; i10<k10; i10++ )
        for( i9=0; i9<k9; i9++ )
        for( i8=0; i8<k8; i8++ )
        for( i7=0; i7<k7; i7++ )
        for( i6=0; i6<k6; i6++ )
        for( i5=0; i5<k5; i5++ )
        for( i4=0; i4<k4; i4++ )
        for( i3=0; i3<k3; i3++ )
        for( i2=0; i2<k2; i2++ )
        for( i1=0; i1<k1; i1++ )
        {
            *Cpr++ = AX(Apr,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12);
        }
    }

}

//-----------------------------------------------------------------------------
//
// Expansion function for 16-bit classes
//
//-----------------------------------------------------------------------------

void expand_16bits(void *cpr, void *cpi, void *apr, void *api, mxComplexity complexity)
{
    bit16 *Cpr = (bit16 *) cpr;
    bit16 *Cpi = (bit16 *) cpi;
    bit16 *Apr = (bit16 *) apr;
    bit16 *Api = (bit16 *) api;
    
    if( complexity == mxCOMPLEX )
    {
        for( i12=0; i12<k12; i12++ )
        for( i11=0; i11<k11; i11++ )
        for( i10=0; i10<k10; i10++ )
        for( i9=0; i9<k9; i9++ )
        for( i8=0; i8<k8; i8++ )
        for( i7=0; i7<k7; i7++ )
        for( i6=0; i6<k6; i6++ )
        for( i5=0; i5<k5; i5++ )
        for( i4=0; i4<k4; i4++ )
        for( i3=0; i3<k3; i3++ )
        for( i2=0; i2<k2; i2++ )
        for( i1=0; i1<k1; i1++ )
        {
            *Cpr++ = AX(Apr,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12);
            *Cpi++ = AX(Api,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12);
        }
    }
    else
    {
        for( i12=0; i12<k12; i12++ )
        for( i11=0; i11<k11; i11++ )
        for( i10=0; i10<k10; i10++ )
        for( i9=0; i9<k9; i9++ )
        for( i8=0; i8<k8; i8++ )
        for( i7=0; i7<k7; i7++ )
        for( i6=0; i6<k6; i6++ )
        for( i5=0; i5<k5; i5++ )
        for( i4=0; i4<k4; i4++ )
        for( i3=0; i3<k3; i3++ )
        for( i2=0; i2<k2; i2++ )
        for( i1=0; i1<k1; i1++ )
        {
            *Cpr++ = AX(Apr,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12);
        }
    }

}

//-----------------------------------------------------------------------------
//
// Expansion function for 8-bit classes
//
//-----------------------------------------------------------------------------

void expand_8bits(void *cpr, void *cpi, void *apr, void *api, mxComplexity complexity)
{
    bit8 *Cpr = (bit8 *) cpr;
    bit8 *Cpi = (bit8 *) cpi;
    bit8 *Apr = (bit8 *) apr;
    bit8 *Api = (bit8 *) api;
    
    if( complexity == mxCOMPLEX )
    {
        for( i12=0; i12<k12; i12++ )
        for( i11=0; i11<k11; i11++ )
        for( i10=0; i10<k10; i10++ )
        for( i9=0; i9<k9; i9++ )
        for( i8=0; i8<k8; i8++ )
        for( i7=0; i7<k7; i7++ )
        for( i6=0; i6<k6; i6++ )
        for( i5=0; i5<k5; i5++ )
        for( i4=0; i4<k4; i4++ )
        for( i3=0; i3<k3; i3++ )
        for( i2=0; i2<k2; i2++ )
        for( i1=0; i1<k1; i1++ )
        {
            *Cpr++ = AX(Apr,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12);
            *Cpi++ = AX(Api,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12);
        }
    }
    else
    {
        for( i12=0; i12<k12; i12++ )
        for( i11=0; i11<k11; i11++ )
        for( i10=0; i10<k10; i10++ )
        for( i9=0; i9<k9; i9++ )
        for( i8=0; i8<k8; i8++ )
        for( i7=0; i7<k7; i7++ )
        for( i6=0; i6<k6; i6++ )
        for( i5=0; i5<k5; i5++ )
        for( i4=0; i4<k4; i4++ )
        for( i3=0; i3<k3; i3++ )
        for( i2=0; i2<k2; i2++ )
        for( i1=0; i1<k1; i1++ )
        {
            *Cpr++ = AX(Apr,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12);
        }
    }

}
