//-------------------------------------------------------------------
#pragma hdrstop
//-------------------------------------------------------------------
//   C-MEX implementation of SUMSKIPNAN - this function is part of the NaN-toolbox. 
//
//   This program is free software; you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation; either version 2 of the License, or
//   (at your option) any later version.
//
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with this program; if not, see <http://www.gnu.org/licenses/>.
//
//
// sumskipnan: sums all non-NaN values
//
// Input:
// - array to sum
// - dimension to sum (1=columns; 2=rows; doesn't work for dim>2!!)
//
// Output:
// - sums
// - count of valid elements (optional)
// - sums of squares (optional)
// - sums of squares of squares (optional)
//
// Author:  Patrick Houweling (phouweling@yahoo.com)
// Version: 1.0
// Date:    17 september 2003
//
//    $Id: sumskipnan_mex.cpp 5590 2009-03-09 11:25:34Z schloegl $
//    Copyright (C) 2009 Alois Schloegl <a.schloegl@ieee.org>
//    This function is part of the NaN-toolbox
//    http://hci.tugraz.at/~schloegl/matlab/NaN/
//
//-------------------------------------------------------------------
//#include <stdlib>
#include <inttypes.h>
#include <math.h>
#include "mex.h"
//-------------------------------------------------------------------

int __sumskipnan2__(double *data, size_t Ni, size_t stride, double *s, size_t *No, char *flag_anyISNAN);
int __sumskipnan3__(double *data, size_t Ni, size_t stride, double *s, double *s2, size_t *No, char *flag_anyISNAN);
int __sumskipnan4__(double *data, size_t Ni, size_t stride, double *s, double *s2, double *s4, size_t *No, char *flag_anyISNAN);
int __sumskipnan2_single__(float *data, size_t Ni, size_t stride, double *s, size_t *No, char *flag_anyISNAN);
int __sumskipnan3_single__(float *data, size_t Ni, size_t stride, double *s, double *s2, size_t *No, char *flag_anyISNAN);
int __sumskipnan4_single__(float *data, size_t Ni, size_t stride, double *s, double *s2, double *s4, size_t *No, char *flag_anyISNAN);



void mexFunction(int POutputCount,  mxArray* POutput[], int PInputCount, const mxArray *PInputs[])
{
    	const int	*SZ;	    
    	double* 	LInput;
    	double* 	LInputI;
    	double* 	LOutputSum;
    	double* 	LOutputSumI;
    	double* 	LOutputCount;
    	double* 	LOutputSum2;
    	double* 	LOutputSum4;
    	double  	x, x2, x2i;
    	unsigned long   LCount, LCountI;
    	double  	LSum, LSum2, LSum4;

    	unsigned	DIM = 0; 
    	unsigned long	D1, D2, D3; 	// NN; 	//  	
    	unsigned    	ND, ND2;	// number of dimensions: input, output
    	unsigned long	ix1, ix2;	// index to input and output
    	unsigned    	j, k, l;	// running indices 
    	int 		*SZ2;		// size of output 	    
	char	 	flag_isNaN = 0;

	// check for proper number of input and output arguments
	if ((PInputCount <= 0) || (PInputCount > 3))
	        mexErrMsgTxt("SumSkipNan.MEX requires 1,2 or 3 arguments.");
	if (POutputCount > 4)
	        mexErrMsgTxt("SumSkipNan.MEX has 1 to 4 output arguments.");

	// get 1st argument
	if(mxIsDouble(PInputs[0]))
		LInput  = mxGetPr(PInputs[0]);
	else 	
		mexErrMsgTxt("First argument must be DOUBLE.");


	if(mxIsComplex(PInputs[0]))
		LInputI = mxGetPi(PInputs[0]);

    	// get 2nd argument
    	if  (PInputCount > 1){
 	       	switch (mxGetNumberOfElements(PInputs[1])) {
		case 0: x = 0.0; 		// accept empty element
			break;
		case 1: x = (mxIsNumeric(PInputs[1]) ? mxGetScalar(PInputs[1]) : -1.0); 
			break;
		default:x = -1.0;		// invalid 
		}
		if ((x < 0) || (x > 65535) || (x != floor(x))) 
			mexErrMsgTxt("Error SUMSKIPNAN.MEX: DIM-argument must be a positive integer scalar");

		DIM = (unsigned)floor(x);	
	}


	// get size 
    	ND = mxGetNumberOfDimensions(PInputs[0]);	
    	// NN = mxGetNumberOfElements(PInputs[0]);
    	SZ = mxGetDimensions(PInputs[0]);		

	// if DIM==0 (undefined), look for first dimension with more than 1 element. 
	for (k = 0; (DIM < 1) && (k < ND); k++) 
		if (SZ[k]>1) DIM = k+1;
	
	if (DIM < 1) DIM=1;		// in case DIM is still undefined 

	ND2 = (ND>DIM ? ND : DIM);	// number of dimensions of output 

	SZ2 = (int*)mxCalloc(ND2, sizeof(int)); // allocate memory for output size

	for (j=0; j<ND; j++)		// copy size of input;  
		SZ2[j] = SZ[j]; 	
	for (j=ND; j<ND2; j++)		// in case DIM > ND, add extra elements 1 
		SZ2[j] = 1; 	

    	for (j=0, D1=1; j<DIM-1; D1=D1*SZ2[j++]); 	// D1 is the number of elements between two elements along dimension  DIM  
	D2 = SZ2[DIM-1];		// D2 contains the size along dimension DIM 	
    	for (j=DIM, D3=1;  j<ND; D3=D3*SZ2[j++]); 	// D3 is the number of blocks containing D1*D2 elements 

	SZ2[DIM-1] = 1;		// size of output is same as size of input but SZ(DIM)=1;

	    // create outputs
	#define TYP mxDOUBLE_CLASS

	if(mxIsComplex(PInputs[0]))
	{	POutput[0] = mxCreateNumericArray(ND2, SZ2, TYP, mxCOMPLEX);
		LOutputSum = mxGetPr(POutput[0]);
		LOutputSumI= mxGetPi(POutput[0]);
    	}
	else 
	{	POutput[0] = mxCreateNumericArray(ND2, SZ2, TYP, mxREAL);
		LOutputSum = mxGetPr(POutput[0]);
    	}

    	if (POutputCount >= 2){
		POutput[1] = mxCreateNumericArray(ND2, SZ2, TYP, mxREAL);
        	LOutputCount = mxGetPr(POutput[1]);
    	}
    	if (POutputCount >= 3){
		POutput[2] = mxCreateNumericArray(ND2, SZ2, TYP, mxREAL);
        	LOutputSum2  = mxGetPr(POutput[2]);
    	}
    	if (POutputCount >= 4){
		POutput[3] = mxCreateNumericArray(ND2, SZ2, TYP, mxREAL);
        	LOutputSum4  = mxGetPr(POutput[3]);
    	}

	mxFree(SZ2);

	if ((POutputCount <3)	&& !mxIsComplex(PInputs[0]))
	{
		// OUTER LOOP: along dimensions > DIM
		for (l = 0; l<D3; l++) 	
	    	{
			ix2 =   l*D1;	// index for output 
			ix1 = ix2*D2;	// index for input 

			// Inner LOOP: along dimensions < DIM
			for (k = 0; k<D1; k++, ix1++, ix2++) 	
			{
				int flag=0;
				size_t count;
				__sumskipnan2__(LInput+ix1, D2, D1, LOutputSum+ix2, &count, &flag_isNaN);
				if (POutputCount > 1)
					LOutputCount[ix2] = double(count);
               		}
               	}		
	}

	else if ((POutputCount == 3) && !mxIsComplex(PInputs[0]))
	{
		// OUTER LOOP: along dimensions > DIM
		for (l = 0; l<D3; l++) 	
    		{
			ix2 =   l*D1;	// index for output 
			ix1 = ix2*D2;	// index for input 

			// Inner LOOP: along dimensions < DIM
			for (k = 0; k<D1; k++, ix1++, ix2++) 	
			{
				size_t count;
				__sumskipnan3__(LInput+ix1, D2, D1, LOutputSum+ix2, LOutputSum2+ix2, &count, &flag_isNaN);
				LOutputCount[ix2]=double(count);
 	       		}
               	}		
	}
	else if ((POutputCount == 4) && !mxIsComplex(PInputs[0]))
	{
		// OUTER LOOP: along dimensions > DIM
		for (l = 0; l<D3; l++) 	
    		{
			ix2 =   l*D1;	// index for output 
			ix1 = ix2*D2;	// index for input 

			// Inner LOOP: along dimensions < DIM
			for (k = 0; k<D1; k++, ix1++, ix2++) 	
			{
				size_t count;
				__sumskipnan4__(LInput+ix1, D2, D1, LOutputSum+ix2, LOutputSum2+ix2, LOutputSum4+ix2, &count, &flag_isNaN);
				LOutputCount[ix2]=double(count);
               		}
               	}		
	}
	else if ((POutputCount < 3) && mxIsComplex(PInputs[0]))
	{
		// OUTER LOOP: along dimensions > DIM
		for (l = 0; l<D3; l++) 	
    		{
			ix2 =   l*D1;	// index for output 
			ix1 = ix2*D2;	// index for input 

			// Inner LOOP: along dimensions < DIM
			for (k = 0; k<D1; k++, ix1++, ix2++) 	
			{
				int flag=0;
				size_t count,countI;
				__sumskipnan2__(LInput+ix1, D2, D1, LOutputSum+ix2, &count, &flag_isNaN);

				__sumskipnan2__(LInputI+ix1, D2, D1, LOutputSumI+ix2, &countI, &flag_isNaN);

				if (count != countI)
		            		mexErrMsgTxt("Number of NaNs is different for REAL and IMAG part");

				if (POutputCount > 1)
					LOutputCount[ix2]=double(count);
			}
		}
	}
	else if ((POutputCount == 3) && mxIsComplex(PInputs[0]))
	{
		// OUTER LOOP: along dimensions > DIM
		for (l = 0; l<D3; l++) 	
    		{
			ix2 =   l*D1;	// index for output 
			ix1 = ix2*D2;	// index for input 

			// Inner LOOP: along dimensions < DIM
			for (k = 0; k<D1; k++, ix1++, ix2++) 	
			{
				size_t count,countI;
				double ssq1,ssq2;
				__sumskipnan3__(LInput+ix1, D2, D1, LOutputSum+ix2, &ssq1, &count, &flag_isNaN);

				__sumskipnan3__(LInputI+ix1, D2, D1, LOutputSumI+ix2, &ssq2, &countI, &flag_isNaN);

				if (count != countI)
		            		mexErrMsgTxt("Number of NaNs is different for REAL and IMAG part");

				LOutputCount[ix2]= double(count);
				LOutputSum2[ix2] = ssq1+ssq2;
			}
		}
	}
	else if ((POutputCount == 4) && mxIsComplex(PInputs[0]))
	{
		// OUTER LOOP: along dimensions > DIM
		for (l = 0; l<D3; l++) 	
    		{
			ix2 =   l*D1;	// index for output 
			ix1 = ix2*D2;	// index for input 

			// Inner LOOP: along dimensions < DIM
			for (k = 0; k<D1; k++, ix1++, ix2++) 	
			{
		        	LCount = 0;
				LSum   = 0.0;
				LSum2  = 0.0;
				LSum4  = 0.0;
	        	    		
				// LOOP  along dimension DIM
		    		for (j=0; j<D2; j++) 	
				{
					x = LInput[ix1 + j*D1];
        	        		if (!mxIsNaN(x))
					{
						LCount++; 
						LSum += x; 
						x2 = x*x;
						LSum2 += x2; 
					}
					x = LInputI[ix1 + j*D1];
       	        			if (!mxIsNaN(x))
					{
						LCountI++; 
						LSum  += x; 
						x2i = x*x;
						LSum2 += x2i; 
						x2 += x2i; 
						LSum4 += x2*x2; 
					}
				}	
				if (LCount != LCountI)
		            		mexErrMsgTxt("Number of NaNs is different for REAL and IMAG part");

				LOutputSum[ix2] = LSum;
               			LOutputCount[ix2] = (double)LCount;
       	        		LOutputSum2[ix2] = LSum2;
               			LOutputSum4[ix2] = LSum4;

			}
		}
	}

    	if  (flag_isNaN && (PInputCount > 2)) {
    		// set FLAG_NANS_OCCURED 
    		switch (mxGetClassID(PInputs[2])) {
    		case mxDOUBLE_CLASS:
    			*(double*)mxGetData(PInputs[2]) = 1.0;
    			break; 
    		case mxSINGLE_CLASS:
    			*(float*)mxGetData(PInputs[2]) = 1.0;
    			break; 
    		case mxLOGICAL_CLASS:
    		case mxCHAR_CLASS:
    		case mxINT8_CLASS:
    		case mxUINT8_CLASS:
    			*(uint8_t*)mxGetData(PInputs[2]) = 1;
    			break; 
    		case mxINT16_CLASS:
    		case mxUINT16_CLASS:
    			*(uint16_t*)mxGetData(PInputs[2]) = 1;
    			break; 
    		case mxINT32_CLASS:
    		case mxUINT32_CLASS:
    			*(uint32_t*)mxGetData(PInputs[2])= 1;
    			break; 
    		case mxINT64_CLASS:
    		case mxUINT64_CLASS:
    			*(uint64_t*)mxGetData(PInputs[2]) = 1;
    			break; 
    		case mxFUNCTION_CLASS:
    		case mxUNKNOWN_CLASS:
    		case mxCELL_CLASS:
    		case mxSTRUCT_CLASS:
    			;
		}
	}
}


int __sumskipnan2__(double *data, size_t Ni, size_t stride, double *s, size_t *No, char *flag_anyISNAN)
{
	register double sum=0; 
	register size_t count=0; 
	register char   flag=0; 
	// LOOP  along dimension DIM
	
	for (size_t j=0; j<Ni; j++, data += stride)
	{
		register double x = *data;
        	if (x==x)
		{
			count++; 
			sum += x; 
		}
		else 
			flag = 1; 
	}

	if (flag && (flag_anyISNAN != NULL)) *flag_anyISNAN = 1; 
	*s  = sum;
        *No = count;

}

int __sumskipnan3__(double *data, size_t Ni, size_t stride, double *s, double *s2, size_t *No, char *flag_anyISNAN)
{
	register double sum=0; 
	register double msq=0; 
	register size_t count=0; 
	register char   flag=0; 
	// LOOP  along dimension DIM
	
	for (size_t j=0; j<Ni; j++, data += stride)
	{
		register double x = *data;
        	if (x==x)
		{
			count++; 
			sum += x; 
			msq += x*x; 
		}
		else 
			flag = 1; 
	}

	if (flag && (flag_anyISNAN != NULL)) *flag_anyISNAN = 1; 
	*s  = sum;
	*s2 = msq; 
        *No = count;
}

int __sumskipnan4__(double *data, size_t Ni, size_t stride, double *s, double *s2, double *s4, size_t *No, char *flag_anyISNAN)
{
	register double _s0=0; 
	register double _s2=0; 
	register double _s4=0; 
	register size_t count=0; 
	register char   flag=0; 
	// LOOP  along dimension DIM
	
	for (size_t j=0; j<Ni; j++, data += stride)
	{
		register double x = *data;
        	if (x==x)
		{
			count++; 
			_s0 += x; 
			x =x*x;
			_s2 += x; 
			_s4 += x*x; 
		}
		else 
			flag = 1; 
	}

	if (flag && (flag_anyISNAN != NULL)) *flag_anyISNAN = 1; 
	*s  = _s0;
	*s2 = _s2; 
	*s4 = _s4; 
        *No = count;
}

int __sumskipnan2_single__(float *data, size_t Ni, size_t stride, double *s, size_t *No, char *flag_anyISNAN)
{
	register double sum=0; 
	register size_t count=0; 
	register char   flag=0; 
	// LOOP  along dimension DIM
	
	for (size_t j=0; j<Ni; j++, data += stride)
	{
		register double x = *data;
        	if (x==x)
		{
			count++; 
			sum += x; 
		}
		else 
			flag = 1; 
	}

	if (flag && (flag_anyISNAN != NULL)) *flag_anyISNAN = 1; 
	*s  = sum;
        *No = count;

}

int __sumskipnan3_single__(float *data, size_t Ni, size_t stride, double *s, double *s2, size_t *No, char *flag_anyISNAN)
{
	register double sum=0; 
	register double msq=0; 
	register size_t count=0; 
	register char   flag=0; 
	// LOOP  along dimension DIM
	
	for (size_t j=0; j<Ni; j++, data += stride)
	{
		register double x = *data;
        	if (x==x)
		{
			count++; 
			sum += x; 
			msq += x*x; 
		}
		else 
			flag = 1; 
	}

	if (flag && (flag_anyISNAN != NULL)) *flag_anyISNAN = 1; 
	*s  = sum;
	*s2 = msq; 
        *No = count;
}

int __sumskipnan4_single__(float *data, size_t Ni, size_t stride, double *s, double *s2, double *s4, size_t *No, char *flag_anyISNAN)
{
	register double _s0=0; 
	register double _s2=0; 
	register double _s4=0; 
	register size_t count=0; 
	register char   flag=0; 
	// LOOP  along dimension DIM
	
	for (size_t j=0; j<Ni; j++, data += stride)
	{
		register double x = *data;
        	if (x==x)
		{
			count++; 
			_s0 += x; 
			x =x*x;
			_s2 += x; 
			_s4 += x*x; 
		}
		else 
			flag = 1; 
	}

	if (flag && (flag_anyISNAN != NULL)) *flag_anyISNAN = 1; 
	*s  = _s0;
	*s2 = _s2; 
	*s4 = _s4; 
        *No = count;
}


#ifdef experimental
	/* x86 assembler code, currently broken, 
		the main advantage would be the support of the extended accuracy  
	*/
	__asm__ ("	movl $0, %eax;\n" 
		 "	fldz;\n" 
		 "	movl _D2, %ebx;\n" 
		 "	movl _ptr, %edx;\n" 
		 "loop3: \n"
		 "	fld (%%edx);\n" 
		 "	fcom (ST0);\n"
		 "	jne end_isnan;"	 
		 "	fadd (ST0)\n"
		 "	inc %eax\n"
		 "	\n"
		 "end_isnan: \n"
		 "	fdecstp;\n" 
		 "	add %edx,_stride\n"
		 "	loop loop3 %ebx;\n" 
		 "	fstp _LSum;\n"
		 "	movl %eax, _LCount;\n"
		 "	nop\n" 
		);
		// #__asm__ ("fldz " : : : "%eax");
#endif 


