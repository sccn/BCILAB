#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /**********************\
     MATLAB
     M = getmat(m)
    \**********************/

    const char *classname;
    const mxArray *m, *m_mat, *m_dmmask;
    
    classname   = mxGetClassName(prhs[0]);
    
    if(strcmp(classname, "ssmat") == 0) m = prhs[0];
    else if(strcmp(classname, "ssdist") == 0) m = mxGetField(prhs[0], 0, "ssmat");
    else if(strcmp(classname, "ssfunc") == 0) m = mxGetField(prhs[0], 0, "ssmat");
    else return;
    
    m_mat       = mxGetField(m, 0, "mat");
    m_dmmask    = mxGetField(m, 0, "dmmask");
    
    if(mxIsEmpty(m_dmmask)) plhs[0] = mxDuplicateArray(m_mat);
    else {
        int i, j, t;
        int dims[3], m1, m2, d, n;
        const double *mat, *dvec;
        double *M;
        const bool *dmmask;
        
        const mxArray *m_dvec = mxGetField(m, 0, "dvec");

        m1  = mxGetM(m_mat);
        m2  = mxGetN(m_mat);
        d   = mxGetM(m_dvec);
        n   = mxGetN(m_dvec);
        
        /*** Get input matrices ***/
        mat     = mxGetPr(m_mat);
        dmmask  = (bool *)mxGetLogicals(m_dmmask);
        dvec    = mxGetPr(m_dvec);
        
        /*** Create output matrices ***/
        dims[0] = m1; dims[1] = m2; dims[2] = n;
        M       = mxGetPr(plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL));
        
        /*** Set output ***/
        j       = 0;
        for(i=0; i<m1*m2; i++) {
            if(dmmask[i]) {
                for(t=0; t<n; t++) M[i+t*m1*m2] = dvec[j+t*d];
                j++;
            }
            else for(t=0; t<n; t++) M[i+t*m1*m2] = mat[i];
        }
    }
}
