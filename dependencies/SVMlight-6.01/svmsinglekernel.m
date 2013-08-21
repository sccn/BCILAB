function model = svmsinglekernel(x,y,model_parm)
% SVMSINGLEKERNEL - compute the kernel value for the given pair of vectors
%
%   SVM-lite is compiled as a DLL and is executed through the
%   MEX interface.  There is a penalty incurred through conversion
%   of the vectors from MATLAB's column-major format for arrays to
%   the C standard row-major format.  In practice, this penalty is
%   inconsequential.  The overall performance is nearly identical to
%   that of SVM-Lite running as a stand-alone program.
    clear fun mexsvmclassify;
    clear fun mexsvmlearn;
    clear fun mexkernel;
    clear fun mexsinglekernel;
    
    result = mexsvmsinglekernel(x,y,model_parm);
    
    clear fun mexsinglekernel;
    