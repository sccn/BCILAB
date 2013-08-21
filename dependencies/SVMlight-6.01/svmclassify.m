function [ err, predictions ] = svmclassify(x,y,model)
% SVMCLASSIFY uses SVM-lite to classify the points stored in the
%   data matrix (x), and compares them to the labels in (y).  The
%   SVM model is a structure containing the support vectors and other
%   necessary parameters that were used to train the SVM.  This is 
%   typically created from a call to SVMLEARN.
%
%   The output includes the error rate (count where sgn(pred(x)) != sgn(y)
%
%   SVM-lite library is compiled as a DLL and is executed through the
%   MEX interface.  There is a penalty incurred through conversion
%   of the vectors from MATLAB's column-major format for arrays to
%   the C standard row-major format.  In practice, this penalty is
%   inconsequential.  The overall performance is nearly identical to
%   that of SVM-Lite running as a stand-alone program.
    clear fun mexsvmclassify;
    clear fun mexsvmlearn;

    [ err, predictions ] = mexsvmclassify(x,y,model);
    
    clear fun mexsvmclassify;
    