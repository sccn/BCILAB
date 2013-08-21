%GMMB_FRACTHRESH    Threshold PDF values according to density quantile.
%
%     MASK = GMMB_FRACTHRESH(pdfmat, histS, thr)
%
%     pdfmat = N x K matrix of PDF values at N points
%              in K different PDFs (the output of gmmb_pdf)
%     histS = the histS structure (1 x K cell array) created with
%             the bayesS structure that was used to compute PDFs.
%     thr = scalar in the range [0, 1], the density quantile
%
%     MASK = N x K logical matrix
%
%     See also GMMB_PDF, GMMB_HIST, GMMB_GENERATEHIST.
%
%  The recommended way to create histS is with gmmb_generatehist.
%  The output is a logical N x K matrix that tells whether point N
%  is an outlier in distribution K, i.e., it does not belong to the
%  thr-quantile of distribution K.
%
% References:
%   [1] Paalanen, P., Kamarainen, J.-K., Ilonen, J., Kälviäinen, H.,
%    Feature Representation and Discrimination Based on Gaussian Mixture Model
%    Probability Densities - Practices and Algorithms, Research Report 95,
%    Lappeenranta University of Technology, Department of Information
%    Technology, 2005.
%
% Author(s):
%    Pekka Paalanen <pekka.paalanen@lut.fi>
%
% Copyright:
%
%   Bayesian Classifier with Gaussian Mixture Model Pdf
%   functionality is Copyright (C) 2003, 2004 by Pekka Paalanen and
%   Joni-Kristian Kamarainen.
%
%   $Name:  $ $Revision: 1.2 $  $Date: 2005/04/14 10:33:34 $
%

function mask = gmmb_fracthresh(pdfmat, histS, thr);

N = size(pdfmat, 1);
K = size(pdfmat, 2);

thresh = gmmb_frac2lhood(histS, thr*ones(1,K));
mask = (pdfmat < repmat(thresh, N, 1));
