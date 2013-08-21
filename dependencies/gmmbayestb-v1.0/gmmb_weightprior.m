%GMMB_WEIGHTPRIOR   Multiply PDF values with constant priors
%
%     P = GMMB_WEIGHTPRIOR(pdfmat, bayesS)
%
%     pdfmat = N x K matrix of PDF values at N points
%              in K different PDFs (the output of gmmb_pdf)
%     bayesS = the bayesS struct used to compute pdfmat,
%              used fields: apriories
%     P = N x K matrix of weighted PDF values.
%
%     See also GMMB_PDF.
%
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
%   $Name:  $ $Revision: 1.1 $  $Date: 2004/11/02 08:32:22 $
%

function p = gmmb_weightprior(pdfmat, bayesS);

N = size(pdfmat, 1);
priors = [bayesS.apriories];

p = pdfmat .* repmat(priors, N, 1);
