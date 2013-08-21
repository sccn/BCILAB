%GMMB_NORMALIZE   Normalize a matrix so that row sums are one.
%
%     p_out = GMMB_NORMALIZE(p_in)
%
%     p_in = N x K matrix
%     p_out = N x K matrix
%
%  If an unnormalized row sum would be zero,
%  the row is left untouched.
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

function p_out = gmmb_normalize(p_in);

K = size(p_in, 2);
divi = sum(p_in, 2);
divi(divi==0) = 1;

p_out = p_in ./  repmat(divi, 1, K);
