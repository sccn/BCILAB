%GMMB_DECIDE   Make decisions; choose index of the max value.
%
%     labels = GMMB_DECIDE(p_in)
%
%     p_in = N x K matrix
%     labels = N x 1 matrix of integers
%
%  The labels will be index of the maximum value on each row,
%  except if the max value is zero, index will also be zero.
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

function labels = gmmb_decide(p_in);

[vals, labels] = max(p_in, [], 2);
labels(vals==0) = 0;
