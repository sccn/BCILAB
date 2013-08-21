%GMMB_HIST Create histS structure from data for PDF-value - density quantile mapping
%
%    histS = GMMB_HIST(data, type, bayesS)
%
%    data, type    are the training data used to create the bayesS.
%
%    This function creates ordered lists of training sample
%    PDF-values for PDF-value - density quantile mapping.
%
%    See gmmb_generatehist, gmmb_lhood2frac, gmmb_frac2lhood, gmmb_fracthresh
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
%    Jarmo Ilonen <jarmo.ilonen@lut.fi>
%    Joni Kamarainen <Joni.Kamarainen@lut.fi>
%
% Copyright:
%
%   Bayesian Classifier with Gaussian Mixture Model Pdf
%   functionality is Copyright (C) 2004 by Pekka Paalanen and
%   Joni-Kristian Kamarainen.
%
%   $Name:  $ $Revision: 1.2 $  $Date: 2005/04/14 10:33:34 $
%

function histS = gmmb_hist(data_, type_, bayesS);

K = size(bayesS,2);

histS = {};

for k = 1:K
	samples = data_(type_==k, :);
	dens = gmmb_pdf( samples, bayesS(k) );
	histS(k) = {sort(dens)};
end
