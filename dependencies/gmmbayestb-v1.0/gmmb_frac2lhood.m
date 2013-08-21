%GMMB_FRAC2LHOOD   Map density quantiles to PDF threshold values
%
%    lhood = GMMB_FRAC2LHOOD(histS, f)
%
%    histS   K-element histS cell array created by gmmb_hist or
%            gmmb_generatehist.
%    lhood   N x K array of likelihood values.
%    f       N x K array of density quantile values
%
%    This function finds the likelihood threshold value corresponding to
%    each density quantile value.
%    For each column k in 1..K, the likelihood value is found from histS{k},
%    so that each column may represent a different distribution.
%
%    See gmmb_hist, gmmb_generatehist, gmmb_lhood2frac, gmmb_fracthresh
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

function lhood = gmmb_frac2lhood(histS, f);

if any(f(:)>1 | f(:)<0)
	error('Density quantile values must be in the range [0,1].');
end

lhood = zeros(size(f));

K = size(f, 2);

for k = 1:K
	v = shiftdim(histS{k});
	len_v = length(v);
	nf = (1-f(:,k)) .* (len_v-1) +1;
	i = floor(nf);

	v(len_v+1) = v(len_v);
	lhood(:, k) = v(i) + ( nf-i ) .* ( v(i+1) - v(i) );
end
