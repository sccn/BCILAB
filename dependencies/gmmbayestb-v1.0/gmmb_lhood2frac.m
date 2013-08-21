%GMMB_LHOOD2FRAC   Map likelihood values to density quantiles
%
%    f = GMMB_LHOOD2FRAC(histS, lhood)
%
%    histS   K-element cell array created by gmmb_hist or gmmb_generatehist
%    lhood   N x K array of likelihood values.
%    f       N x K array of density quantile values
%
%    This function finds the corresponding density quantile value of each
%    likelihood value in the "lhood" array.
%    For each column k in 1..K, the density quantile is found from histS{k},
%    so that each column may represent a different distribution.
%
%    See gmmb_hist, gmmb_generatehist, gmmb_frac2lhood, gmmb_fracthresh
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

function fracmat = gmmb_lhood2frac(histS, lhood);

fracmat = zeros(size(lhood));

K = size(lhood, 2);

for k = 1:K
	fracmat(:, k) = my_l2f_bin(histS{k}, lhood(:,k));
end


% binary search, linear interpolating version
function f = my_l2f_bin(v, lhoods);

f = zeros(size(lhoods));
len_v = length(v);

for i = 1:length(lhoods)
	x = lhoods(i);
	a = 1;
	b = len_v+1;
	
	if x < v(a)
		f(i) = 1;
		continue;
	end
	
	if x >= v(b-1)
		f(i) = 0;
		continue;
	end
	
	while (b-a) > 1
		m = floor((a+b)/2);
		if( v(m) > x)
			b = m;
		else
			a = m;
		end
	end
	
	h = v(b) - v(a);
	if h ~= 0
		sf = (x - v(a)) / h;
	else
		sf = 0.5;
	end
	
	f(i) = 1 - (a+sf-1) / (len_v-1);
end


