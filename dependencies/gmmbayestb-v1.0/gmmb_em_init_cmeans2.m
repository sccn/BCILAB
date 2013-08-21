% GMMB_EM_INIT_CMEANS2
%
% initS = gmmb_em_init_cmeans1(data, C)
%
% Create an initialization structure for EM,
% called from gmmb_em, see gmmb_em.
%
% C-means clustering means, cluster weight and covariance
%
% Author(s):
%    Pekka Paalanen <pekka.paalanen@lut.fi>
%
% Copyright:
%
%   Bayesian Classifier with Gaussian Mixture Model Pdf
%   functionality is Copyright (C) 2004 by Pekka Paalanen and
%   Joni-Kristian Kamarainen.
%
%   $Name:  $ $Revision: 1.2 $  $Date: 2004/11/02 09:00:18 $


function initS = gmmb_em_init_cmeans2(data, C);

N = size(data,1);	% number of points
D = size(data,2);	% dimensions

if C>1
	[lbl, mu] = gmmb_cmeans(data, C, 15);
	% initialization has random nature, results will vary
else
	lbl = ones(N, 1);
	mu = mean(data, 1);
end

% covariances initialization
sigma = zeros(D,D,C);
weight = ones(C,1);
for c = 1:C
	sigma(:,:,c) = gmmb_covfixer(cov(data(lbl==c, :)));
	weight(c) = sum(lbl==c) / N;
end

initS = struct(...
	'mu', mu.', ...
	'sigma', sigma, ...
	'weight', weight ...
	);

