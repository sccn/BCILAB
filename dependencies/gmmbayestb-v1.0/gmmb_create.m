%GMMB_CREATE - Construct new Bayesian classifier with Gaussian mixture model pdf
%
%     S = GMMB_CREATE(data, class, method [, parameters]) Generates a
%     Bayesian classifier for one or several classes having GMM
%     distribution with estimated mean values, variances and
%     apriories. Classifier is returned in bayesS struct S.
%
%     method can be 'EM', 'FJ' or 'GEM'.
%       EM and FJ can work with complex numbers.
%
%     See also readme.txt, GMMB_HIST, GMMB_GENERATEHIST.
%
%     Parameters are delegated directly to underlaying GMM estimation
%     function (gmmb_em, gmmb_fj, gmmb_gem). See also them.
%
% Examples:
%   bayesS = gmmb_create(data, class, 'EM', 'components', 5, 'thr', 1e-8);
%   bayesS = gmmb_create(data, class, 'FJ', 'Cmax', 50, 'thr', 1e-9);
%   bayesS = gmmb_create(data, class, 'GEM', 'Cmax', 10, 'verbose', true);
%
%   The bayesS struct is documented in readme.txt.
%
% References:
%   [1] Duda, R.O., Hart, P.E, Stork, D.G, Pattern Classification,
%   2nd ed., John Wiley & Sons, Inc., 2001.
%
% Author(s):
%    Joni Kamarainen <Joni.Kamarainen@lut.fi>
%    Pekka Paalanen <pekka.paalanen@lut.fi>
%
% Copyright:
%
%   Bayesian Classifier with Gaussian Mixture Model Pdf is
%   Copyright (C) 2003, 2004 by Pekka Paalanen and Joni-Kristian
%   Kamarainen.
%
%   $Name:  $ $Revision: 1.1 $  $Date: 2004/11/02 08:32:22 $
%

function [bayesS, varargout] = gmmb_create(data, cl, method, varargin);

K = max(cl);

mu ={};
sigma = {};
weight = {};
prior = {};
stats = {};

for k = 1:K
	cvals = data(cl == k, :);
	N = size(cvals,1);	% points
	
	switch method
	case 'EM'
		[estim stat] = gmmb_em(cvals, varargin{:});
	case 'FJ'
		[estim stat] = gmmb_fj(cvals, varargin{:});
	case 'GEM'
		[estim stat] = gmmb_gem(cvals, varargin{:});
	otherwise
		error('Unknown method');
	end

	mu{k} = estim.mu;
	sigma{k} = estim.sigma;
	weight{k} = estim.weight;
	prior{k} = N/size(data,1);
	stats{k} = stat;
end

bayesS = struct('mu', mu,...
		'sigma', sigma,...
		'apriories', prior,...
		'weight', weight);

if nargout > 1
	varargout{1} = stats;
end

