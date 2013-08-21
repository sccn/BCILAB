%GMMB_GEM    - Greedy EM estimated GMM parameters
% Produces a bayesS struct without 'apriories'
% This is just a wrapper for the Vlassis Greedy EM algorithm implementation.
%
% estimate = GMMB_GEM(data[, parameters])
% [estimate,stats] = GMMB_GEM(...)
%
% Parameters (default):
%   verbose	print some progress numbers (false)
%   animate	plot data and ellipses during algorithm evaluation (false)
%   Cmax	the maximum number of GMM components
%   ncand	number of candidate locations for each new component (10)
% At least Cmax should be set explicitly.
% example:
%    estS = gmmb_gem(data, 'Cmax', 10, 'animate', true);
%
% References:
%   [1] Vlassis, N., Likas, A., A Greedy EM Algorithm for Gaussian Mixture
%   Learning, Neural Processing Letters 15, Kluwer Academic Publishers, 2002.
%   http://carol.wins.uva.nl/~vlassis/research/learning/index_en.html
%
% Author(s):
%    Pekka Paalanen <pekka.paalanen@lut.fi>
%
% Copyright:
%
%   Bayesian Classifier with Gaussian Mixture Model Pdf
%   functionality is Copyright (C) 2003 by Pekka Paalanen and
%   Joni-Kristian Kamarainen.
%
%   $Name:  $ $Revision: 1.1 $  $Date: 2004/11/02 08:32:22 $
%
%
% Logging
%   parameters
%
%      logging   What kind of logging to do:
%        0 - no logging
%        1 - normal logging
%        2 - extra logging: store all intermediate mixtures
%      If the 'stats' output parameter is defined, then 'logging'
%      defaults to 1, otherwise it is forced to 0.
%
%  the 'stats' struct:
%      iterations: EM iteration count
%      loglikes:   iterations long vector of the log-likelihood
%                NOTE: mean(log(p)), not sum(log(p)) as it should(?)
%    extra logging: (not supported yet)
%      initialmix: parameters for the initial mixture
%      mixtures:   parameters for all intermediate mixtures
%

function [estimate, varargout] = gmmb_gem(data, varargin);

[N, D] = size(data);	% number of points (n), dimensions (d)

% defaults
conf = struct(...
	'verbose', 0, ...
	'animate', 0, ...
	'Cmax', ceil(min(50, N/(D*D)/3)), ...
	'ncand', 10, ...
	'logging', 0 ...
	);

if nargout>1
	conf.logging = 1;
	varargout{1} = [];
end

conf = getargs(conf, varargin);

C = conf.Cmax;

N_limit = (D+D*(D+1)/2+1)*3;
if N < N_limit
	warning_wrap('gmmb_gem:data_amount', ...
	   ['Training data may be insufficient. ' ...
	    'Have: ' num2str(N) ', recommended: >' num2str(N_limit) ...
	    ' points.']);
end


if nargout<2
	conf.logging=0;
end

[W, M, R, stats] = gmmbvl_em(data, C, conf.ncand, conf.animate, ...
                             conf.verbose, conf.logging);

Cfinal = size(R,1);
sigma = zeros(D, D, Cfinal);
for c = 1:Cfinal
	Rk = reshape(R(c,:),D,D);
	sigma(:,:,c) = Rk' * Rk;
end

estimate = struct('mu', M.',...
	'sigma', sigma,...
	'weight', W);

if(conf.logging>0)
	varargout{1} = stats;
end
