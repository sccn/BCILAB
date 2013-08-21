%GMMB_EM     - EM estimated GMM parameters
%
% estS = gmmb_em(data)
% estS = gmmb_em(data, <params>...)
% [estS, stats] = gmmb_em(...)
%
% This version works with complex numbers too.
%
% data = N x D matrix
% params can be a list of 'name', value -pairs.
% stats is a matrix, row (cov fixes, loops, final log-likelihood)
%
% Parameters (default value):
%
%   maxloops	maximum number of loops allowed (100)
%   thr		convergence threshold; relative log-likelihood change (1e-6)
%		set negative to use only maxloops condition
%   components	number of components in GMM (3)
%   verbose	print progress messages (false)
%   init        name of the initialization method to use ('fcm1')
%
% Init methods:
%
%  fcm1         Fuzzy C-means clustering, requires the Fuzzy Logic Toolbox.
%               This is the original init method from GMMBayes Toolbox v0.1
%  cmeans1      C-means clustering for means, uniform weigths and covariances
%  cmeans2      C-means clustering for means, weigths and covariances
%
% Example:
%   estS = gmmb_em(data, 'init', 'fcm1', 'components', 5, 'thr', 1e-8)
%
% References:
%   [1] Duda, R.O., Hart, P.E, Stork, D.G, Pattern Classification,
%   2nd ed., John Wiley & Sons, Inc., 2001.
%   [2] Bilmes, J.A., A Gentle Tutorial of the EM Algorithm and its
%    Application to Parameter Estimation for Gaussian Mixture and Hidden
%    Markov Models
%   International Computer Science Institute, 1998
%
% Author(s):
%    Joni Kamarainen <Joni.Kamarainen@lut.fi>
%    Pekka Paalanen <pekka.paalanen@lut.fi>
%
% Copyright:
%
%   Bayesian Classifier with Gaussian Mixture Model Pdf
%   functionality is Copyright (C) 2003 by Pekka Paalanen and
%   Joni-Kristian Kamarainen.
%
%   $Name:  $ $Revision: 1.2 $  $Date: 2004/11/02 09:00:18 $
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
%      covfixer2:  iterations-by-C matrix of gmmb_covfixer fix round counts
%      loglikes:   iterations long vector of the log-likelihood
%    extra logging:
%      initialmix: parameters for the initial mixture
%      mixtures:   parameters for all intermediate mixtures
%


function [estimate, varargout] = gmmb_em(data, varargin);

% default parameters
conf = struct(...
	'maxloops', 100, ...
	'thr', 1e-6, ...
	'verbose', 0, ...
	'components', 3, ...
	'logging', 0, ...
	'init', 'fcm1' ...
	);

if nargout>1
	conf.logging = 1;
	varargout{1} = [];
end

conf = getargs(conf, varargin);

if nargout<2
	conf.logging=0;
end

% for logging
log_covfixer2 = {};
log_loglikes = {};
log_initialmix = {};
log_mixtures = {};


% --- initialization ---

N = size(data,1);	% number of points
D = size(data,2);	% dimensions
C = conf.components;

% the number of free parameters in a Gaussian
if isreal(data)
	Nparc = D+D*(D+1)/2;
else
	Nparc = 2*D + D*D;
end
N_limit = (Nparc+1)*3*C;
if N < N_limit
	warning_wrap('gmmb_em:data_amount', ...
	   ['Training data may be insufficient for selected ' ...
	    'number of components. ' ...
	    'Have: ' num2str(N) ', recommended: >' num2str(N_limit) ...
	    ' points.']);
end

switch lower(conf.init)
	case 'fcm1'
		initS = gmmb_em_init_fcm1(data, C, conf.verbose);
	case 'cmeans1'
		initS = gmmb_em_init_cmeans1(data, C);
	case 'cmeans2'
		initS = gmmb_em_init_cmeans2(data, C);
	otherwise
		error(['Unknown initializer method: ' conf.init]);
end


if any(initS.weight == 0)
	error('Initialization produced a zero weight.');
end

mu = initS.mu;
sigma = initS.sigma;
weight = initS.weight;


log_initialmix = initS;
fixerloops = zeros(1, C);


% old values for stopping condition calculations
old_loglike = -realmax;

loops=1;
fixing_cycles = 0;

tulo = gmmcpdf(data, mu, sigma, weight);

while 1
	% one EM cycle
	pcompx = tulo ./ (sum(tulo,2)*ones(1,C));
	
	if ~all( isfinite(pcompx(:))  )
		error('Probabilities are no longer finite.');
	end
	
	for c = 1:C
		% calculate new estimates
		psum = sum(pcompx(:,c));
		
		% weight
		weight(c) = 1/N*psum;
	
		% mean
		nmu = sum(data.*(pcompx(:,c)*ones(1,D)), 1).' ./ psum;
		mu(:,c) = nmu;
		
		% covariance
		moddata = (data - ones(N,1)*(nmu.')) .* (sqrt(pcompx(:,c))*ones(1,D));
		% sqrt(pcompx) is because it will be squared back
		nsigma = (moddata' * moddata) ./ psum;
		
		% covariance matrix goodness assurance
		[sigma(:,:,c), fixerloops(1,c)] = gmmb_covfixer(nsigma);
		% covfixer may change the matrix so that log-likelihood
		% decreases. So, if covfixer changes something,
		% disable the stop condition. If going into infinite
		% fix/estimate -loop, quit.
	end
	
	% finish test
	tulo = gmmcpdf(data, mu, sigma, weight);
	loglike = sum(log(sum(tulo, 2)));
	
	if conf.verbose ~= 0
		disp([ 'log-likelihood diff ' num2str(loglike-old_loglike)  ' on round ' num2str(loops) ]);
	end
	
	if conf.logging>0
		log_covfixer2{loops} = fixerloops;
		log_loglikes{loops} = loglike;
	end
	if conf.logging>1
		log_mixtures{loops} = struct(...
			'weight', weight, ...
			'mu', mu, ...
			'sigma', sigma);
	end

	if any(fixerloops ~= 0)
		% if any cov's were fixed, increase count and
		% do not evaluate stopping threshold.
		fixing_cycles = fixing_cycles +1;
		if conf.verbose ~= 0
			disp(['fix cycle ' num2str(fixing_cycles) ...
			      ', fix loops ' num2str(fixerloops)]);
		end
	else
		% no cov's were fixed this round, reset the counter
		% and evaluate threshold.
		fixing_cycles = 0;
		if (abs(loglike/old_loglike -1) < conf.thr)
			break;
		end
	end
	
	if fixing_cycles > 20
		warning_wrap('gmmb_em:fixing_loop', ...
		       ['A covariance matrix has been fixed repeatedly' ...
		        ' too many times, quitting EM estimation.']);
		break;
	end
	
	if loops >= conf.maxloops
		break;
	end

	loops = loops +1;
	old_loglike = loglike;
end



estimate = struct('mu', mu,...
		'sigma', sigma,...
		'weight', weight);

if conf.logging>1
	varargout{1} = struct(...
		'iterations', {loops}, ...
		'covfixer2', {cat(1,log_covfixer2{:})}, ...
		'loglikes', {cat(1,log_loglikes{:})}, ...
		'initialmix', {log_initialmix}, ...
		'mixtures', {log_mixtures});
end
if conf.logging == 1
	varargout{1} = struct(...
		'iterations', {loops}, ...
		'covfixer2', {cat(1,log_covfixer2{:})}, ...
		'loglikes', {cat(1,log_loglikes{:})} );
end


% ------------------------------------------

function tulo = gmmcpdf(data, mu, sigma, weight);
N = size(data, 1);
C = size(weight,1);

pxcomp = zeros(N,C);
for c = 1:C
	pxcomp(:,c) = gmmb_cmvnpdf(data, mu(:,c).', sigma(:,:,c));
end
tulo = pxcomp.*repmat(weight.', N,1);

