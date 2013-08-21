%GMMB_FJ     - Figueiredo-Jain estimated GMM parameters
% Produces a bayesS struct without 'apriories'.
%
% Works with complex numbers directly.
%
% estimate = GMMB_FJ(data[, parameters])
% [estimate, stats] = GMMB_FJ(...)
%
% Parameters (default):
%   maxloops	maximum number of loops per a CEM run (500)
%   Cmax	the maximum number of GMM components
%		set to -1 to use all data points as component means
%   Cmin	the minimum number of GMM components tried (1)
%   verbose	print some progress numbers (false)
%   thr		CEM2 threshold, relative loglikelihood change (1e-6)
%   animate	plot data and ellipses as the algorithm runs (false)
%   covtype	Covariance matrix structure: 1=diagonal, other=free (0)
%   broken	With complex data, calculate no. of component parameters
%		as with real data (true).
% At least Cmax should be set explicitly.
% Example:
%   estS = gmmb_fj(data, 'Cmax', 50, 'thr', 1e-9)
%
% References:
%   [1] Duda, R.O., Hart, P.E, Stork, D.G, Pattern Classification,
%   2nd ed., John Wiley & Sons, Inc., 2001.
%   [2] Bilmes, J.A., A Gentle Tutorial of the EM Algorithm and its
%    Application to Parameter Estimation for Gaussian Mixture and Hidden
%    Markov Models
%   International Computer Science Institute, 1998
%   [3] Figueiredo, M.A.T., Jain, A.K., Unsupervised Learning on
%    Finite Mixture Models, IEEE transactions of pattern analysis and
%    machine intelligence, vol.24, no3, March 2002
%
% This code is directly based on [3] and code published on
% Figueiredo homepage: http://www.lx.it.pt/~mtf/
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
%   $Name:  $ $Revision: 1.2 $  $Date: 2004/11/02 09:00:18 $
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
%      iterations: CEM (full) iteration count
%      costs:      iterations long vector of the cost value
%      annihilations: component annihilation log [annih, deletion]*iters
%      covfixer2:  iterations-by-C matrix of gmmb_covfixer fix round counts
%      loglikes:   iterations long vector of the log-likelihood
%    extra logging:
%      initialmix: parameters for the initial mixture
%      mixtures:   parameters for all intermediate mixtures
%

function [estimate, varargout] = gmmb_fj(data, varargin);

[N, D] = size(data);	% number of points (n), dimensions (d)

% defaults
conf = struct(...
	'maxloops', 500, ...
	'Cmax', ceil(min(50, N/(D*D)/3)), ...
	'Cmin', 1, ...
	'verbose', 0, ...
	'thr', 1e-6, ...
	'animate', 0, ...
	'covtype', 0, ...
	'broken', 1, ...
	'logging', 0 ...
	);

if nargout>1
	conf.logging = 1;
	varargout{1} = [];
end

conf = getargs(conf, varargin);

C = conf.Cmax;

if nargout<2
	conf.logging=0;
end

% for logging
log_covfixer2 = {};
log_loglikes = {};
log_costs = {};
log_annih = {};
log_initialmix = {};
log_mixtures = {};


if (C<1) | (C>N)
	C = N;
	mu = data.';
else
	% initialize mu as random points from data
	permi = randperm(N);
	mu = data(permi(1:C),:).';  % D x C
end


% initialize sigma
s2 = max(diag(gmmb_covfixer(cov(data,1))/10));
sigma = repmat(s2*eye(D), [1 1 C]);

% weights initialization
alpha = ones(1,C) * (1/C);


log_initialmix = struct(...
	'weight', alpha, ...
	'mu', mu, ...
	'sigma', sigma);


% the number of free parameters in a Gaussian
if isreal(data) | (conf.broken ~= 0)
	if conf.covtype == 1
		Nparc = D+D;	% (N)
	else
		Nparc = D+D*(D+1)/2;	% (N)
	end
else
	% data is complex valued
	if conf.covtype == 1
		Nparc = 2*D + D;	% (N)
	else
		Nparc = 2*D + D*D;	% (N)
	end
end
Nparc2 = Nparc/2;

N_limit = (Nparc+1)*3*conf.Cmin;
if N < N_limit
	warning_wrap('gmmb_fj:data_amount', ...
	   ['Training data may be insufficient for selected ' ...
	    'minimum number of components. ' ...
	    'Have: ' num2str(N) ', recommended: >' num2str(N_limit) ...
	    ' points.']);
end


if conf.animate ~= 0
	aniH = my_plot_init;
	my_plot_ellipses(aniH, data, mu, sigma, alpha);
end


t = 0;
Cnz = C;	% (k_nz) k = kmax
Lmin = NaN;

u = zeros(N,C);	% semi_indic.'
for c = 1:C
	u(:,c) = gmmb_cmvnpdf(data, mu(:,c).', sigma(:,:,c));
end
indic = u .* repmat(alpha, N,1);

old_loglike = sum(log(sum(realmin+indic, 2)));
old_L = Nparc2*sum(log(alpha)) + (Nparc2+0.5)*Cnz*log(N) - old_loglike;


while Cnz >= conf.Cmin
	repeating = 1;
	
	fixing_cycles = 0;
	loops = 0;
	while repeating
		t = t+1;
		loops = loops +1;
		
		fixed_on_this_round = 0;
		log_covfixer2{t,1} = 0;
		c = 1;
		while c <= C
			indic = u .* repmat(alpha, N,1);
			normindic = indic ./ (realmin + repmat(sum(indic,2), 1,C));
			
			normf = 1/sum(normindic(:,c));
			aux = repmat(normindic(:,c), 1,D) .* data;
			
			nmu = normf * sum(aux,1);
			mu(:,c) = nmu.';
			
			if conf.covtype == 1
				nsigma =  normf*diag(sum(aux .* conj(data), 1)) - diag(nmu.*conj(nmu));
			else
				nsigma =  normf*(aux' * data) - nmu'*nmu;
			end
			[sigma(:,:,c) log_fixcount] = gmmb_covfixer(nsigma);
			% covfixer may change the matrix so that log-likelihood
			% decreases. So, if covfixer changes something,
			% disable the stop condition. If going into infinite
			% fix/estimate -loop, quit.

			if conf.logging>0
				% the component indexes are not constants,
				% cannot record component-wise fix counts
				log_covfixer2{t,1} = ...
				   log_covfixer2{t,1} + log_fixcount;
			end
			
			alpha(c) = max(0, sum(normindic(:,c))-Nparc2) / N;
			alpha = alpha / sum(alpha);
			
			if ~all( isfinite(alpha(:)) )
				% something went wrong
				% probably there is not enough data to
				%support estimation
				warning_wrap('gmmb_fj:weight_finity', 'Mixture weights are no longer finite, aborting estimation.');
				alpha(:) = 0;
				Cnz = 0;
				repeating = 0;
			end

			if alpha(c) == 0
				Cnz = Cnz -1;
			else
				if log_fixcount ~= 0
					% mark the covariance fix only,
					% if the component is not annihilated
					fixed_on_this_round = 1;
				end
				try
					u(:,c) = gmmb_cmvnpdf( data, ...
					   mu(:,c).', sigma(:,:,c) );
				catch
					disp('covariance went bzrk !!!');
					sigma(:,:,c)
					%keyboard
					Cnz = 0;
				end
			end
			c=c+1;
			
			if Cnz <= 0
				% number of components fell to zero
				% nothing can be done
				error('Estimation failed, number of components fell to zero. Not enough training data?');
			end

		end % while c <= C

		% purge alpha == 0 if necessary
		annihilated_count = length(find(alpha==0));
		if annihilated_count > 0
			nz = find(alpha>0);
			alpha = alpha(nz);
			mu = mu(:,nz);
			sigma = sigma(:,:,nz);
			u = u(:,nz);
			C = length(nz);
		end

		if conf.animate ~= 0
			my_plot_ellipses(aniH, data, mu, sigma, alpha);
		end
		
		u = zeros(N,C);	% semi_indic.'
		for c = 1:C
			u(:,c) = gmmb_cmvnpdf(data, mu(:,c).', sigma(:,:,c));
		end
		indic = u .* repmat(alpha, N,1);
		
		loglike = sum(log(realmin+sum(indic, 2)));
		L = Nparc2*sum(log(alpha)) + (Nparc2+0.5)*Cnz*log(N) - loglike;

		
		if conf.verbose ~= 0
			disp(['Cnz=' num2str(Cnz) ' t=' num2str(t) ' '...
			   num2str(abs(loglike - old_loglike)) ...
			   ' <? ' num2str(conf.thr*abs(old_loglike))]);
			disp(['t=' num2str(t) ' L= ' num2str(L)]);
		end
		
		if conf.logging>0
			log_loglikes{t} = loglike;
			log_costs{t} = L;
			log_annih{t} = [annihilated_count, 0];
		end
		if conf.logging>1
			log_mixtures{t} = struct(...
				'weight', alpha, ...
				'mu', mu, ...
				'sigma', sigma);
		end

		if fixed_on_this_round ~= 0
			% if any cov's were fixed, increase count and
			% do not evaluate stopping threshold.
			fixing_cycles = fixing_cycles +1;
			if conf.verbose ~= 0
				disp(['fix cycle ' num2str(fixing_cycles)]);
			end
		else
			% no cov's were fixed this round, reset the counter
			% and evaluate threshold.
			fixing_cycles = 0;
			if (abs(loglike/old_loglike -1) < conf.thr)
				repeating = 0;
			end
		end
		
		old_L = L;
		old_loglike = loglike;
		
		if fixing_cycles > 20
			repeating = 0;
		end
		if loops > conf.maxloops
			repeating = 0;
		end
	end % while repeating
	
	if isnan(Lmin) | (L <= Lmin)
		Lmin = L;
		estimate = struct('mu', mu,...
			'sigma', sigma,...
			'weight', alpha.');
	end
	if conf.verbose ~= 0
		disp(['Cnz = ' num2str(Cnz)]);
	end

	% annihilate the least probable component
	m = find(alpha == min(alpha(alpha>0)));
	alpha(m(1)) = 0;
	Cnz = Cnz -1;
	% alpha doesn't need to be normalized here, even if it would seem logical to do so.
	
	if conf.logging > 0
		log_annih{t}(2) = 1;
	end
	
	if Cnz > 0
		alpha = alpha / sum(alpha);
	
		% purge alpha == 0 if necessary
		if length(find(alpha==0)) > 0
			nz = find(alpha>0);
			alpha = alpha(nz);
			mu = mu(:,nz);
			sigma = sigma(:,:,nz);
			u = u(:,nz);
			C = length(nz);
		end
		
		u = zeros(N,C);	% semi_indic.'
		for c = 1:C
			u(:,c) = gmmb_cmvnpdf(data, mu(:,c).', sigma(:,:,c));
		end
		indic = u .* repmat(alpha, N,1);
		
		old_loglike = sum(log(realmin+sum(indic, 2)));
		old_L = Nparc2*sum(log(alpha)) + (Nparc2+0.5)*Cnz*log(N) - old_loglike;
	end
end



if conf.logging>1
	varargout{1} = struct(...
		'iterations', {t}, ...
		'costs', {cat(1,log_costs{:})}, ...
		'annihilations', {sparse(cat(1,log_annih{:}))}, ...
		'covfixer2', {cat(1,log_covfixer2{:})}, ...
		'loglikes', {cat(1,log_loglikes{:})}, ...
		'initialmix', {log_initialmix}, ...
		'mixtures', {log_mixtures});
end
if conf.logging == 1
	varargout{1} = struct(...
		'iterations', {t}, ...
		'costs', {cat(1,log_costs{:})}, ...
		'annihilations', {sparse(cat(1,log_annih{:}))}, ...
		'covfixer2', {cat(1,log_covfixer2{:})}, ...
		'loglikes', {cat(1,log_loglikes{:})} ...
		);
end


% purge alpha==0
e = estimate;
inds = find(e.weight>0);
estimate.mu = e.mu(:,inds);
estimate.sigma = e.sigma(:,:,inds);
estimate.weight = e.weight(inds);

if conf.animate ~= 0
	my_plot_ellipses(aniH, data, estimate.mu, estimate.sigma, estimate.weight);
end

%disp(['Cfinal = ' num2str(length(inds))]);

% -----------------------------------------------------------

function h = my_plot_init;
h = figure;
figure(h);
title('Distribution of x_1 and x_2 values','FontSize',14);
xlabel('x_1 value','FontSize',14);
ylabel('x_2 value','FontSize',14);
zlabel('weight','FontSize',14);
view(2)
tic;

function my_plot_ellipses(h, data, mu, sigma, weight);
dtime = 0.3;

D = size(mu, 1);

if D ~= 2
	error('Can plot only 2D objects.');
end

[x,y,z] = cylinder([2 2], 40);
xy = [ x(1,:) ; y(1,:) ];

figure(h);

plot(data(:,1), data(:,2), 'rx');

hold on
C = size(mu, 2);
for c = 1:C
	mxy = chol(sigma(:,:,c))' * xy;
	x = mxy(1,:) + mu(1,c);
	y = mxy(2,:) + mu(2,c);
	z = ones(size(x))*weight(c);
	plot3(x,y,z, 'k-');
end
drawnow;
hold off

t = toc;
if t+0.01<dtime
	pause(dtime-t);
end
tic

