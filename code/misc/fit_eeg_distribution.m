function [mu,sig,alpha,beta] = fit_eeg_distribution(X,min_clean_fraction,max_dropout_fraction,quants,step_sizes,beta)
% Estimate the mean and standard deviation of clean EEG from contaminated data.
% [Mu,Sigma,Alpha,Beta] = fit_eeg_distribution(X,MinCleanFraction,MaxDropoutFraction,FitQuantiles,StepSizes,ShapeRange)
%
% This function estimates the mean and standard deviation of clean EEG from a sample of amplitude
% values (that have preferably been computed over short windows) that may include a large fraction
% of contaminated samples. The clean EEG is assumed to represent a generalized Gaussian component in
% a mixture with near-arbitrary artifact components. By default, at least 25% (MinCleanFraction) of
% the data must be clean EEG, and the rest can be contaminated. No more than 10%
% (MaxDropoutFraction) of the data is allowed to come from contaminations that cause lower-than-EEG
% amplitudes (e.g., sensor unplugged). There are no restrictions on artifacts causing
% larger-than-EEG amplitudes, i.e., virtually anything is handled (with the exception of a very
% unlikely type of distribution that combines with the clean EEG samples into a larger symmetric
% generalized Gaussian peak and thereby "fools" the estimator). The default parameters should be
% fine for a wide range of settings but may be adapted to accomodate special circumstances.
%
% The method works by fitting a truncated generalized Gaussian whose parameters are constrained by
% MinCleanFraction, MaxDropoutFraction, FitQuantiles, and ShapeRange. The alpha and beta parameters
% of the gen. Gaussian are also returned. The fit is performed by a grid search that always finds a
% close-to-optimal solution if the above assumptions are fulfilled.
%
% In:
%   X : vector of amplitude values of EEG, possible containing artifacts
%       (coming from single samples or windowed averages)
%
%   MinCleanFraction : Minimum fraction of values in X that needs to be clean
%                      (default: 0.25)
%
%   MaxDropoutFraction : Maximum fraction of values in X that can be subject to
%                        signal dropouts (e.g., sensor unplugged) (default: 0.1)
%
%   FitQuantiles : Quantile range [lower,upper] of the truncated generalized Gaussian distribution
%                  that shall be fit to the EEG contents (default: [0.022 0.6])
%
%   StepSizes : Step size of the grid search; the first value is the stepping of the lower bound
%               (which essentially steps over any dropout samples), and the second value
%               is the stepping over possible scales (i.e., clean-data quantiles)
%               (default: [0.01 0.01])
%
%   ShapeRange : Range that the clean EEG distribution's shape parameter beta may take (default:
%                1.7:0.15:3.5)
%
% Out:
%   Mu : estimated mean of the clean EEG distribution
%
%   Sigma : estimated standard deviation of the clean EEG distribution
%
%   Alpha : estimated scale parameter of the generalized Gaussian clean EEG distribution (optional)
%
%   Beta : estimated shape parameter of the generalized Gaussian clean EEG distribution (optional)
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2013-08-15

% Copyright (C) Christian Kothe, SCCN, 2013, christiankothe@gmail.com
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License as published by the Free Software Foundation; either version 2 of the
% License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
% even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License along with this program; if not,
% write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
% USA

% assign defaults
if ~exist('min_clean_fraction','var') || isempty(min_clean_fraction)
    min_clean_fraction = 0.25; end
if ~exist('max_dropout_fraction','var') || isempty(max_dropout_fraction)
    max_dropout_fraction = 0.1; end
if ~exist('quants','var') || isempty(quants)
    quants = [0.022 0.6]; end
if ~exist('step_sizes','var') || isempty(step_sizes)
    step_sizes = [0.01 0.01]; end
if ~exist('beta','var') || isempty(beta)
    beta = 1.7:0.15:3.5; end

% sanity checks
if ~isvector(quants) || numel(quants) > 2
    error('Fit quantiles needs to be a 2-element vector (support for matrices deprecated).'); end
if any(quants(:)<0) || any(quants(:)>1)
    error('Unreasonable fit quantiles.'); end
if any(step_sizes<0.0001) || any(step_sizes>0.1)
    error('Unreasonable step sizes.'); end
if any(beta>=7) || any(beta<=1)
    error('Unreasonable shape range.'); end

% sort data so we can access quantiles directly
X = double(sort(X(:)));
n = length(X);

% calc z bounds for the truncated standard generalized Gaussian pdf and pdf rescaler
for b=1:length(beta)    
    zbounds{b} = sign(quants-1/2).*gammaincinv(sign(quants-1/2).*(2*quants-1),1/beta(b)).^(1/beta(b)); %#ok<*AGROW>
    rescale(b) = beta(b)/(2*gamma(1/beta(b)));
end

% determine the quantile-dependent limits for the grid search
lower_min = min(quants);                    % we can generally skip the tail below the lower quantile
max_width = diff(quants);                   % maximum width is the fit interval if all data is clean
min_width = min_clean_fraction*max_width;   % minimum width of the fit interval, as fraction of data

% get matrix of shifted data ranges
X = X(bsxfun(@plus,(1:round(n*max_width))',round(n*(lower_min:step_sizes(1):lower_min+max_dropout_fraction))));
X1 = X(1,:); X = bsxfun(@minus,X,X1);

opt_val = Inf;
% for each interval width...
for m = round(n*(max_width:-step_sizes(2):min_width))
    % scale and bin the data in the intervals
    nbins = round(3*log2(1+m/2));
    H = bsxfun(@times,X(1:m,:),nbins./X(m,:));
    logq = log(histc(H,[0:nbins-1,Inf]) + 0.01);
    
    % for each shape value...
    for b=1:length(beta)
        bounds = zbounds{b};
        % evaluate truncated generalized Gaussian pdf at bin centers
        x = bounds(1)+(0.5:(nbins-0.5))/nbins*diff(bounds);
        p = exp(-abs(x).^beta(b))*rescale(b); p=p'/sum(p);
        
        % calc KL divergences
        kl = sum(bsxfun(@times,p,bsxfun(@minus,log(p),logq(1:end-1,:)))) + log(m);
        
        % update optimal parameters
        [min_val,idx] = min(kl);
        if min_val < opt_val
            opt_val = min_val;
            opt_beta = beta(b);
            opt_bounds = bounds;
            opt_lu = [X1(idx) X1(idx)+X(m,idx)];
        end
    end
end

% recover distribution parameters at optimum
alpha = (opt_lu(2)-opt_lu(1))/diff(opt_bounds);
mu = opt_lu(1)-opt_bounds(1)*alpha;
beta = opt_beta;

% calculate the distribution's standard deviation from alpha and beta
sig = sqrt((alpha^2)*gamma(3/beta)/gamma(1/beta));
