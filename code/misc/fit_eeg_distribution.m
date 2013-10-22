function [mu,sig] = fit_eeg_distribution(x,min_clean_fraction,max_dropout_fraction,fit_quantiles,step_sizes,num_bins,ensemble_quantile) %#ok<INUSD>
% Fit a truncated Gaussian to possibly contaminated EEG data.
% [Mu,Sigma] = fit_eeg_distribution(X,MinCleanFraction,MaxDropoutFraction,FitQuantiles,StepSizes,NumBins,EnsembleQuantile)
%
% This function assumes that the observations are EEG amplitude values that can be characterized
% as a Gaussian component (the clean data) with a heavy and potentially peaky tail, and possibly a
% relatively small fraction of samples with amplitude below that of the Gaussian component.
%
% The method works by finding a quantile range that best fits a truncated Gaussian (in terms of KL
% divergence), using a grid search over a range of the data that is restricted by MinCleanFraction
% and MaxDropoutFraction.
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
%   FitQuantiles : Quantile range [lower,upper] of the truncated Gaussian distribution
%                  that shall be fit to the EEG contents, or ensemble of quantile ranges under 
%                  consideration (Nx2 matrix) (default: [0.022*ones(21,1) (0.5:0.01:0.7)'])
%
%   StepSizes : Step size of the grid search; the first value is the stepping of the lower bound
%               (which essentially steps over any dropout samples), and the second value
%               is the stepping over possible scales (i.e., clean-data quantiles)
%               (default: [0.01 0.01])
%
%   NumBins : Deprecated argument (was the number of bins, by now auto-determined).
%
%   EnsembleQuantile : If multiple FitQuantiles are passed in, this is the quantile of the resulting
%                      mean estimates that shall be returned (together with its associated sigma);
%                      for example, a value of 0.5 would yield the median fit of the ensemble, while
%                      a value of 1.0 would yield the fit with largest mean among the ensemble 
%                      (default: 0.8)
%
% Out:
%   Mu : estimated mean of the distribution
%
%   Sigma : estimated standard deviation of the distribution
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2013-08-15

% Copyright (C) Christian Kothe, SCCN, 2013, christian@sccn.ucsd.edu
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
if ~exist('fit_quantiles','var') || isempty(fit_quantiles)
    fit_quantiles = [0.022*ones(21,1) (0.5:0.01:0.7)']; end
if ~exist('step_sizes','var') || isempty(step_sizes)
    step_sizes = [0.01 0.01]; end
if ~exist('ensemble_quantile','var') || isempty(ensemble_quantile)
    ensemble_quantile = 0.8; end

% sort data so we can access quantiles directly
x = sort(x(:));
n = length(x);

% iterate over ensemble of quantile ranges
all_quantiles = fit_quantiles;
for k=1:size(all_quantiles,1)
    fit_quantiles = all_quantiles(k,:);
    
    % calculate bounds for the truncated Gaussian pdf (in standard deviations)
    bounds = -sqrt(2)*erfcinv(2*[min(fit_quantiles) max(fit_quantiles)]);
    
    % determine the limits for the grid search
    lower_min = min(fit_quantiles);                     % we can generally skip the tail below the lower quantile
    min_width = min_clean_fraction*diff(fit_quantiles); % minimum width of the fit interval, as fraction of data
    max_width = diff(fit_quantiles);                    % maximum width is the fit interval if all data is clean
    
    opt_val = Inf;
    % for each interval width
    for width = min_width : step_sizes(2) : max_width
        inds = (1:round(n*width))';
        offsets = round(n*(lower_min : step_sizes(1) : lower_min+max_dropout_fraction));
        
        % get data in shifted intervals
        T = x(bsxfun(@plus,inds,offsets));
        
        % calculate binned histograms
        bins = round(3*log2(1+length(inds)/2));
        p = exp(-0.5*(bounds(1)+(0.5:(bins-0.5))/bins*diff(bounds)).^2)/(sqrt(2*pi)); p=p'/sum(p);
        q = histc(bsxfun(@times,bsxfun(@minus,T,T(1,:)),1./(T(end,:)-T(1,:))),[0:bins-1,Inf]/bins) + 0.01;
        
        % calc KL divergences
        kl = sum(bsxfun(@times,p,log(bsxfun(@rdivide,p,q(1:end-1,:))))) + log(length(inds));
        
        % update optimal range
        [min_val,idx] = min(kl);
        if min_val < opt_val
            opt_val = min_val;
            opt_lu = T([1,end],idx);
        end
    end
    
    % recover mu and sigma from optimal bounds
    sig(k) = (opt_lu(2)-opt_lu(1))/diff(bounds);
    mu(k) = opt_lu(1)-bounds(1)*sig(k);
end

% select ensemble representative
if ~isvector(all_quantiles)
    [dummy,idx] = sort(mu,'ascend'); %#ok<ASGLU>
    pick = idx(round(ensemble_quantile*end));
    mu = mu(pick);
    sig = sig(pick);
end

