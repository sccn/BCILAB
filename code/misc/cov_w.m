function wcov = cov_w(X, varargin)
% X is a matrix with a variable in each column, and one observation in each row. The vector weights contains weights for each row,
% and must thus have the same length as the number of rows of X.
%
% Dependencies: hlp_varargin2struct
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-03

% Copyright (C) Christian Kothe, SCCN, 2010, christian@sccn.ucsd.edu
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


n = size(X,1); % number of observations, i.e. rows of X

opts = hlp_varargin2struct(varargin, 'weights', ones(n, 1)/n);

if any(opts.weights<0)
    opts.weights(opts.weights<0) = 0;
    warning('cov_w:neg_weights','cov_w.m: Negative weights were set to zero')
end

if length(opts.weights) ~= size(X, 1)
    error('cov_w.m: The number of weights does not equal the number of observations')
end

if n < 2
   warning('cov_w:one_row','cov_w.m: X contains only one row, so meaningful covariances cannot be computed') 
end

opts.weights = opts.weights(:)'; % make sure that opts.weights is a row vector

% normalize weights
opts.weights = opts.weights./sum(opts.weights(:));
weights_mat = repmat(opts.weights', 1, size(X, 2));
% calculate weighted mean

weighted_X = weights_mat.*X;

col_means = sum(weighted_X, 1); % since the weights are normalized to add up to one, summing the weigthed rows gives the weighted mean

X_wmeans_subtracted = X - repmat(col_means, size(X, 1), 1);

sum_squared_weights = sum(opts.weights.^2);

wcov = X_wmeans_subtracted'*(weights_mat.*X_wmeans_subtracted);

wcov = wcov./(1-sum_squared_weights);
