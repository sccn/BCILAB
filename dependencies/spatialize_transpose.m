function [X,dims] = spatialize_transpose(X)
% Collapse all extra dimensions of a time series matrix into the spatial (channel) dimension and
% transpose.
dims = size(X);
X = permute(X,[2 1 3:length(dims)]);
X = X(:,:);