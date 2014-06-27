function [X,dims] = spatialize(X)
% Collapse all extra dimensions of a time series matrix into the spatial (channel) dimension.
dims = size(X);
X = reshape(permute(X,[2 1 3:length(dims)]),dims(2),[])';
