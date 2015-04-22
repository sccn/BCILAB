function [wmean,wcov]=weightmecov(data,weights)

%WEIGHTMECOV computes the reweighted mean and covariance matrix of multivariate data.
% 
% Required input arguments:
%      data : data matrix
%   weights : weights of the observations
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
%Last update: 29 july 2005 

n = size(data,1);
nvar = size(data,2);

if ~(isempty(find(weights<0)))
    error('The weights are negative');
end

if size(weights,1)==1
   weights=weights';
end

wmean=sum(diag(weights/sum(weights))*data);
wcov=((data - repmat(wmean,n,1))'*diag(weights)*(data - repmat(wmean,n,1)))/(sum(weights.^2)-1);
