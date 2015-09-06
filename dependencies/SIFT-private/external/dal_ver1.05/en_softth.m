% en_softth - soft threshold function for the Elastic-net regularization
%
% Copyright(c) 2009 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt

function [vv,ss]=en_softth(vv,lambda,info)

n = size(vv,1);
theta  = info.theta;

if theta<1
  I=find(abs(vv)>lambda*theta);
  vv=sparse(I,1,(abs(vv(I))-lambda*theta).*sign(vv(I))/(1+lambda*(1-theta)),n,1);
else
  vv=l1_softth(vv,lambda,info);
end

ss=en_spec(vv,theta);