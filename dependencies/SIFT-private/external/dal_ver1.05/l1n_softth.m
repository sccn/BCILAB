% l1n_softth - soft threshold function for non-negative L1 regularization
%
% Copyright(c) 2009-2011 Ryota Tomioka
%                   2011 Shigeyuki Oba
% This software is distributed under the MIT license. See license.txt

function [vv,ss]=l1n_softth(vv,lambda,info)

n = size(vv,1);

I=find(vv>lambda);

vv=sparse(I,1,[vv(I)-lambda],n,1);

ss=abs(vv);