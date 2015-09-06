% l1_softth - soft threshold function for L1 regularization
%
% Copyright(c) 2009 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt

function [vv,ss]=l1_softth(vv,lambda,info)

n = size(vv,1);

Ip=find(vv>lambda);
In=find(vv<-lambda);

vv=sparse([Ip;In],1,[vv(Ip)-lambda;vv(In)+lambda],n,1);

ss=abs(vv);