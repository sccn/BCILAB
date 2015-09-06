% l1n_spec - spectrum function for the non-negative L1 regularizer
% 
% Copyright(c) 2009-2011 Ryota Tomioka
%                   2011 Shigeyuki Oba
% This software is distributed under the MIT license. See license.txt

function ss=l1n_spec(ww)

n=size(ww,1);

Ip=find(ww>0); lenp=length(Ip);
In=find(ww<0); lenn=length(In);

ss=sparse([Ip;In],1,[ww(Ip);inf*ones(lenn,1)],n,1);
