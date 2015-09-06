% hessMultdall1 - function that computes H*x for DAL with L1
%                 regularization
%
% Copyright(c) 2009 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt

function yy = hessMultdall1(xx, A, eta, Hinfo)

hloss=Hinfo.hloss;
AF=Hinfo.AF;
I=Hinfo.I;
n=Hinfo.n;
len=length(I);
yy = hloss*xx + eta(1)*(AF*(AF'*xx));

B=Hinfo.B;
if ~isempty(B)
  yy = yy + eta(2)*(B*(B'*xx));
end


