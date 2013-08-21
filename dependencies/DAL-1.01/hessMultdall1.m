% hessMultdall1 - function that computes H*x for DAL with L1
%                 regularization
%
% Copyright(c) 2009 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt

function yy = hessMultdall1(xx, A, AT, eta, Hinfo)

hloss=Hinfo.hloss;
AF=Hinfo.AF;
I=Hinfo.I;
n=Hinfo.n;
len=length(I);
% F = sparse(I,ones(len,1), (xx'*AF)', n, 1);
yy = hloss*xx + eta*(AF*(AF'*xx));

B=Hinfo.B;
if ~isempty(B)
  yy = yy + eta*(B*(B'*xx));
end


