% randsparse - generates a random sparse vector or a column-wise
%              sparse matrix
%
% Example:
%  ww = randsparse(64, 8);
%  ww = randsparse([64, 64], 8);
%
% Copyright(c) 2009 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt
function ww = randsparse(n, k)

if length(n)==1
  I=randperm(n);
  ww=zeros(n,1);
  ww(I(1:k))=randn(k,1);
else
  ww=zeros(n);
  I=randperm(n(2));
  ww(:,I(1:k))=randn(n(1),k);
end
