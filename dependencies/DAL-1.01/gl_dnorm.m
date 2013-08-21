% gl_dnorm - conjugate of the grouped L1 regularizer
%
% Copyright(c) 2009 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt
function [nm,ishard]=gl_dnorm(ww,blks)

nm=0;
ix0=0;
for kk=1:length(blks)
  I=ix0+(1:blks(kk));
  ix0=I(end);
  nm=max(nm, norm(ww(I)));
end
ishard=1;