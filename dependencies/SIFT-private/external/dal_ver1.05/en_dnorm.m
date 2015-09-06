% en_dnorm - conjugate of the Elastic-net regularizer
%
% Copyright(c) 2009 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt
function [nm,ishard]=en_dnorm(ww,lambda,theta)

if theta<1
  ishard=0;
  nm = 0.5*sum(max(0,abs(ww)-lambda*theta).^2)/(lambda*(1-theta));
else
  ishard=1;
  nm = max(abs(ww));
end

