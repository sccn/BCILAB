% loss_sqd - conjugate squared loss function
%
% Syntax:
% [floss, gloss, hloss, hmin]=loss_sqd(aa, bb)
% 
% Copyright(c) 2009 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt
function varargout = loss_sqd(aa, bb)

gloss = aa-bb;
floss = 0.5*sum(gloss.^2)-0.5*sum(bb.^2);
hloss = spdiag(ones(size(aa)));
hmin  = 1;
  
if nargout<=3
  varargout = {floss, gloss, hmin};
else
  varargout = {floss, gloss, hloss, hmin};
end

