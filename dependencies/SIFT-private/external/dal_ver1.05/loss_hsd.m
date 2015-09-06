% loss_hsd - conjugate hyperbolic secant loss function
% 
% Syntax:
% [floss, gloss, hloss, hmin]=loss_hsd(aa, yy)
%
% Copyright(c) 2009-2011 Ryota Tomioka
%              2009      Stefan Haufe
% This software is distributed under the MIT license. See license.txt
function varargout=loss_hsd(zz, bb)

m=length(bb);

gloss = atanh(zz)-bb;
% floss = sum(zz.*(atanh(zz)-bb)+0.5*log(1-zz.^2)-log(pi));
floss = sum(-zz.*bb)+m*log(2/pi)+negent((1-zz)/2);

hmin = 2;

if nargout<=3
  varargout={floss, gloss, hmin};
else
  hloss = spdiag(1./(2*(1-zz.^2)));
  varargout={floss, gloss, hloss, hmin};
end

function out=negent(pp)
I = find(0<pp & pp<1);
out = sum((1-pp(I)).*log(1-pp(I))+pp(I).*log(pp(I)));
