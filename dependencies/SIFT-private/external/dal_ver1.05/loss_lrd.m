% loss_lrd - conjugate logistic loss function
% 
% Syntax:
% [floss, gloss, hloss, hmin]=loss_lrd(aa, yy)
%
% Copyright(c) 2009 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt
function varargout = loss_lrd(aa, yy)

mm=length(aa);

gloss=nan*ones(mm,1);

ya = aa.*yy;

I = find(0<ya & ya<1);


floss = sum((1-ya(I)).*log(1-ya(I))+ya(I).*log(ya(I)));
gloss(I) = yy(I).*log(ya(I)./(1-ya(I)));

hmin  = 4;

if nargout<=3
  varargout={floss, gloss, hmin};
else
  hloss = spdiag(1./(ya.*(1-ya)));
  varargout={floss, gloss, hloss, hmin};
end


