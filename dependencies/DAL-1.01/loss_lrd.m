% loss_lrd - conjugate logistic loss function
%
% Syntax:
% [floss, gloss, hloss, hmin]=loss_lrd(aa, yy)
%
% Copyright(c) 2009 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt
function varargout = loss_lrd(aa, yy)
% in the dual formulation, we return te funtion value, gradient, (optionally) hessian, and some other value

mm=length(aa);

gloss=nan*ones(mm,1);

ya = aa.*yy;

I = find(0<ya & ya<1);

% this may have something to do with the Fenchel dual of the logistic loss.??
% (note: understand this derivation!)
floss = sum((1-ya(I)).*log(1-ya(I))+ya(I).*log(ya(I)));
gloss(I) = yy(I).*log(ya(I)./(1-ya(I)));

hmin  = 4;

% Hessian is optional
if nargout<=3
    varargout={floss, gloss, hmin};
else
    hloss = spdiag(1./(ya.*(1-ya)));
    varargout={floss, gloss, hloss, hmin};
end


