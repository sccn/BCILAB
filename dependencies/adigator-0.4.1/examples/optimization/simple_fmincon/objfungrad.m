function [f,G] = objfungrad(x)
% Wrapper file for the objective function gradient - if fmincon wants the
% gradient, then we evaluate the adigator created file, otherwise we
% evaluate the objfun file.

if nargout == 1
  f = objfun(x);
else
  X.f = x;          % Create input to deriv file
  X.dx = [1;1];
  F = objfun_x(X);  % Call deriv file
  f = F.f;
  % In this case the gradient is full, so we can just transpose the
  % derivative output
  G = F.dx.';
end