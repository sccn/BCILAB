function [c,ceq,Jc,Jceq] = confunjac(x)
% Wrapper file for the constraint function gradient - if fmincon wants the
% jacobian, then we evaluate the adigator created file, otherwise we
% evaluate the confun file.

if nargout == 2
  [c,ceq] = confun(x);
else
  X.f  = x;              % Create input to deriv file
  X.dx = [1;1]; 
  [C,Ceq] = confun_x(X); % Call deriv file
  c   = C.f;
  ceq = Ceq.f;
  % In this case the Jacobian of C wrt x is full, so we can simply reshape
  % the output - Also transpose it since fmincon takes weird jacobians
  Jc   = reshape(C.dx,[2 2]).';
  Jceq = Ceq.dx;
end