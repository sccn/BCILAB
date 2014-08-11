function H = laghess1(x,lambda)
X.f  = x;         % Set Inputs
X.dx = [1;1];
L = lagrangian1_xx(X,lambda.ineqnonlin);
% Hessian for this is full so we can just reshape L.dxdx
H = reshape(L.dxdx,L.dxdx_size);