function H = laghess2(x,lambda)
X.f  = x;         % Set Inputs
X.dx = [1;1];
GL = laggrad_x(X,lambda.ineqnonlin);
% Hessian for this is full so we can just reshape L.dxdx
H = reshape(GL.dx,GL.dx_size);