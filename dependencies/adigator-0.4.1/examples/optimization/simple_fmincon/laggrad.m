function GL = laggrad(x,lambda)

% Build input for x
X.f  = x;
X.dx = [1;1];
% Get objective gradient
F  = objfun_x(X);
GO = F.dx.';
% Get constraint jacobian
C  = confun_x(X);
JC = reshape(C.dx,C.dx_size);

GL = GO + lambda.'*JC;