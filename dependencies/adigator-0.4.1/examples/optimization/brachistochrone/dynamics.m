function Xdot = dynamics(X,U)
% Brachistochrone Dynamics
Xdot = zeros(size(X));
X3 = X(:,3);
Xdot(:,1) = X3.*sin(U);
Xdot(:,2) = -X3.*cos(U);
Xdot(:,3) = 9.81.*cos(U);