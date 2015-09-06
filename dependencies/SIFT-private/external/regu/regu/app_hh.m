function A = app_hh(A,beta,v)
%APP_HH Apply a Householder transformation.
%
% A = app_hh(A,beta,v)
%
% Applies the Householder transformation, defined by
% vector v and scaler beta, to the matrix A; i.e.
%     A = (eye - beta*v*v')*A .

% Per Christian Hansen, IMM, 03/11/92.

A = A - (beta*v)*(v'*A);
