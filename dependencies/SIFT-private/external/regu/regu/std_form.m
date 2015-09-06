function [A_s,b_s,L_p,K,M] = std_form(A,L,b,W)
%STD_FORM Transform a general-form reg. problem into one in standard form.
%
% [A_s,b_s,L_p,K,M] = std_form(A,L,b)      (method 1)
% [A_s,b_s,L_p,x_0] = std_form(A,L,b,W)    (method 2)
%
% Transforms a regularization problem in general form
%    min { || A x - b ||^2 + lambda^2 || L x ||^2 }
% into one in standard form
%    min { || A_s x_s - b_s ||^2 + lambda^2 || x_s ||^2 } .
%
% Two methods are available.  In both methods, the regularized
% solution to the original problem can be written as
%    x = L_p*x_s + d
% where L_p and d depend on the method as follows:
%    method = 1: L_p = pseudoinverse of L, d  = K*M*(b - A*L_p*x_s)
%    method = 2: L_p = A-weighted pseudoinverse of L, d = x_0.
%
% The transformation from x_s back to x can be carried out by means
% of the subroutine gen_form.

% Reference: P. C. Hansen, "Rank-Deficient and Discrete Ill-PosedProblems.
% Numerical Aspects of Linear Inversion", SIAM, Philadelphia, 1997.

% Per Christian Hansen, IMM, April 8, 2001.

% Nargin determines which method.
if (nargin==3)

  % Initialization for method 1.
  [m,n] = size(A); [p,np] = size(L);
  if (np~=n), error('A and L must have the same number of columns'), end

  % Special treatment of the case where L is square.
  if (p==n)
    L_p = inv(L); K = []; M = []; A_s = A/L; b_s = b;
    return
  end

  % Compute a QR factorization of L'.
  [K,R] = qr(full(L')); R = R(1:p,:);

  % Compute a QR factorization of A*K(:,p+1:n)).
  [H,T] = qr(A*K(:,p+1:n)); T = T(1:n-p,:);

  % Compute the transformed quantities.
  L_p = K(:,1:p)/R';   %(R\(K(:,1:p)'))';
  K   = K(:,p+1:n);
  M   = T\(H(:,1:n-p)');
  A_s = H(:,n-p+1:m)'*A*L_p;
  b_s = H(:,n-p+1:m)'*b;

else

  % Initialization for method 2.
  [m,n] = size(A); [p,nl] = size(L); nu = n-p;
  if (nl~=n), error('A and L must have the same number of columns'), end

  % Special treatment of the case where L is square.
  if (p==n)
    L_p = inv(L); A_s = A/L; b_s = b;
    x_0 = zeros(n,1); K = x_0; % Fix output name.
    return
  end

  % Compute T and x_0;
  [T,x_0] = pinit(W,A,b);
  b_s = b - A*x_0;

  % Compute the transformed quantities.
  L1 = inv(L(:,1:p));
  L_p = [L1;zeros(nu,p)] - W*(T(:,1:p)*L1);
  A_s = A*L_p;

  % Fix output name.
  K = x_0;

end