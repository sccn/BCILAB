function x = splsqr(A,b,lambda,Vsp,maxit,tol,reorth)
%SPLSQR Subspace preconditioned LSQR for discrete ill-posed problems.
%
% x = splsqr(A,b,lambda,Vsp,maxit,tol,reorth)
%
% Subspace preconditioned LSQR (SP-LSQR) for solving the Tikhonov problem
%    min { || A x - b ||^2 + lambda^2 || x ||^2 }
% with a preconditioner based on the subspace defined by the columns of
% the matrix Vsp.  While not necessary, we recommend to use a matrix Vsp
% with orthonormal columns.
%
% The output x holds all the solution iterates as columns, and the last
% iterate x(:,end) is the best approximation to x_lambda.
%
% The parameter maxit is the maximum allowed number of iterations (default
% value is maxit = 300).  The parameter tol is used a stopping criterion
% for the norm of the least squares residual relative to the norm of the
% right-hand side (default value is tol = 1e-12).
%
% A seventh input parameter reorth ~= 0 enforces MGS reorthogonalization
% of the Lanczos vectors.

% This is a model implementation of SP-LSQR.  In a real implementation the
% Householder transformations should use LAPACK routines, only the final
% iterate should be returned, and reorthogonalization is not used.  Also,
% if Vsp represents a fast transformation (such as the DCT) then explicit
% storage of Vsp should be avoided.  See the reference for details.

% Reference: M. Jacobsen, P. C. Hansen and M. A. Saunders, "Subspace pre-
% conditioned LSQR for discrete ill-posed problems", BIT 43 (2003), 975-989.

% Per Christian Hansen and Michael Jacobsen, IMM, July 29, 2007.

% Input check.
if nargin < 5, maxit = 300; end
if nargin < 6, tol = 1e-12; end
if nargin < 7, reorth = 0; end
if maxit < 1, error('Number of iterations must be positive'); end;

% Prepare for SP-LSQR algorithm.
[m,n] = size(A);
k = size(Vsp,2);
z = zeros(n,1);
if reorth
    UU = zeros(m+n,maxit);
    VV = zeros(n,maxit);
end

% Initial QR factorization of [A;lamnda*eye(n)]*Vsp;
QQ = qr([A*Vsp;lambda*Vsp]);

% Prepare for LSQR iterations.
u = app_house_t(QQ,[b;z]);
u(1:k) = 0;
beta = norm(u);
u = u/beta;

v = app_house(QQ,u);
v = A'*v(1:m) + lambda*v(m+1:end);
alpha = norm(v);
v = v/alpha; 

w = v;
Wxw = zeros(n,1);

phi_bar = beta;
rho_bar = alpha;

if reorth, UU(:,1) = u; VV(:,1) = v; end;

for i=1:maxit

    % beta*u = A*v - alpha*u;
    uu = [A*v;lambda*v]; 
    uu = app_house_t(QQ,uu);
    uu(1:k) = 0;
    u = uu - alpha*u;
    if reorth
        for j=1:i-1, u = u - (UU(:,j)'*u)*UU(:,j); end
    end
    beta = norm(u);
    u = u/beta;

    % alpha * v = A'*u - beta*v;
    vv = app_house(QQ,u);
    v = A'*vv(1:m) + lambda*vv(m+1:end) - beta*v;
    if reorth
        for j=1:i-1, v = v - (VV(:,j)'*v)*VV(:,j); end
    end
    alpha = norm(v);
    v = v/alpha;

    if reorth, UU(:,i) = u; VV(:,i) = v; end;

    % Update LSQR parameters.
    rho = norm([rho_bar beta]);
    c = rho_bar/rho;
    s = beta/rho;
    theta = s*alpha;
    rho_bar = -c*alpha;
    phi = c*phi_bar;
    phi_bar = s*phi_bar;

    % Update the LSQR solution.
    Wxw = Wxw + (phi/rho)*w;
    w = v - (theta/rho)*w;

    % Compute residual and update the SP-LSQR iterate.
    r = [b - A*Wxw ; -lambda*Wxw];
    r = app_house_t(QQ,r);
    r = r(1:k);
    xv = triu(QQ(1:k,:))\r;
    x(:,i) = Vsp*xv + Wxw;

    % Stopping criterion.
    if phi_bar*alpha*abs(c) < tol*norm(b), break, end

end

%-----------------------------------------------------------------

function Y = app_house(H,X)
% Y = app_house(H,X)
% Input:  H = matrix containing the necessary information of the
%             Householder vectors v in the lower triangle and R in
%             the upper triangle; e.g., computed as H = qr(A).
%         X = matrix to be multiplied with orthogonal matrix.
% Output: Y = Q*X

[n,p] = size(H);
Y = X;
for k = p:-1:1
	v = ones(n+1-k,1);
	v(2:n+1-k) = H(k+1:n,k);
    beta = 2/(v'*v);
    Y(k:n,:) = Y(k:n,:) - beta*v*(v'*Y(k:n,:));
end

%-----------------------------------------------------------------

function Y = app_house_t(H,X)
% Y = app_house_t(H,X)
% Input:  H = matrix containing the necessary information of the
%             Householder vectors v in the lower triangle and R in
%             the upper triangle; e.g., computed as H = qr(A).
%         X = matrix to be multiplied with transposed orthogonal matrix.
% Output: Y = Q'*X

[n,p] = size(H);
Y = X;
for k = 1:p
	v = ones(n+1-k,1);
	v(2:n+1-k) = H(k+1:n,k);
    beta = 2/(v'*v);
    Y(k:n,:) = Y(k:n,:) - beta*v*(v'*Y(k:n,:));
end