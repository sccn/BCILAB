function [X,rho,eta] = rrgmres(A,b,k)
%RRGMRES Range-restricted GMRES algorithm for square inconsistent systems
%
% [X,rho,eta] = rrgmres(A,b,k)
%
% RRGMRES is a variant of the GMRES algorithm for square linear systems
% A x = b, with starting vector A*b (instead of b as in GMRES).  This
% function returns all k iterates, stored as columns of the matrix X.
% The solution norm and residual norm are returned in eta and rho, resp.
%
% For symmetric matrices, use the function mr2 instead.

% Reference: D. Calvetti, B. Lewis & L. Reichel, "GMRES-type methods for
% inconsistent systems", Lin. Alg. Appl. 316 (2000), 157-169.

% Per Christian Hansen, IMM, September 19, 2007.

% Check input arguments.
if (k < 1), error('Number of steps k must be positive'), end
[m,n] = size(A);
if (m ~= n), error('A must be square'), end
if nargin < 3, error('Too few input arguments'), end

% Allocate space.
X = zeros(n,k); % Matrix of solutions.
V = zeros(n,k); % Orthonormal vectors spanning Krylov subspace
h = zeros(k,1); % New column of Hessenberg matrix
Q = zeros(k+1); % H = Q*T, Q orthogonal
T = zeros(k);   %          T upper triangular
W = zeros(n,k); % W = V*inv(T)
if (nargout>1), rho = zeros(k,1); end
if (nargout>2), eta = zeros(k,1); end

% Initialize variables.
r      = A*b;
alpha  = norm(r);
V(:,1) = r/alpha;  % Initial vector of Krylov subspace is A*b.
Q(1,1) = 1;
x      = zeros(n,1);
beta   = V(:,1)'*b;

% Begin iterations.
for i=1:k
    
    r = A*V(:,i);
    
    % Modified Gram-Schmidt on the new vector.
    for j=1:k
        h(j) = V(:,j)'*r;
        r = r - V(:,j)*h(j);
    end
    alpha = norm(r);
    
    % Store new Arnoldi vector and update projected rhs.
    V(:,i+1) = r/alpha;
    beta = [beta;V(:,i+1)'*b];
    
    % Apply previous rotations to h.
    T(1:i,i) = Q(1:i,1:i)'*h(1:i);
    
    % Compute Givens rotation parameters.
    rc = T(i,i);
    if alpha == 0
        c = 1; s = 0;
    elseif abs(alpha) > abs(rc)
        tau = -rc/alpha;
        s = 1 / sqrt(1 + abs(tau)^2);
        c = s*tau;
    else
        tau = -alpha/rc;
        c = 1 / sqrt(1 + abs(tau)^2);
        s = c*tau;
    end
    
    % Apply givens rotations.
    T(i,i) = c'*rc - s'*alpha;
    Q(1:i,[i,i+1]) = Q(1:i,i)*[c s];
    Q(i+1,[i,i+1]) = [-s c];
        
    if abs(T(i,i)) <= eps
        disp('Hession matrix is (numerically) singular')
    end
    
    % Update W = V*inv(T);
    W(:,i) = (V(:,i) - W(:,1:i-1)*T(1:i-1,i))/T(i,i);
    
    % Update solution.
    x = x + (Q(1:i+1,i)'*beta)*W(:,i);
    
    % Update output variables.
    X(:,i) = x;
    if nargout>1
        if i==1
            rho(i) = sqrt( b'*b - abs(Q(1:i+1,i)'*beta)^2 );
        else
            rho(i) = sqrt( rho(i-1)^2 - abs(Q(1:i+1,i)'*beta)^2 );
        end
    end
    if nargout==3, eta(i) = norm(x); end
    
end