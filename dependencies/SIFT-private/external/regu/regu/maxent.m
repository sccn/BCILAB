function [x_lambda,rho,eta,data,X] = maxent(A,b,lambda,w,x0)
%MAXENT Maximum entropy regularization.
%
% [x_lambda,rho,eta] = maxent(A,b,lambda,w,x0)
%
% Maximum entropy regularization:
%    min { || A x - b ||^2 + lambda^2*x'*log(diag(w)*x) } ,
% where -x'*log(diag(w)*x) is the entropy of the solution x.
% If no weights w are specified, unit weights are used.
%
% If lambda is a vector, then x_lambda is a matrix such that
%    x_lambda = [x_lambda(1), x_lambda(2), ... ] .
%
% This routine uses a nonlinear conjugate gradient algorithm with "soft"
% line search and a step-length control that insures a positive solution.
% If the starting vector x0 is not specified, then the default is
%    x0 = norm(b)/norm(A,1)*ones(n,1) .

% Per Christian Hansen, IMM and Tommy Elfving, Dept. of Mathematics,
% Linkoping University, 06/10/92.

% Reference: R. Fletcher, "Practical Methods for Optimization",
% Second Edition, Wiley, Chichester, 1987.

% Set defaults.
flat = 1e-3;     % Measures a flat minimum.
flatrange = 10;  % How many iterations before a minimum is considered flat.
maxit = 150;     % Maximum number of CG iterations.
minstep = 1e-12; % Determines the accuracy of x_lambda.
sigma = 0.5;     % Threshold used in descent test.
tau0 = 1e-3;     % Initial threshold used in secant root finder.

% Initialization.
[m,n] = size(A); x_lambda = zeros(n,length(lambda)); F = zeros(maxit,1);
if (min(lambda) <= 0)
  error('Regularization parameter lambda must be positive')
end
if (nargin ==3), w  = ones(n,1); end
if (nargin < 5), x0 = ones(n,1); end

% Treat each lambda separately.
for j=1:length(lambda);

  % Prepare for nonlinear CG iteration.
  l2 = lambda(j)^2;
  x  = x0; Ax = A*x;
  g  = 2*A'*(Ax - b) + l2*(1 + log(w.*x));
  p  = -g;
  r  = Ax - b;

  % Start the nonlinear CG iteration here.
  delta_x = x; dF = 1; it = 0; phi0 = p'*g;
  while (norm(delta_x) > minstep*norm(x) & dF > flat & it < maxit & phi0 < 0)
    it = it + 1;

    % Compute some CG quantities.
    Ap = A*p; gamma = Ap'*Ap; v = A'*Ap;

    % Determine the steplength alpha by "soft" line search in which
    % the minimum of phi(alpha) = p'*g(x + alpha*p) is determined to
    % a certain "soft" tolerance.
    % First compute initial parameters for the root finder.
    alpha_left = 0; phi_left = phi0;
    if (min(p) >= 0)
      alpha_right = -phi0/(2*gamma);
      h = 1 + alpha_right*p./x;
    else
      % Step-length control to insure a positive x + alpha*p.
      I = find(p < 0);
      alpha_right = min(-x(I)./p(I));
      h = 1 + alpha_right*p./x; delta = eps;
      while (min(h) <= 0)
        alpha_right = alpha_right*(1 - delta);
        h = 1 + alpha_right*p./x;
        delta = delta*2;
      end
    end
    z = log(h);
    phi_right = phi0 + 2*alpha_right*gamma + l2*p'*z;
    alpha = alpha_right; phi = phi_right;

    if (phi_right <= 0)

      % Special treatment of the case when phi(alpha_right) = 0.
      z = log(1 + alpha*p./x);
      g_new = g + l2*z + 2*alpha*v; t = g_new'*g_new;
      beta = (t - g'*g_new)/(phi - phi0);

    else

      % The regular case: improve the steplength alpha iteratively
      % until the new step is a descent step.
      t = 1; u = 1; tau = tau0;
      while (u > -sigma*t)

        % Use the secant method to improve the root of phi(alpha) = 0
        % to within an accuracy determined by tau.
        while (abs(phi/phi0) > tau)
          alpha = (alpha_left*phi_right - alpha_right*phi_left)/...
                  (phi_right - phi_left);
          z = log(1 + alpha*p./x);
          phi = phi0 + 2*alpha*gamma + l2*p'*z;
          if (phi > 0)
            alpha_right = alpha; phi_right = phi;
          else
            alpha_left  = alpha; phi_left  = phi;
          end
        end

        % To check the descent step, compute u = p'*g_new and
        % t = norm(g_new)^2, where g_new is the gradient at x + alpha*p.
        g_new = g + l2*z + 2*alpha*v; t = g_new'*g_new;
        beta = (t - g'*g_new)/(phi - phi0);
        u = -t + beta*phi;
        tau = tau/10;

      end  % End of improvement iteration.

    end  % End of regular case.

    % Update the iteration vectors.
    g = g_new; delta_x = alpha*p;
    x = x + delta_x;
    p = -g + beta*p;
    r = r + alpha*Ap;
    phi0 = p'*g;

    % Compute some norms and check for flat minimum.
    rho(j,1) = norm(r); eta(j,1) = x'*log(w.*x);
    F(it) = rho(j,1)^2 + l2*eta(j,1);
    if (it <= flatrange)
      dF = 1;
    else
      dF = abs(F(it) - F(it-flatrange))/abs(F(it));
    end

    data(it,:) = [F(it),norm(delta_x),norm(g)];
    X(:,it) = x;

  end  % End of iteration for x_lambda(j).

  x_lambda(:,j) = x;

end