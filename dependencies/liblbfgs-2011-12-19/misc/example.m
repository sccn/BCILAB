%
% Script to exemplify the use of the Matlab port of the liblbfgs-library
%

% These are all available options
options = ...         
struct('m',6,...                % The number of corrections to approximate the inverse hessian matrix.
    'epsilon',1e-7,...          % A minimization terminates when ||g|| < epsilon * max(1, ||x||)
    'past',0,...                % Distance for delta-based convergence test
    'delta',1e-5,...            % Delta for convergence test.
    'MaxIter',500,...           %The maximum number of iterations.
    'linesearch','more_thuente',...  % The line search algorithm (other options are 'backtracking_armijo', 'backtracking', 'backtracking_wolfe' or 'backtracking_strong_wolfe')
    'max_linesearch',40,...     % The maximum number of trials for the line search.
    'min_step',1e-20,...        % The maximum step of the line search.
    'max_step',1e20,...         % The minimum step of the line search routine.
    'ftol',1e-4,...             % A parameter to control the accuracy of the line search routine.
    'wolfe',0.9,...             % A coefficient for the Wolfe condition.
    'gtol',0.9,...              % A parameter to control the accuracy of the line search routine.
    'xtol',1e-16,...            % The machine precision for floating-point values.
    'orthantwise_c',0,...   	% Coefficient for the L1 norm of variables.
    'orthantwise_start',0,...   % First index for the parameters subject to L1-penalty
    'orthantwise_end',-1,...    % Last index for the parameters subject to L1-penalty
    'DerivativeCheck','off',... % Derivative check using finite differences  ('on','off')
    'Display','iter');          % Available options are 'final','iter' or 'none'

N = 10000;           % number of variables
x0 = randn(N,1);      % initial guess
[x,fval,msg] = liblbfgs(@objective,x0,options);

% Example of syntax for passing additional arguments to the objective function
% [x,fval,msg] = lbfgs(@(x) objective(x,par),x0,options);
