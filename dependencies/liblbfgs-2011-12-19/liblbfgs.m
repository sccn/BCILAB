function [x,fx,msg] = liblbfgs(f,x0,varargin)
% LIBLBGFS - Limited memory quasi-Newton method for solving
%         problems of the form
%
% min_{x} f(x) + c*|\tilde{x}|,
%
% where x are the parameters, f is a scalar valued real function, c is a
% positive scalar (default value=0), |.| denotes the 1-norm, and \tilde{x}
% is some subset of the parameters x.
%
% Inputs:  1. f - function handle specifying the objective function.
%                 The function must take x as an argument
%                 and return the function value and gradient at x.
%          2. x0 - initial guess
%          3. options (name-value pairs)
%             'm':                  The number of corrections to approximate the inverse hessian matrix 
%                                   (default: 6)
%             'epsilon':            A minimization terminates when ||g|| < \ref epsilon * max(1, ||x||) 
%                                   (default: 1e-5)
%             'past':               Distance for delta-based convergence test (default: 0)
%             'delta':              Delta for convergence test (default: 1e-5)
%             'MaxIter':            The maximum number of iterations (default: 10)
%             'linesearch':         The line search algorithm; can be one of the following:
%                                   - 'more_thuente' (default),
%                                   - 'backtracking_armijo',
%                                   - 'backtracking_wolfe',
%                                   - 'backtracking_strong_wolfe'
%             'max_linesearch':     The maximum number of trials for the line search (default: 40)
%             'min_step':           The minimum step of the line search (default: 1e-20)
%             'max_step':           The minimum step of the line search routine (default: 1e20)
%             'ftol':               A parameter to control the accuracy of the line search routine (default: 1e-4)
%             'wolfe':              A coefficient for the Wolfe condition (default: 0.9)
%             'gtol':               A parameter to control the accuracy of the line search routine (default: 0.9)
%             'xtol':               The machine precision for floating-point values (default: 1e-16)
%             'orthantwise_c':      Coefficient for the L1 norm of variables (default: 0)
%             'orthantwise_start':	First index for the parameters subject to L1-penalty (default: 0)
%             'orthantwise_end':    Last index for the parameters subject to L1-penalty (default: -1)
%             'DerivativeCheck':    Derivative check using finite differences (default: 'off')
%             'Display':            Options for displaying progress (default: 'none')
%
% Outputs: 1. x - solution vector
%          2. fx - objective function value at x
%          3. msg - status message (scalar). See code below for
%                   interpretation.
%
% See "example.m" for an example. References are found in "readme.txt"
%
% NOTE: When using the 1-norm penalization (i.e. c>0), the 'linesearch'-option must be set to
% 'backtracking' or 'backtracking_wolfe' (default value='more_thuente'). Also note that if
% f(x) is not convex, the algorithm used to handle the 1-norm is not guaranteed to converge to a local minima
%

if nargin<2
    error('lbfgs requires at least two input arguments'); end
if ~isa(f,'function_handle')
    error('First argument must be a function handle'); end
if ~isnumeric(x0)
    error('Second argument must be a numeric array'); end

% Parse options
vals = [{ ...
    'm',6,...                   % The number of corrections to approximate the inverse hessian matrix.
    'epsilon',1e-5,...          % A minimization terminates when ||g|| < \ref epsilon * max(1, ||x||)
    'past',0,...                % Distance for delta-based convergence test
    'delta',1e-5,...            % Delta for convergence test.
    'MaxIter',10,...            % The maximum number of iterations.
    'linesearch','more_thuente',...  % The line search algorithm.
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
    'DerivativeCheck','off',... % Derivative check using finite differences
    'Display','none'...         % Options for displaying progress
    } varargin];

if ~isempty(varargin)
    % get names & values
    names = vals(1:2:end);
    values = vals(2:2:end);
    % use only the last assignment for each name
    [s,indices] = sort(names);
    indices(strcmp(s(1:end-1),s(2:end))) = [];
    % build a struct
    defaults = cell2struct(values(indices),names(indices),2);
else
    defaults = cell2struct(vals(2:2:end),vals(1:2:end),2);
end

% Set linesearch algorithm
if defaults.orthantwise_c
    switch defaults.linesearch
        case {'backtracking','backtracking_wolfe'}
            defaults.linesearch = 2;
        otherwise
            error(['options.linesearch must be set to "backtracking"' ...
                'or "backtracking_wolfe" when using the orthantwise method.']);
    end
else
    switch defaults.linesearch
        case {'more_thuente','default'}
            defaults.linesearch = 0;
        case 'backtracking_armijo'
            defaults.linesearch = 1;
        case {'backtracking','backtracking_wolfe'}
            defaults.linesearch = 2;
        case 'backtracking_strong_wolfe'
            defaults.linesearch = 3;
        otherwise
            error('Unknown line search algorithm specified. Available options are ''more_thuente'', ''backtracking_armijo'', ''backtracking'', ''backtracking_wolfe'' or ''backtracking_strong_wolfe''.');
    end
end

defaults.Display = lower(defaults.Display);

% Correct for 1-based indexing in matlab and 0-based indexing in c
if defaults.orthantwise_start>0
    defaults.orthantwise_start = defaults.orthantwise_start-1;
end
x0 = x0(:);

% Check derivatives with finite differences
if strcmpi(defaults.DerivativeCheck,'on')
    [unused,grad] = f(x0);
    grad_fd = finite_differences(f,x0);
    disp('Analytical gradient:  Finite difference gradient:');
    if numel(grad)<30
        fprintf('%14.8f             %1.8f\n',[grad grad_fd]');
    else
        fprintf('%14.8f             %1.8f\n',[grad(1:30) grad_fd(1:30)]');
        fprintf('(Displaying the first 30 components only)\n\n');
    end
    [err,ind] = max(abs(grad-grad_fd));
    fprintf('\nMax difference %1.8f in component %d.\n\n',err,ind);
    fprintf('Press any key to continue.\n\n')
    pause
end

% Call mex function
[x,fx,msg] = lbfgs_(f,x0,defaults);

if ~strcmp(defaults.Display,'none')
    % Interpret and display exit message
    switch msg
        case {0,1}
            disp('Success!');
        case 2
            disp('Initial variables already minimze the objective function.');
        case -1024
            disp('Unknown error.');
        case -1023
            disp('Logic error.');
        case -1022
            disp('Out of memory');
        case -1021
            disp('The minimization process has been canceled.');
        case -1020
            disp('Invalid number of variables specified.');
        case -1019
            disp('Invalid number of variables (for SSE) specified.');
        case -1018
            disp('The array x must be aligned to 16 (for SSE).');
        case -1017
            disp('Invalid parameter epsilon specified.');
        case -1016
            disp('Invalid parameter past specified.');
        case -1015
            disp('Invalid parameter delta specified.');
        case -1014
            disp('Invalid parameter linesearch specified.');
        case -1013
            disp('Invalid parameter min_step specified.');
        case -1012
            disp('Invalid parameter max_step specified.');
        case -1011
            disp('Invalid parameter ftol specified.');
        case -1010
            disp('Invalid parameter wolfe specified.');
        case -1009
            disp('Invalid parameter gtol specified.');
        case -1008
            disp('Invalid parameter xtol specified.');
        case -1007
            disp('Invalid parameter max_linesearch specified.');
        case -1006
            disp('Invalid parameter orthantwise_c specified.');
        case -1005
            disp('Invalid parameter orthantwise_start specified.');
        case -1004
            disp('Invalid parameter orthantwise_end specified.');
        case -1003
            disp('The line-search step went out of the interval of uncertainty.');
        case -1002
            disp('A logic error occurred; alternatively, the interval of uncertainty became too small.');
        case -1001
            disp('A rounding error occurred; alternatively, no line-search step satisfies the sufficient decrease and curvature conditions.');
        case -1000
            disp('The line-search step became smaller than min_step.');
        case -999
            disp('The line-search step became larger than max_step.');
        case -998
            disp('The line-search routine reaches the maximum number of evaluations.');
        case -997
            disp(['The maximum number of iterations (' num2str(options.MaxIter) ') was reached.']);
        case -996
            disp('Relative width of the interval of uncertainty is at most xtol.');
        case -995
            disp('A logic error (negative line-search step) occurred.');
        case -994
            disp('The current search direction increases the objective function value.');
        otherwise
            disp(['Unknown exit code: ' msg]);
    end
end
end


% Finite difference calculation of gradient using forward differences
function grad = finite_differences(func,x)

N = numel(x);
h = sqrt(eps)*nonzerosign(x).*max(abs(x));
f0 = func(x);
grad = zeros(N,1);
for k = 1:N
    x_temp = x(k);
    x(k) = x(k)+h(k);
    grad(k,1) = (func(x)-f0)/h(k);
    x(k) = x_temp;
end

end
