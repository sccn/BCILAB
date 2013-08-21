function [sol, iter, objectiv] = solve_tvdn(y, epsilon, A, At, param)
%SOLVE_TVDN Solve TVDN problem
%   Usage: sol = solve_tvdn(y, epsilon, A, At, param)
%          sol = solve_tvdn(y, epsilon, A, At)
%          [sol,iter,objectiv] = solve_tvdn(...)
%
%   Input parameters:
%         y     : Measurements
%         epsilon: Radius of the L2 ball
%         A     : Operator
%         At    : Adjoint of A
%         param : Optional parameter
%   Output parameters:
%         sol   : Solution
%         iter  : Current iteration
%         objectiv: vector (evaluation of the objectiv function each iteration)
%
%   sol = SOLVE_TVDN(Y, epsilon, A, At, PARAM) solves:
%
%      sol arg min |X|TV  s.t.  ||y-A x||_2 < epsilon
%
%
%   Y contains the measurements. A is the forward measurement operator and
%   At the associated adjoint operator. PARAM a Matlab structure containing
%   the following fields:
%
%   General parameters:
% 
%    param.verbose : 0 no log, 1 print main steps, 2 print all steps.
%
%    param.maxit : max. nb. of iterations (default: 200).
%
%    param.useGPU : Use GPU to compute the TV prox operator. Please prior 
%     call init_gpu and free_gpu to launch and release the GPU library (default: 0).
%
%    param.tol : is stop criterion for the loop. The algorithm stops if
%
%         (  n(t) - n(t-1) )  / n(t) < tol,
%      
%
%     where  n(t) = |(X)|{TV} is the objective function at iteration t
%     by default, tol=10e-4.
%
%    param.gamma : control the converge speed (default: 1e-1).
% 
% 
%   Projection onto the L2-ball :
%
%    param.tight_b2 : 1 if A is a tight frame or 0 if not (default = 1)
% 
%    param.nu_b2 : bound on the norm of the operator A, i.e.
%
%        ` ||A x||^2 <= nu  ||x||^2 
%
%
%    param.tol_b2 : tolerance for the projection onto the L2 ball (default: 1e-3):
%
%        epsilon/(1-tol) <= |Y - A Z|2 <= epsilon/(1+tol)
%
%    
%    param.maxit_b2 : max. nb. of iterations for the projection onto the L2
%     ball (default 200).
% 
% 
%   Proximal TV operator:
%
%    param.maxit_tv : Used as stopping criterion for the proximal TV
%     operator. Maximum number of iterations.
%
%   The problem is solved thanks to a Douglas-Rachford splitting
%   algorithm.
%
%   Demos: demo_tvdn
%
%   References:
%     P. Combettes and J. Pesquet. A douglas-rachford splitting approach to
%     nonsmooth convex variational signal recovery. Selected Topics in Signal
%     Processing, IEEE Journal of, 1(4):564-574, 2007.
%
%   Url: http://unlocbox.sourceforge.net/doc/solver/solve_tvdn.php

% Copyright (C) 2012 LTS2-EPFL, by Nathanael Perraudin, Gilles Puy,
% David Shuman, Pierre Vandergheynst.
% This file is part of UnLocBoX version 1.1.70
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% Author: Gilles Puy, Nathanael Perraudin
% Date: Nov. 1, 2012


% Optional input arguments
if nargin<5, param=struct; end


% Optional input arguments
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'tol'), param.tol = 1e-4; end
if ~isfield(param, 'maxit'), param.maxit = 200; end
if ~isfield(param, 'gamma'), param.gamma = 1e-2; end
if ~isfield(param, 'useGPU'), param.useGPU = 0; end

% Input arguments for projection onto the L2 ball
param_b2.A = A; param_b2.At = At;
param_b2.y = y; param_b2.epsilon = epsilon;
param_b2.verbose = param.verbose;
if isfield(param, 'nu_b2'), param_b2.nu = param.nu_b2; end
if isfield(param, 'tol_b2'), param_b2.tol = param.tol_b2; end
if isfield(param, 'tight_b2'), param_b2.tight = param.tight_b2; end
if isfield(param, 'maxit_b2')
    param_b2.maxit = param.maxit_b2;
end

% Input arguments for prox TV
param_tv.verbose = param.verbose; param_tv.tol = param.tol;
param_tv.useGPU = param.useGPU;
if isfield(param, 'maxit_tv')
    param_tv.maxit= param.maxit_tv;
end

% Initialization
xhat = At(y); 
[~,~,prev_norm,iter,objectiv,~] = convergence_test(tv_norm(y));
% Main loop
while 1
    
    %
    if param.verbose>=1
        fprintf('Iteration %i:\n', iter);
    end
    
    % Projection onto the L2-ball
    [sol, param_b2.u] = fast_proj_b2(xhat, NaN, param_b2);
    
    % Global stopping criterion
    curr_norm = tv_norm(sol);
    [stop,rel_norm,prev_norm,iter,objectiv,crit] = convergence_test(curr_norm,prev_norm,iter,objectiv,param);
    if stop
        break;
    end
    if param.verbose >= 1
        fprintf('  ||x||_TV = %e, rel_norm = %e\n', ...
            curr_norm, rel_norm);
    end

    
    % Proximal L1 operator
    xhat = 2*sol - xhat;
    temp = prox_tv(xhat, param.gamma, param_tv);
    xhat = temp + sol - xhat;
    
    
end

% Log
if param.verbose>=1
    % L1 norm
    fprintf('\n Solution found:\n');
    fprintf(' Final TV norm: %e\n', curr_norm);
    
    % Residual
    temp = A(sol);
    fprintf(' epsilon = %e, ||y-Ax||_2=%e\n', epsilon, ...
        norm(y(:)-temp(:)));
    
    % Stopping criterion
    fprintf(' %i iterations\n', iter);
    fprintf(' Stopping criterion: %s \n\n', crit);
    
end

end