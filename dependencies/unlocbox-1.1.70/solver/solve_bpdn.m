function [sol,iter,objectiv] = solve_bpdn(y, epsilon, A, At, Psi, Psit, param)
%SOLVE_BPDN Solve BPDN (basis pursuit denoising) problem
%   Usage: sol = solve_bpdn(y, epsilon, A, At, Psi, Psit, param)
%          sol = solve_bpdn(y, epsilon, A, At, Psi, Psit)
%          [sol,iter,objectiv] = solve_bpdn(...)
%
%   Input parameters:
%         y     : Measurements
%         epsilon: Radius of the L2 ball
%         A     : Operator
%         At    : Adjoint of A
%         Psi   : Operator
%         Psit  : Adjoint of Psi
%         param : Optional parameter
%   Output parameters:
%         sol   : Solution
%         iter  : Current iteration
%         objectiv: vector (evaluation of the objectiv function each iteration)
%
%   sol = solve_BPDN(y, A, At, Psi, Psit, param) solves:
%
%      sol arg min |PSI X|1   s.t.  ||y-A x||_2 < epsilon
%
%
%   Y contains the measurements. A is the forward measurement operator and
%   At the associated adjoint operator. Psit is a sparfying transform and Psi
%   its adjoint. PARAM a Matlab structure containing the following fields:
%
%   General parameters:
% 
%    param.verbose : 0 no log, 1 print main steps, 2 print all steps.
%
%    param.maxit : max. nb. of iterations (default: 200).
%
%    param.tol : is stop criterion for the loop. The algorithm stops if
%
%         (  n(t) - n(t-1) )  / n(t) < tol,
%      
%
%     where  n(t) = ||Psi(x)|| is the objective function at iteration t
%     by default, tol=10e-4.
%
%    param.gamma : control the converge speed (default: 1e-1).
% 
% 
%   Projection onto the L2-ball :
%
%    param.tight_b2 : 1 if A is a tight frame or 0 if not (default = 1)
% 
%    nu_b2 : bound on the norm of the operator A, i.e.
%
%        ` ||A x||^2 <= nu  ||x||^2 
%
%
%    tol_b2 : tolerance for the projection onto the L2 ball (default: 1e-3):
%
%      epsilon/(1-tol) <= |Y - A Z|2 <= epsilon/(1+tol)
%
%    
%    maxit_b2 : max. nb. of iterations for the projection onto the L2
%     ball (default 200).
% 
% 
%   Proximal L1 operator:
%
%    tol_l1 : Used as stopping criterion for the proximal L1
%     operator. Min. relative change of the objective value between two
%     successive estimates.
%
%    maxit_l1 : Used as stopping criterion for the proximal L1
%     operator. Maximum number of iterations.
% 
%    param.nu_l1 : bound on the norm^2 of the operator Psi, i.e.
%
%        ` ||Psi x||^2 <= nu  ||x||^2 
%
% 
%    param.tight_l1 : 1 if Psit is a tight frame or 0 if not (default = 1)
% 
%    param.weights : weights (default = 1) for a weighted L1-norm defined
%     as:
%
%        sum_i{weights_i.*abs(x_i)}
%
%
%   The problem is solved thanks to a Douglas-Rachford splitting
%   algorithm.
%
%   Demos: demo_weighted_l1
%
%   References:
%     P. Combettes and J. Pesquet. A douglas-rachford splitting approach to
%     nonsmooth convex variational signal recovery. Selected Topics in Signal
%     Processing, IEEE Journal of, 1(4):564-574, 2007.
%
%   Url: http://unlocbox.sourceforge.net/doc/solver/solve_bpdn.php

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
%

% Optional input arguments
if nargin<7, param=struct; end

% Optional input arguments
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'tol'), param.tol = 1e-4; end
if ~isfield(param, 'maxit'), param.maxit = 200; end
if ~isfield(param, 'gamma'), param.gamma = 1e-2; end
if ~isfield(param, 'pos_l1'), param.pos_l1 = 0; end

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

% Input arguments for prox L1
param_l1.A = Psi; param_l1.At = Psit; param_l1.pos = param.pos_l1;
param_l1.verbose = param.verbose; param_l1.tol = param.tol;
if isfield(param, 'nu_l1')
    param_l1.nu = param.nu_l1;
end
if isfield(param, 'tight_l1')
    param_l1.tight = param.tight_l1;
end
if isfield(param, 'maxit_l1')
    param_l1.maxit = param.maxit_l1;
end
if isfield(param, 'tol_l1')
    param_l1.tol = param.tol_l1;
end
if isfield(param, 'weights')
    param_l1.weights = param.weights;
else
    param_l1.weights = 1;
end

% Initialization
xhat = At(y); 
[~,~,prev_norm,iter,objectiv,~] = convergence_test();

% Main loop
while 1
    
    %
    if param.verbose >= 1
        fprintf('Iteration %i:\n', iter);
    end
    
    % Projection onto the L2-ball
    [sol, param_b2.u] = fast_proj_b2(xhat, NaN, param_b2);
    
    % Global stopping criterion
    dummy = Psit(sol);
    curr_norm = sum(param_l1.weights(:).*abs(dummy(:)));    
    [stop,rel_norm,prev_norm,iter,objectiv,crit] = convergence_test(curr_norm,prev_norm,iter,objectiv,param);
    if stop
        break;
    end
    if param.verbose >= 1
        fprintf('  ||x||_1 = %e, rel_norm = %e\n', ...
            curr_norm, rel_norm);
    end

    
    % Proximal L1 operator
    xhat = 2*sol - xhat;
    temp = prox_l1(xhat, param.gamma, param_l1);
    xhat = temp + sol - xhat;
    
    
end

% Log
if param.verbose>=1
    % L1 norm
    fprintf('\n Solution found:\n');
    fprintf(' Final L1 norm: %e\n', curr_norm);
    
    % Residual
    dummy = A(sol); res = norm(y(:)-dummy(:), 2);
    fprintf(' epsilon = %e, ||y-Ax||_2=%e\n', epsilon, res);
    
    % Stopping criterion
    fprintf(' %i iterations\n', iter);
    fprintf(' Stopping criterion: %s \n\n', crit);
    
end

end