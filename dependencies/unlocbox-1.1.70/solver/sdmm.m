function [x_n,iter,hist] = sdmm(F, param)
%SDMM Simultaneous-direction method of multipliers algorithm
%   Usage: sol = sdmm(F,param);
%          sol = sdmm(F);
%          [sol,iter,objectiv] = sdmm(...);
%
%   Input parameters:
%         F     : Array of function to minimize
%         param : Optional parameter
%   Output parameters:
%         sol   : Solution
%         iter  : Current iteration
%         objectiv: vector (evaluation of the objectiv function each iteration)
%
%   SDMM, from simultaneous-direction method of multipliers solves:
%
%      sol = argmin sum(f_i( L_i x))
%
%
%   where x belong to R^N, L_i are linear operators and x_i are the
%   minimization variables.
%
%    F is an array of structure representing the functions.
%
%     In the function F(i), there have to be:
%
%      F(i).eval(x_i) : an operator to evaluate the function
%      F(i).prox(x_i, gamma) : an operator to evaluate the prox of the function
%      F(i).x0 : vector of initial value
%
%     Optionally you can define
%
%      F(i).L : linear operator, given in a matrix form or an operator (default identity)
%
%    param is a Matlab structure containing the following fields:
%
%     General parameters:
%
%      param.tol : is stop criterion for the loop. The algorithm stops if
%
%           (  n(t) - n(t-1) )  / n(t) < tol,
%
%
%       where  n(t) = sum W_i*f_i(x) is the objective function at iteration t*
%       by default, tol=10e-4.
%
%      param.maxit : is the maximum number of iteration. By default, it is 200.
%
%      param.verbose : 0 no log, 1 print main steps, 2 print all steps.
%
%      param.gamma : convergence parameter (default 1)
%
%      param.abs_tol : If activated, this stopping critterion is the
%       objectiv function smaller than param.tol. By default 0.
%
%   See also:  admm forward_backward douglas_rachford
%
%   Demos:  demo_sdmm
%
%   References:
%     P. Combettes and J. Pesquet. A douglas-rachford splitting approach to
%     nonsmooth convex variational signal recovery. Selected Topics in Signal
%     Processing, IEEE Journal of, 1(4):564-574, 2007.
%
%   Url: http://unlocbox.sourceforge.net/doc/solver/sdmm.php

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


% Author:  Nathanael Perraudin
% Date: fev 23 2012
%

% Relationship to ADMM:
%
%       minimize sum  f_i(A_i*x)
%           s.t. A_i*x_i - z = 0
%
% L_i = A_i                                 L_i is the A_i in ADMM in a generalized consensus form
% F(i).y_n = A_i*x_i                        y_n's are in a space transformed by A
% f in F(i).prox = f(A_i*x)                 f's are implicitly applied to A-transformed x's
% F(i).z_n = -A_i*u_i = -A_i*(1/rho)*u_i    z_n's are negated and scaled dual variables after transformation into A
% x_n = z                                   x_n is the consensus variable
%
% B = -eye                          this comes from the Ax - z = 0
% c = 0

% handle optional arguments
if nargin<2
    param=struct(); end
if ~isfield(param,'tol')
    param.tol = 10e-4; end
if ~isfield(param,'maxit')
    param.maxit = 200; end
if ~isfield(param,'verbose')
    param.verbose = 1; end
if ~isfield(param,'gamma')
    param.gamma = 1; end

% sanity-check the functions
for i=1:length(F)
    if ~isfield(F(i),'L') || isempty(F(i).L)
        F(i).L = eye(numel(F(i).x0)); end
    F(i).x0 = F(i).x0(:);
    F(i).y_n=F(i).x0;                                           % this is utter crap & must be fixed...
    F(i).z_n=F(i).x0;
end

% initialize the (most likely diagonal) rescaling matrix Q -- because L'*L must be diagonal for them to be "tight frames"
Q = 0;
for i=1:m
    Q = Q + F(i).L' * F(i).L; end
Q_inv = inv(Q);

% initialize the consensus variable z from the primal and scaled dual variables A_i*x_i and B_i*u_i
x_n = 0;
for i=1:m
    x_n = x_n + F(i).L' * (F(i).y_n - F(i).z_n); end
x_n = Q_inv * x_n;

% Main loop
[~,~,prev_obj,iter,hist,~] = convergence_test(objective(x_n,F));
while 1
    if param.verbose >= 1
        fprintf('Iteration %i:\n', iter); end
    
    for i=1:m
        % primal variable update
        s_n = F(i).L*x_n;                                       % Ai*z
        F(i).y_n = F(i).prox(s_n + F(i).z_n,param.gamma);       % Ai*xi = argmin_{x} 0.5*||Ai*z + (-Ai*ui) - Ai*x||_2^2 + gamma*f(Ai*x) <=> Ai * (argmin_{x} 0.5*||x - z + ui||_2^2 + gamma*f(x))
                                                                % the prox operator works on variables that have been mapped by Ai into the appropriate space
        % dual variable update
        F(i).z_n = F(i).z_n + s_n - F(i).y_n;                   % -Ai*ui = -Ai*ui + Ai*z - Ai*xi <=> ui = ui + xi - z
    end
    
    % consensus variable update
    x_n = 0;                                                    
    for i=1:length(F)
        x_n = x_n + F(i).L'*(F(i).y_n-F(i).z_n); end            % z = 1/Q * Ai'(Ai*xi - (-Ai*ui)) = 1/Q * Ai'(Ai*xi + Ai*ui)
    x_n=Q_inv*x_n;
    
    % check stopping conditions
    cur_obj = objective(x_n,F);
    [stop,rel_norm,prev_obj,iter,hist,crit] = convergence_test(cur_obj,prev_obj,iter,hist,param);
    if stop
        break; end
    
    if param.verbose >= 1
        fprintf('  ||f|| = %e, rel_norm = %e\n', curr_obj, rel_norm); end
end

function obj = objective(x_n,F)
obj = 0;
for i=1:length(F)
    obj = obj + F(i).eval(F(i).L*x_n); end
