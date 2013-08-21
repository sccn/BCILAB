function [sol,iter,objectiv] = rlr(x_0,f,A,At, param)
%RLR Regularized Linear Regression 
%   Usage: sol = rlr(x_0,f,A,At, param)
%          sol = rlr(x_0,f,A,At)
%          [sol,iter,objectiv] = rlr(..,)
%
%   Input parameters:
%         x_0   : Starting point of the algorithm
%         f     : Function to minimize
%         A     : Operator
%         At    : Adjoint operator
%         param : Optional parameter
%   Output parameters:
%         sol   : Solution
%         iter  : Current iteration
%         objectiv: vector (evaluation of the objectiv function each iteration)
%
%   This function solve minimization problem using forward-backward splitting
%
%   sol = RLR(x_0,f,A,At, param) solves:
%
%      sol = argmin |X_0-AX|2^2 + f2(x)      for x belong to R^N
%
%
%   where x is the variable. 
%
%    x_0 is the starting point.
%
%    f is a structure representing a convex function. Inside the structure, there
%     have to be the prox of the function that can be called by f.prox and 
%     the  norm of the function that can be called by f.eval. If no function
%     f is defined, f is by default the squared L^1 norm
%
%    A is the operator
%
%    At is the adjoint operator of A
%
%    param a Matlab structure containing the following fields:
%
%     General parameters:
%
%      param.gamma : is the step size. Watch out, this parameter is bounded. It should
%       be below 2/||A||_2. By default, it's 10e-1
%
%      param.tol : is stop criterion for the loop. The algorithm stops if
%
%           (  n(t) - n(t-1) )  / n(t) < tol,
%      
%
%       where  n(t) = f_1(Lx)+f_2(x) is the objective function at iteration t
%       by default, tol=10e-4.
%
%      param.method : is the method used to solve the problem. It can be 'FISTA' or
%       'ISTA'. By default, it's 'FISTA'.
%
%      param.lambda*: is the weight of the update term in ISTA method. By default 1.
%       This should be between 0 and 1. It's set to 1 for FISTA.
%
%      param.maxit : is the maximum number of iteration. By default, it is 200.
% 
%      param.verbose : 0 no log, 1 print main steps, 2 print all steps.
%
%      param.abs_tol : If activated, this stopping critterion is the
%       objectiv function smaller than param.tol. By default 0.
%            
%   See also:  forward_backward douglas_rachford admm
%
%   Demos: demo_rlr 
%
%   References:
%     P. Combettes and J. Pesquet. A douglas-rachford splitting approach to
%     nonsmooth convex variational signal recovery. Selected Topics in Signal
%     Processing, IEEE Journal of, 1(4):564-574, 2007.
%
%   Url: http://unlocbox.sourceforge.net/doc/solver/rlr.php

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

% test the evaluate function
[f] = test_eval(f);

% Author: Nathanael Perraudin, Gilles Puy
% Date: sept 30 2011
%



% Optional input arguments
if nargin<5, param=struct; end

if ~isfield(param, 'gamma'), param.gamma = 1; end
if ~isfield(param, 'tol'), param.tol=10e-4 ; end
if ~isfield(param, 'maxit'), param.maxit=200; end
if ~isfield(param, 'verbose'), param.verbose=1 ; end
if ~isfield(param, 'lambda'), param.lambda=1 ; end
if ~isfield(param, 'method'), param.method='FISTA' ; end


% setting the function f2 
f2.grad=@(x) 2*At((A(x)-x_0));
f2.eval=@(x) (norm(A(x)-x_0,'fro'))^2;

[sol,iter,objectiv]=forward_backward(x_0,f,f2,param);
