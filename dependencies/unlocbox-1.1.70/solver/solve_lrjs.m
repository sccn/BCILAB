function [sol, iter, objectiv] = solve_lrjs(x0,y, epsilon, A, At, param)
%SOLVE_LRJS Solve LRJS (low rank joint sparsity) problem
%   Usage: sol = solve_lrjs(y, epsilon, A, At, param)
%          sol = solve_lrjs(y, epsilon, A, At)
%          [sol,iter,objectiv] = solve_lrjs(...)
%
%   Input parameters:
%         x0    : a starting point (matrix of the good size)
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
%   sol = solve_tvdn(Y, epsilon, A, At, PARAM) solves:
%
%      sol arg min |X| + alpha*||x||_12  s.t.  ||y-A x||_2 < epsilon
%
%
%   Y contains the measurements. A is the forward measurement operator and
%   At the associated adjoint operator. PARAM a Matlab structure containing
%   the following fields:
%
%   General parameters:
% 
%    param.alpha : regularization parameter (default: 1).
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
%     where  n(t) = |(X)|{TV} is the objective function at iteration t
%     by default, tol=10e-4.
%
%    param.gamma : control the converge speed (default: 1).
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
%   Proximal L!2 operator:
%
%    param.g_d , param.g_t are the group vectors. 
%     (default param.g_d=1:length(x);, param.g_t=ones(size(x));*)
%
%     param.g_d contains the indices of the elements to be grouped and param.g_t the size of the different groups.
%
%     Warning: param.g_d and param.g_t have to be row vector!     
%     
%     Example: suppose x=[x1 x2 x3 x4 x5 x6] and Group 1: [x1 x2 x4 x5] group 2: [x3 x6]
%              
%     In matlab: 
%
%           param.g_d=[1 2 4 5 3 6]; param.g_t=[4 2];
%
%     Or this is also possible:
%
%           param.g_d=[4 5 3 6 1 2]; param.g_t=[2 4]; 
%
%    param.maxit_tv : Used as stopping criterion for the proximal TV
%     operator. Maximum number of iterations.
%
%   The problem is solved thanks to the ppxa solver
%
%   Demos: demo_lrjs
%
%   References:
%     M. Golbabaee and P. Vandergheynst. Compressed sensing of simultaneous
%     low-rank and joint-sparse matrices. In Infoscience - EPFL, 2012.
%     
%
%   Url: http://unlocbox.sourceforge.net/doc/solver/solve_lrjs.php

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

% Author: Nathanael Perraudin
% Date: Jan. 21, 2013


% Optional input arguments
if nargin<6, param=struct; end

% get size parameter
[n1,n2]=size(x0);

% Optional input arguments
if ~isfield(param, 'verbose'), param.verbose = 2; end
if ~isfield(param, 'tol'), param.tol = 1e-4; end
if ~isfield(param, 'maxit'), param.maxit = 200; end
if ~isfield(param, 'gamma'), param.gamma = 1; end
if ~isfield(param, 'alpha'), param.alpha = 1; end

% Input arguments for projection onto the L2 ball
param_b2.A = A; param_b2.At = At;
param_b2.y = y; param_b2.epsilon = epsilon;
param_b2.verbose = param.verbose-1;
if isfield(param, 'nu_b2'), param_b2.nu = param.nu_b2; end
if isfield(param, 'tol_b2'), param_b2.tol = param.tol_b2; end
if isfield(param, 'tight_b2'), param_b2.tight = param.tight_b2; end
if isfield(param, 'maxit_b2')
    param_b2.maxit = param.maxit_b2;
end

% Input arguments for prox L12
if isfield(param, 'g_d'), param_l12.g_d = param.g_d; end
if isfield(param, 'g_t'), param_l12.g_t = param.g_t; end
param_l12.verbose=param.verbose-1;




% projection  
F(1).prox=@(x,T) proj_b2_test(x, y, A, At, param_b2);
F(1).eval=@(x) 0;
    
    
%prox l12
F(2).prox=@(x,T) prox_l12(x,param.alpha*T,param_l12);
F(2).eval=@(x) param.alpha*norm_l12(x,param_l12.g_d,param_l12.g_t);
    
    
%prox nuclearnorm
param_svt.verbose=param.verbose-1;
F(3).prox=@(x,T) reshape(prox_nuclearnorm(reshape(x,n1,n2),T,param_svt),n1*n2,1);
F(3).eval=@(x) sum((svd(reshape(x,n1,n2))));
    

  
%solve the problem
[sol, iter, objectiv]=ppxa(x0(:),F,param);
sol=reshape(sol,n1,n2);

end