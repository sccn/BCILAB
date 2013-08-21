function [sol, iter,objectiv] = admm(x_0,f1, f2, L, param)
%ADMM alternating-direction method of multipliers
%   Usage: sol = admm(x_0,f1,f2,L,param);
%          sol = admm(x_0,f1,f2,L);
%          [sol,iter,objectiv] = admm(...);
%
%   Input parameters:
%         x_0   : Starting point of the algorithm
%         f1    : First function to minimize
%         f2    : Second function to minimize
%         L     : Linear operation on x
%         param : Optional parameter
%   Output parameters:
%         sol   : Solution
%         iter  : Current iteration
%         objectiv: vector (evaluation of the objectiv function each iteration)
%
%   ADMM (using alternating-direction method of multipliers) solves:
%
%      sol = argmin f1(x) + f2(y) such that y=Lx
%
%   
%   where
%   x is the variable.
%
%    x_0 is the starting point.
%
%    f1 is a structure representing a convex function. Inside the structure, there
%     have to be the prox of the function that can be called by f1.prox and 
%     the  norm of the function that can be called by f1.eval. 
%     WARNING !!!  The prox of f1 is not the usual prox! But the solution to this problem:
%
%        prox_{f1, gamma }^L(z)=min_x  1/2 |LX-Z|2^2 + gamma f1(x)
%
%
%    f2 is a structure representing a convex function. Inside the structure, there
%     have to be the prox of the function that can be called by f2.prox and 
%     the  norm of the function that can be called by f2.eval. If no function
%     f2 is defined, f2 is by default the squared L^1 norm.
%     The prox of f2 is the usual prox:
%
%        prox_{f2, gamma }(z)=min_x  1/2 |X-Z|2^2 + gamma f2(x)
%
%
%    L is a linear operator or a matrix (be careful with operator, they might not be frame)
%           
%
%    param a Matlab structure containing the following fields:
%
%     General parameters:
%
%      param.gamma : is the convergence parameter. By default, it's 1. (greater than 0)  
%                    % this is the augmented Lagrangian parameter (rho in Boyd)
%
%      param.tol : is stop criterion for the loop. The algorithm stops if
%
%           (  n(t) - n(t-1) )  / n(t) < tol,
%      
%
%       where  n(t) = f_1(Lx)+f_2(x) is the objective function at iteration t
%       by default, tol=10e-4.
%
%      param.maxit : is the maximum number of iteration. By default, it is 200.
% 
%      param.verbose : 0 no log, 1 print main steps, 2 print all steps.
%
%   See also:  sdmm, ppxa, generalized_forward_backward
%
%   Demos:  demo_admm
%
%   References:
%     P. Combettes and J. Pesquet. Proximal splitting methods in signal
%     processing. Fixed-Point Algorithms for Inverse Problems in Science and
%     Engineering, pages 185-212, 2011.
%
%   Url: http://unlocbox.sourceforge.net/doc/solver/admm.php

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
% Date: 22 october 2012
%




% Optional input arguments
if nargin<4, param=struct; end

if ~isfield(param, 'tol'), param.tol=10e-4 ; end
if ~isfield(param, 'maxit'), param.maxit=200; end
if ~isfield(param, 'verbose'), param.verbose=1 ; end
if ~isfield(param, 'lambda'), param.lambda=1 ; end     % ck: unused here
if ~isfield(param, 'gamma'), param.gamma=1 ; end

if nargin<3 
    f2.prox=@(x) prox_L1(x, 1, param);                 % ck: bug 
    f2.eval=@(x) norm(x(:),1);   
end

% test the evaluate function
[f1] = test_eval(f1);
[f2] = test_eval(f2);

if isa(L,'numeric')
   OpL= @(x) L*x;
else
   OpL= L;
end

% Initialization

curr_norm = f1.eval(x_0)+f2.eval(OpL(x_0));  
[~,~,prev_norm,iter,objectiv,~] = convergence_test(curr_norm);
y_n = x_0;
z_n = zeros(size(x_0));


% Main loop
while 1
    
    %
    if param.verbose >= 1
        fprintf('Iteration %i:\n', iter);
    end
    
    % note: this does not implement over-relaxation
    
    % Algorithm
    
    x_n=f1.prox(y_n-z_n,param.gamma);  % ck: (x_n is x in Boyd; params.gamma is rho in Boyd)
    s_n=OpL(x_n);                      % ck: (s_n is x_hat in Boyd)
    y_n=f2.prox(s_n+z_n,param.gamma);  % ck: (y_n is z in Boyd)
    z_n=z_n+s_n-y_n ;                  % ck: (z_n is u in Boyd)
    sol=x_n; 

    
    % Global stopping criterion
    curr_norm = f1.eval(sol)+f2.eval(OpL(sol));  
    [stop,rel_norm,prev_norm,iter,objectiv,crit] = convergence_test(curr_norm,prev_norm,iter,objectiv,param);
    if stop
        break;
    end
    if param.verbose >= 1
        fprintf('  ||f|| = %e, rel_norm = %e\n', ...
            curr_norm, rel_norm);
    end  
    
end

% Log
if param.verbose>=1
    % L1 norm
    fprintf('\n Solution found:\n');
    fprintf(' Final relative norm: %e\n', rel_norm );
    
    
    % Stopping criterion
    fprintf(' %i iterations\n', iter);
    fprintf(' Stopping criterion: %s \n\n', crit);
    
end

end
