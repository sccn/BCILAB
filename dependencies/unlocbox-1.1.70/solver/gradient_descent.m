function [sol, iter,objectiv] = gradient_descent(x_0,f, param)
%GRADIENT_DESCENT Forward-backward splitting algorithm
%   Usage: sol = gradient_descent(x_0,f, param);
%          sol = gradient_descent(x_0,f);
%          [sol,iter,objective] = gradient_descent(...);
%
%   Input parameters:
%         x_0   : Starting point of the algorithm
%         f    : First function to minimize
%         param : Optional parameter
%   Output parameters:
%         sol   : Solution
%         iter  : Current iteration
%         objectiv: vector (evaluation of the objective function)
%
%   GRADIENT_DESCENT solves:
%
%      sol = argmin f(x)     for x belong to R^N
%
%
%   where x is the variable.
%
%    x_0 is the starting point.
%
%    f is a structure representing a convex function, with a  beta 
%     Lipschitz  continuous gradient. Inside the structure, there
%     have to be the gradient of the function that can be called by f.grad 
%     and the  norm of the function that can be called by f.eval.
%     If no function f is defined, f is by default the L^2 norm.
%
%    param a Matlab structure containing the following fields:
%
%     General parameters:
%
%      param.gamma : is the step size. Watch out, this parameter is bounded. It should
%       be below 1/beta (*f2 is beta Lipchitz continuous). By default, it's 10e-1
%
%      param.tol : is stop criterion for the loop. The algorithm stops if
%
%           (  n(t) - n(t-1) )  / n(t) < tol,
%      
%
%       where  n(t) = f(x) is the objective function at iteration t
%       by default, tol=10e-4.
%
%      param.method : is the method used to solve the problem. It can be 'FISTA' or
%       'ISTA'. By default, it's 'ISTA'. 
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
%   See also:  forward_backward generalized_forward_backward
%
%   Demos: demo_gradient_descent
%
%   References:
%     A. Beck and M. Teboulle. A fast iterative shrinkage-thresholding
%     algorithm for linear inverse problems. SIAM Journal on Imaging
%     Sciences, 2(1):183-202, 2009.
%
%   Url: http://unlocbox.sourceforge.net/doc/solver/gradient_descent.php

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
% Date: 1 nov 2012



% Optional input arguments
if nargin<3, param=struct; end

if ~isfield(param, 'gamma'), param.gamma = 1; end
if ~isfield(param, 'tol'), param.tol=10e-4 ; end
if ~isfield(param, 'maxit'), param.maxit=200; end
if ~isfield(param, 'verbose'), param.verbose=1 ; end
if ~isfield(param, 'lambda'), param.lambda=1 ; end
if ~isfield(param, 'method'), param.method='ISTA' ; end
if ~isfield(param, 'abs_tol'), param.abs_tol=0 ; end

if nargin<2 
    error(' Not enought input argument! Please specify a function to minimise')
end

% test the evaluate function
[f] = test_eval(f);

% Initialization

[~,rel_norm,prev_norm,iter,objectiv,crit] = convergence_test(f.eval(x_0));
    
if param.verbose>=1, 
    if strcmp(param.method,'ISTA'),
        fprintf('Algorithm selected: ISTA\n');
        
    else
        fprintf('Algorithm selected: FISTA\n');

    end
end
% ISTA
x_n = x_0;

%FISTA
u_n=x_0;
sol=x_0;
tn=1;

% Main loop
while 1
    
    %
    if param.verbose >= 1
        fprintf('Iteration %i:\n', iter);
    end
    
    if strcmp(param.method,'ISTA'),
        % ISTA algorithm
        y_n=x_n-param.gamma*f.grad(x_n);
        sol=x_n+param.lambda*(y_n-x_n);
        x_n=sol; % updates
    else
        % FISTA algorithm
        x_n=u_n-param.gamma*f.grad(u_n);
        tn1=(1+sqrt(1+4*tn^2))/2;
        u_n=x_n+(tn-1)/tn1*(x_n-sol);
        %updates
        sol=x_n;
        tn=tn1;
    end
    
    % Global stopping criterion
    curr_norm = f.eval(sol)+eps;  
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
    % Print norm
    fprintf('\n Solution found:\n');
    fprintf(' Final relative norm: %e\n', rel_norm );
    
    
    % Stopping criterion
    fprintf(' %i iterations\n', iter);
    fprintf(' Stopping criterion: %s \n\n', crit);
    
end

end