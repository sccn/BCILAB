function [sol, iter,objectiv] = generalized_forward_backward(x_0, F, f , param)
%GENERALIZED_FORWARD_BACKWARD Generalized forward backward algorithm
%   Usage: sol = generalized_forward_backward(x_0,F, f2, param);
%          sol = generalized_forward_backward(x_0,F, f2);
%          [sol,iter,objectiv] = generalized_forward_backward(...);
%
%   Input parameters:
%         x_0   : Starting point of the algorithm
%         F     : Array of structure representing the functions to minimize
%         f2    : Another function to minimize with a known gradient
%         param : Optional parameter
%   Output parameters:
%         sol   : Solution
%         iter  : Current iteration
%         objectiv: vector (evaluation of the objectiv function each iteration)
% 
%   GENERALIZED_FORWARD_BACKWARD solves:
%
%      sol = argmin_{z} f2(z) + Sum_{i} wi Fi(z)     for z belong to R^N
%
%
%   With z the variable and wi the weight accorded to every term of the sum
%
%    x_0 : is the starting point.
%
%    F : a structure array with all the prox operator inside and eventually 
%     the norm if no norm is defined, the L^1 norm is used the prox: F(i).prox and 
%     norm: F(i).eval are defined in the same way as in the
%     Forward-backward and Douglas-Rachford algorithms
%
%    f2 is a structure representing a convex function, with a  beta 
%     Lipschitz  continuous gradient. Inside the structure, there
%     have to be the gradient of the function that can be called by f2.grad
%     and the  norm of the function that can be called by f2.eval.
%     If no function f2 is defined, f2 is by default the L^2 norm.
%
%    param is a Matlab structure containing the following fields:
%
%     General parameters:
%
%      param.gamma : is the step size. Watch out, this parameter is bounded. It should
%       be below 1/beta (*f2 is beta Lipschitz continuous). By default, it's 10e-1
%
%      param.weights : weights of different functions (default = 1/N,
%        where N is the total number of function) 
%
%      param.tol : is stop criterion for the loop. The algorithm stops if
%
%           (  n(t) - n(t-1) )  / n(t) < tol,
%      
%
%       where  n(t) = f_1(Lx)+f_2(x) is the objective function at iteration t
%       by default, tol=10e-4.
%
%      param.abs_tol : If activated, this stopping critterion is the
%       objectiv function smaller than param.tol. By default 0.
%
%      param.lambda*: is the weight of the update term. By default 1.
%       This should be between 0 and 1.
%
%      param.maxit : is the maximum number of iteration. By default, it is 200.
% 
%      param.verbose : 0 no log, 1 print main steps, 2 print all steps.     
%
%   See also:  douglas_rachford ppxa admm
%
%   Demos: demo_generalized_forward_backward
%
%   References:
%     H. Raguet, J. Fadili, and G. Peyré. Generalized forward-backward
%     splitting. arXiv preprint arXiv:1108.4404, 2011.
%     
%
%   Url: http://unlocbox.sourceforge.net/doc/solver/generalized_forward_backward.php

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
% Date: Oct 24 2012
%

% % Reshape x if vector line
% if (size(x,2)>size(x,1))
%     x=x';
% end

% Optional input arguments
if nargin<4, param=struct; end

if ~isfield(param, 'weights'), param.weights=ones(size(F,2),1); end
if ~isfield(param, 'verbose'), param.verbose=1 ; end
if ~isfield(param, 'maxit'), param.maxit=100 ; end
if ~isfield(param, 'tol'), param.tol=1e-3 ; end
if ~isfield(param, 'gamma'), param.gamma=1 ; end
if ~isfield(param, 'lambda'), param.lambda=1 ; end

if nargin<3 
    f.grad=@(x) 2*x;
    f.eval=@(x) norm(x(:)-x_0(:),2)^2;  
end


% Normalizing the weights
param.weights= param.weights/sum(param.weights);


% Number of function
l=size(F,2);



% test the evaluate function
[f] = test_eval(f);
[F] = test_eval(F);


% Definition of the gradient function
grad = @(y) f.grad(y);

% Algorithm - Initialisation
z=zeros([l,size(x_0)]);


for i=1:l
    z(i,:,:)=x_0;
end


sol=x_0;
curr_norm = f.eval(sol)+norm_sumg(sol,F);
[~,~,prev_norm,iter,objectiv,~] = convergence_test(curr_norm);

% Algorithm - Loop

while 1
    
    %
    if param.verbose >= 1
        fprintf('Iteration %i:\n', iter);
    end
    i=1;
    for g=F
       temp=2*sol-reshape(z(i,:,:),size(x_0))-param.lambda*grad(sol);
       z(i,:,:) = reshape(z(i,:,:),size(x_0))+ param.gamma*(g.prox(temp,1/param.weights(i)*param.lambda)-sol);
       i=i+1;
    end
    
    sol=zeros(size(x_0));
    for i=1:l
        sol=sol+param.weights(i)* reshape(z(i,:,:),size(x_0));
    end
    
    % Global stopping criterion

    curr_norm = f.eval(sol)+norm_sumg(sol,F);
    [stop,rel_norm,prev_norm,iter,objectiv,crit] = convergence_test(curr_norm,prev_norm,iter,objectiv,param);
    if stop
        break;
    end
    if param.verbose >= 1
        fprintf(' Norm of the general objectiv function: ||f|| = %e, rel_norm = %e\n', ...
            curr_norm, rel_norm);
    end


  
    
end





% Calculation of the norm
norm_G=f.eval(sol)+norm_sumg(sol,F);

% Log after the calculous of the prox
if param.verbose >= 1
    fprintf('  Generalized forward backward: Sum_k ||G_k(x)|| = %e\n', norm_G);


    % Stopping criterion
    fprintf(' %i iterations\n', iter);
    fprintf(' Stopping criterion: %s \n\n', crit);

end
end