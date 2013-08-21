function [sol,iter] = prox_l2(x, gamma, param)
%PROX_L2 Proximal operator with L2 norm
%   Usage:  sol=prox_l2(x, gamma)
%           sol=prox_l2(x, gamma, param)
%           [sol, iter]=prox_l2(x, gamma, param)
%
%   Input parameters:
%         x     : Input signal.
%         gamma : Regularization parameter.
%         param : Structure of optional parameters.
%   Output parameters:
%         sol   : Solution.
%         iter  : Number of iterations at convergence.
%
%   PROX_L2(x, gamma, param) solves:
%
%      sol = argmin_{z} 0.5*||x - z||_2^2 + gamma  |A W Z - Y|2^2
%
%
%   where w are some weights.
%
%   param is a Matlab structure containing the following fields:
%
%    param.weights : weights for a weighted L2-norm (default = 1)
%
%    param.y : measurements (default: 0).
%
%    param.A : Forward operator (default: Id).
%
%    param.At : Adjoint operator (default: Id).
%
%    param.tight : 1 if A is a tight frame or 0 if not (default = 1)
%
%    param.nu : bound on the norm of the operator A (default: 1), i.e.
%
%        ` ||A x||^2 <= nu  ||x||^2 
%
%   
%    param.tol : is stop criterion for the loop. The algorithm stops if
%
%         (  n(t) - n(t-1) )  / n(t) < tol,
%      
%
%     where  n(t) = f(x)+ 0.5 X-Z2^2 is the objective function at iteration t
%     by default, tol=10e-4.
%
%    param.maxit : max. nb. of iterations (default: 200).
%
%    param.verbose : 0 no log, 1 a summary at convergence, 2 print main
%     steps (default: 1)
%
%   See also:  proj_b2 prox_l1
%
%   Url: http://unlocbox.sourceforge.net/doc/prox/prox_l2.php

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
% Date: Nov 2012
%

% Optional input arguments
if nargin<3, param=struct; end

% Optional input arguments
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'tight'), param.tight = 1; end
if ~isfield(param, 'nu'), param.nu = 1; end
if ~isfield(param, 'tol'), param.tol = 1e-3; end
if ~isfield(param, 'maxit'), param.maxit = 200; end
if ~isfield(param, 'At'), param.At = @(x) x; end
if ~isfield(param, 'A'), param.A = @(x) x; end
if ~isfield(param, 'y'), param.y = zeros(size(x)); end
if ~isfield(param, 'weights'), param.weights = ones(size(x)); end

siz = size(x);

% test the parameters
gamma=test_gamma(gamma);
param.weights=test_weights(param.weights);

% Projection
if param.tight % TIGHT FRAME CASE
    
    sol=(x+gamma*2*param.At(param.y).*param.weights)./(gamma*2*param.nu*param.weights.^2+1);
    
    
    % Infos for log...
    % L2 norm of the estimate
    dummy = param.A(param.weights.*sol)-param.y;
    norm_l2 = .5*norm(x(:) - sol(:), 2)^2 + gamma * norm(dummy(:))^2;
    % stopping criterion
    crit_L2 = 'REL_OB'; 
    % number of iteration
    iter_L2=0;
else % NON TIGHT FRAME
    
    % Initializations
    u_n=x;
    sol=x;
    tn=1;
    prev_l2 = 0; iter_L2 = 0;
    % stepsize
    stepsize=1/(2*max(param.weights).^2*param.nu+1);
    
    % gradient
    grad= @(z) z-x+gamma*2*param.weights.*param.At(param.A(param.weights.*z)-param.y);
    
    % Init
    if param.verbose > 1
        fprintf('  Proximal l2 operator:\n');
    end
    while 1
        
        % L2 norm of the estimate
        dummy = param.A(param.weights.*sol)-param.y;
        norm_l2 = .5*norm(x(:) - sol(:), 2)^2 + gamma * norm(dummy(:))^2;
        rel_l2 = abs(norm_l2-prev_l2)/norm_l2;
        
        % Log
        if param.verbose>1
            fprintf('   Iter %i, ||A w x- y||_2^2 = %e, rel_l2 = %e\n', ...
                iter_L2, norm_l2, rel_l2);
        end
        
        % Stopping criterion
        if (rel_l2 < param.tol)
            crit_L2 = 'REL_OB'; break;
        elseif iter_L2 >= param.maxit
            crit_L2 = 'MAX_IT'; break;
        end
        
        % FISTA algorithm
        x_n=u_n-stepsize*grad(u_n);
        tn1=(1+sqrt(1+4*tn^2))/2;
        u_n=x_n+(tn-1)/tn1*(x_n-sol);
        %updates
        sol=x_n;
        tn=tn1;
        
 
        % Update
        prev_l2 = norm_l2;
        iter_L2 = iter_L2 + 1;
        
    end
end

% Log after the projection onto the L2-ball
if param.verbose >= 1
    fprintf(['  prox_L2: ||A w x- y||_2^2 = %e,', ...
        ' %s, iter = %i\n'], norm_l2, crit_L2, iter_L2);
end
iter=iter_L2;

sol = reshape(sol,siz);
end

