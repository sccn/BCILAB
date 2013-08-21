function [sol, u, iter] = fast_proj_b2(x, ~, param)
%FAST_PROJ_B2 Projection onto a L2-ball
%   Usage:  sol=fast_proj_b2(x, ~, param)
%           [sol, u, iter]=fast_proj_b2(x, ~, param)
%
%   Input parameters:
%         x     : Input signal.
%         param : Structure of parameters.
%   Output parameters:
%         sol   : Solution.
%         u     : Residue
%         iter  : Number of iterations at convergence
%
%   FAST_PROJ_B2(x,~,param) solves:
%
%      sol = argmin_{z} |X - Z|2^2   s.t.  ||y - A z||_2 < epsilon
%
%
%   Remark: the projection is the proximal operator of the indicative function of
%   Y - A Z2 < epsilon. So it can be written:
%
%      prox_{f, gamma }(x)      where       f= i_c(||y - A z||_2 < epsilon)
%
%
%
%   param is a Matlab structure containing the following fields:
%
%    param.y : measurements (default: 0).
%
%    param.A : Forward operator (default: Id).
%
%    param.At : Adjoint operator (default: Id).
%
%    param.epsilon : Radius of the L2 ball (default = 1e-3).
%
%    param.tight : 1 if A is a tight frame or 0 if not (default = 1)
%
%    param.nu : bound on the norm of the operator A (default: 1), i.e.
%
%      ` ||A x||^2 <= nu  ||x||^2 
%
%
%    param.tol : tolerance for the projection onto the L2 ball  (default: 1e-3) . The algorithms
%     stops if
%   
%      epsilon/(1-tol) <= |Y - A Z|2 <= epsilon/(1+tol)
%
%    param.maxit : max. nb. of iterations (default: 200).
%
%    param.verbose : 0 no log, 1 a summary at convergence, 2 print main
%     steps (default: 1)
%
%   Rem: The input "~" is useless but needed for compatibility issue.
%
%   See also:  proj_b2 proj_b1
%
%   References:
%     M. Fadili and J. Starck. Monotone operator splitting for optimization
%     problems in sparse recovery. In Image Processing (ICIP), 2009 16th IEEE
%     International Conference on, pages 1461-1464. IEEE, 2009.
%     
%     A. Beck and M. Teboulle. A fast iterative shrinkage-thresholding
%     algorithm for linear inverse problems. SIAM Journal on Imaging
%     Sciences, 2(1):183-202, 2009.
%
%   Url: http://unlocbox.sourceforge.net/doc/prox/fast_proj_b2.php

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



% Author: Gilles Puy, Perraudin Nathanael
% Date: Dec. 1, 2010
%

% Optional input arguments
if ~isfield(param, 'y'), param.y = 0; end
if ~isfield(param, 'A'), param.A = @(x) x; param.At = @(x) x; end
if ~isfield(param, 'epsilon'), param.epsilon = 1e-3; end
if ~isfield(param, 'tight'), param.tight = 1; end
if ~isfield(param, 'tol'), param.tol = 1e-3; end
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'nu'), param.nu = 1; end
if ~isfield(param, 'maxit'), param.maxit = 200; end
if ~isfield(param, 'u'), param.u = zeros(size(param.y)); end

% Useful functions for the projection
sc = @(z) z*min(param.epsilon/norm(z(:)), 1); % scaling

% Projection
if (param.tight) % TIGHT FRAME CASE   
    
    temp = param.A(x) - param.y;
    sol = x + 1/param.nu * param.At(sc(temp)-temp);
    crit_B2 = 'TOL_EPS'; iter = 0;
    u = 0;
    
else % NON TIGHT FRAME CASE
    
    % Initializations
    sol = x; u = param.u; v = u;
    iter = 1; true = 1; told = 1;
    
    % Tolerance onto the L2 ball
    epsilon_low = param.epsilon/(1+param.tol);
    epsilon_up = param.epsilon/(1-param.tol);
    
    % Check if we are in the L2 ball
    dummy = param.A(sol);
    norm_res = norm(param.y(:)-dummy(:), 2);
    if norm_res <= epsilon_up
        crit_B2 = 'TOL_EPS'; true = 0;
    end
    
    % Projection onto the L2-ball
    % Init
    if param.verbose > 1
        fprintf('  Proj. B2:\n');
    end
    while true
        
        % Residual
        res = param.A(sol) - param.y; norm_res = norm(res(:), 2);
        
        % Scaling for the projection
        res = u*param.nu + res; norm_proj = norm(res(:), 2);
        
        % Log
        if param.verbose>1
            fprintf('   Iter %i, epsilon = %e, ||y - Ax||_2 = %e\n', ...
                iter, param.epsilon, norm_res);
        end
        
        % Stopping criterion
        if (norm_res>=epsilon_low && norm_res<=epsilon_up)
            crit_B2 = 'TOL_EPS'; break;
        elseif iter >= param.maxit
            crit_B2 = 'MAX_IT'; break;
        end
        
        % Projection onto the L2 ball
        t = (1+sqrt(1+4*told^2))/2;
        ratio = min(1, param.epsilon/norm_proj);
        u = v;
        v = 1/param.nu * (res - res*ratio);
        u = v + (told-1)/t * (v - u);
        
        % Current estimate
        sol = x - param.At(u);
        
        % Update number of iteration
        told = t;
        
        % Update number of iterations        
        iter = iter + 1;
        
    end
end

% Log after the projection onto the L2-ball
if param.verbose >= 1
    temp = param.A(sol);
    fprintf(['  Proj. B2: epsilon = %e, ||y-Ax||_2 = %e,', ...
        ' %s, iter = %i\n'], param.epsilon, norm(param.y(:)-temp(:)), ...
        crit_B2, iter);
end


end
