function [sol, iter] = prox_tv(b, gamma, param)
%PROX_TV Total variation proximal operator
%   Usage:  sol=prox_tv(x, gamma)
%           sol=prox_tv(x, gamma,param)
%           [sol, iter]=prox_tv(...)
%
%   Input parameters:
%         x     : Input signal.
%         gamma : Regularization parameter.
%         param : Structure of optional parameters.
%   Output parameters:
%         sol   : Solution.
%         iter  : Number of iterations at convergence.
%
%   PROX_TV(y, gamma, param) solves:
%
%      sol = argmin_{z} 0.5*||x - z||_2^2 + gamma  |X|TV
%
%
%   param is a Matlab structure containing the following fields:
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
%    param.useGPU : Use GPU to compute the TV prox operator. Please prior 
%     call init_gpu and free_gpu to launch and release the GPU library (default: 0).
%
%    param.verbose : 0 no log, 1 a summary at convergence, 2 print main
%     steps (default: 1)
%
%   See also:  prox_l1 prox_tv3d
%
%
%   References:
%     A. Beck and M. Teboulle. Fast gradient-based algorithms for constrained
%     total variation image denoising and deblurring problems. Image
%     Processing, IEEE Transactions on, 18(11):2419-2434, 2009.
%
%   Url: http://unlocbox.sourceforge.net/doc/prox/prox_tv.php

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


% Author: Gilles Puy Nathanael Perraudin
% Date: October 15, 2010
%

% Optional input arguments

if nargin<3, param=struct; end

if ~isfield(param, 'tol'), param.tol = 10e-4; end
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'maxit'), param.maxit = 200; end
if ~isfield(param, 'useGPU'), param.useGPU = 0; end

siz = size(b);

% handle the shape property
if isfield(param,'shape') && ~isempty(param.shape)
    b = reshape(b,param.shape(1),param.shape(2),[]); end
if ndims(b) > 2 %#ok<ISMAT>
    param = rmfield(param,'shape');
    sol = zeros(size(b));
    for d = 1:size(b,3)
        sol(:,:,d) = prox_tv(b(:,:,d),gamma,param); end
    sol = sol(:);
else
    
    % Test of gamma
    gamma=test_gamma(gamma);
    
    % Initializations
    [r, s] = gradient_op(b*0);
    pold = r; qold = s;
    told = 1; prev_obj = 0;
    
    % Main iterations
    if param.verbose > 1
        fprintf('  Proximal TV operator:\n');
    end
    
    if param.useGPU
        r_ptr = libpointer('singlePtr', single(r));
        s_ptr = libpointer('singlePtr', single(s));
        if sum(imag(b(:)))
            warning('Data are converted to real values\n');
            b = real(b);
        end
        b_ptr = libpointer('singlePtr', single(b));
        sol = zeros(size(b, 1), size(b, 2));
        sol_ptr = libpointer('singlePtr', single(sol));
        Depth = 1;
        Width  = size(sol, 1);
        Height = size(sol, 2);
        obj = 0; rel_obj = 0; crit_TV = 'test';
        
        calllib('gpuImgLib','ProxTV', sol_ptr, b_ptr, r_ptr, s_ptr, Height, Width, Depth, gamma, param.maxit, param.tol);
        sol = reshape(sol_ptr.Value, [size(b, 1) size(b, 2) 1]);
        sol = double(sol);
        iter = param.maxit;
    else
        for iter = 1:param.maxit
            
            % Current solution
            sol = b - gamma*div_op(r, s);
            
            % Objective function value
            obj = .5*norm(b(:)-sol(:), 2)^2 + gamma * tv_norm(sol);
            rel_obj = abs(obj-prev_obj)/obj;
            prev_obj = obj;
            
            % Stopping criterion
            if param.verbose>1
                fprintf('   Iter %i, obj = %e, rel_obj = %e\n', ...
                    iter, obj, rel_obj);
            end
            if rel_obj < param.tol
                crit_TV = 'TOL_EPS'; break;
            end
            
            % Udpate divergence vectors and project
            [dx, dy] = gradient_op(sol);
            r = r - 1/(8*gamma) * dx; s = s - 1/(8*gamma) * dy;
            weights = max(1, sqrt(abs(r).^2+abs(s).^2));
            p = r./weights; q = s./weights;
            
            % FISTA update
            t = (1+sqrt(4*told^2))/2;
            r = p + (told-1)/t * (p - pold); pold = p;
            s = q + (told-1)/t * (q - qold); qold = q;
            told = t;
            
        end
    end
    
    % Log after the minimization
    if ~exist('crit_TV', 'var'), crit_TV = 'MAX_IT'; end
    if param.verbose >= 1
        fprintf(['  Prox_TV: obj = %e, rel_obj = %e,' ...
            ' %s, iter = %i\n'], obj, rel_obj, crit_TV, iter);
    end
    
end

sol = reshape(sol,siz);

end
