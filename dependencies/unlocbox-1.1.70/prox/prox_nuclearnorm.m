function [sol] = prox_nuclearnorm(x, gamma, param)
%PROX_NUCLEARNORM Proximal operator with the nuclear norm
%   Usage:  sol=prox_nuclearnorm(x, gamma)
%           sol=prox_nuclearnorm(x, gamma, param)
%
%   Input parameters:
%         x     : Input signal.
%         gamma : Regularization parameter.
%         param : Structure of optional parameters.
%   Output parameters:
%         sol   : Solution.
%
%   prox_NuclearNorm(x, gamma, param) solves:
%
%      sol = min_{z} 0.5*||x - z||_2^2 + gamma  |X|*
%
%   
%   param is a Matlab structure containing the following fields:
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
%    param.verbose : 0 no log, 1 a summary at convergence, 2 print main
%     steps (default: 1)
%
%   See also:  prox_l1 proj_b1 prox_tv
%
%   Url: http://unlocbox.sourceforge.net/doc/prox/prox_nuclearnorm.php

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
% Date: June 2012 EPFL
%

if nargin<3, param=struct; end

% Optional input arguments
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'At'), param.At = @(x) x; end
if ~isfield(param, 'A'), param.A = @(x) x; end
if ~isfield(param, 'nu'), param.nu = 1; end

siz = size(x);

x = param.A(x);

% handle the shape property
if isfield(param,'shape') && ~isempty(param.shape)
    x = reshape(x,param.shape(1),param.shape(2),[]); end
if ndims(x) > 2 %#ok<ISMAT>
    param = rmfield(param,'shape');
    sol = zeros(size(x));
    for d = 1:size(x,3)
        sol(:,:,d) = prox_nuclearnorm(x(:,:,d),gamma,param); end
    sol = sol(:);
else
    % Test of gamma
    gamma=test_gamma(gamma);
    
    % Useful functions
    soft = @(z, T) sign(z).*max(abs(z)-T, 0);
     
    [U,Sigma,V] =  svd(x,'econ');
        
    % Shrink:
    sigma = diag(Sigma);
    sigma = soft(sigma,gamma);
    r = sum(sigma > 0);
    U = U(:,1:r); V = V(:,1:r); Sigma = diag(sigma(1:r)); nuclearNorm = sum(diag(Sigma));
    sol = U*Sigma*V.';
        
    if param.verbose >= 1
        fprintf('  prox nuclear norm: rank= %i, |x|_* = %e \n', r, nuclearNorm);
    end
end

% finally map back
sol = 1/param.nu*param.At(sol);
sol = reshape(sol,siz);

end
