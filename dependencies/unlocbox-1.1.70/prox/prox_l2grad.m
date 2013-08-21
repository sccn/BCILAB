function [sol,iter] = prox_l2grad(x, gamma, param)
%PROX_L2grad Proximal operator of the 2 norm of the gradient in 1 dimension
%   Usage:  sol=prox_l2grad(x, gamma)
%           sol=prox_l2grad(x, gamma, param)
%           [sol, iter]=prox_l2grad(x, gamma, param)
%
%   Input parameters:
%         x     : Input signal.
%         gamma : Regularization parameter.
%         param : Structure of optional parameters.
%   Output parameters:
%         sol   : Solution.
%         iter  : Number of iterations at convergence.
%
%   This function compute the 1 dimensional proximal operator of x. For
%   matrices, the function is applied to each column. For N-D
%   arrays, the function operates on the first 
%   dimension.
%
%   PROX_L2GRAD(x, gamma, param) solves:
%
%      sol = argmin_{z} 0.5*||x - z||_2^2 + gamma  |GRAD(Z)|2^2
%
%
%   param is a Matlab structure containing the following fields:
%
%    param.abasis : to use another basis than the DFT (default: 0). To be
%                     done -- Not working yet
%
%    param.weights : weights if you use a an array.
%
%    param.verbose : 0 no log, 1 a summary at convergence, 2 print main
%     steps (default: 1)
%
%    param.deriveorder : Order ot the derivative default 1
%
%   See also:  proj_b1 prox_l1inf prox_l12 prox_tv
%
%   Url: http://unlocbox.sourceforge.net/doc/prox/prox_l2grad.php

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
if ~isfield(param, 'weights'), param.weights = 1; end
if ~isfield(param, 'deriveorder'), param.deriveorder = 1; end
if ~isfield(param, 'abasis'), param.abasis = 0; end

warning=0;
% test the parameters
test_gamma(gamma,warning);
test_weights(param.weights,warning);

if param.abasis
    error('This option is not currently supported, please contact the developer')
end

p=param.deriveorder;

% useful function
h=@(t) 1./(1+param.weights(:)*gamma*t').^p;

% size of the signal

L=size(x,1);
l=(0:L-1)';
lambda=2-2*cos(2*pi*l/L);

y=fft(x);

sol=ifft(y.*h(lambda)');

% one iteration
iter=1;

curr_norm=norm(gradient(x))^2;

% Summary
if param.verbose>=1
   fprintf('  Prox_l2grad: %i iteration(s), ||grad(x)||^2=%g\n',iter,curr_norm);
end



