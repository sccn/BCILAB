function [sol] = proj_b1(x, ~, param)
%PROJ_B1 Projection onto a L1-ball
%   Usage:  sol=proj_b1(x, ~, param)
%
%   Input parameters:
%         x     : Input signal.
%         param : Structure of parameters.
%   Output parameters:
%         sol   : Solution.
%         iter  : Number of iterations at convergence
%
%   PROJ_B1(x,~,param) solves:
%
%      sol = argmin_{z} |X - Z|2^2   s.t.  ||w.*z||_1 < epsilon
%
%
%   Remark: the projection is the proximal operator of the indicative function of
%   |W.*Z|1 < epsilon. So it can be written:
%
%      prox_{f, gamma }(x)      where       f= i_c(||w.*z||_1 < epsilon)
%
%
%   param is a Matlab structure containing the following fields:
%
%    param.epsilon : Radius of the L1 ball (default = 1e-3).
%
%    param.w : contain the weights (default ones).
%
%    param.verbose : 0 no log, 1 a summary at convergence, 2 print main
%     steps (default: 1)
%
%   Rem: The input "~" is useless but needed for compatibility issue.
%
%   This function just recalls the function oneProjector from spgl toolbox.
%   It's has been done for compatibility issues.
%
%   See also:  fast_proj_b2 proj_b1 prox_l1
%
%   Url: http://unlocbox.sourceforge.net/doc/prox/proj_b1.php

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

%
% Author: Gilles Puy, Nathanael Perraudin
% Date: October 2011
%

% Optional input arguments
if nargin<3, param=struct; end
if ~isfield(param, 'epsilon'), param.epsilon = 1e-3; end
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'w'), param.w = ones(size(x)); end

% Call the function one Projector
sol= oneProjector(x,param.w,param.epsilon);

% Log after the projection onto the L2-ball
if param.verbose >= 1
    fprintf('  Proj. B1: epsilon = %e, ||x||_2 = %e,\n', param.epsilon, norm(sol,1));
end

end