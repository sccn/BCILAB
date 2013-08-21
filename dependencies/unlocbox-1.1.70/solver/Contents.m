% Unlocbox - Solvers
%
%  General solvers
%    ADMM               -  Alternating-direction method of multipliers algorithm.
%    DOUGLAS_RACHFORD   -  Foward backward splitting algorithm.
%    FORWARD_BACKWARD   -  Foward backward splitting algorithm.
%    GENERALIZED_FORWARD_BACKWARD - Generaliyed foward backward splitting algorithm.
%    PPXA               -  Parallel Proximal algorithm.
%    SDMM               -  Simultaneous-direction method of multipliers algorithm.
%    GRADIENT_DESCENT   -  Simple gradient descent solver.
%    BACKWARD_BACKWARD  -  Backward backward algorithm
%    POCS               -  Projection onto convex setss
%
%  Composed solver
%    RLR                -  Regularized Linear Regresssion solver (special case of forward-backward)
%    SOLVE_BPDN         -  Solve a BPDN (Basis Pursuit denoising) problem
%    SOLVE_TVDN         -  Solve a TVDN (TV denoising) problem
%    SOLVE_LRJS         -  Solve a LRJS (low rank joint sparsity) problem
%
%  For help, bug reports, suggestions etc. please send email to
%  unlocbox-help@lists.sourceforge.net
%
%   Url: http://unlocbox.sourceforge.net/doc/solver/Contents.php

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
