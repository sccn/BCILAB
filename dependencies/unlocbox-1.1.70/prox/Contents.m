% Unlocbox - Proximal operators
%
%   The proximal operator of function f evaluated in z is the solution of 
%   the folowing equation:   
%
%      prox_{f, gamma }(x)=min_x  1/2 |X-Z|2^2 + gamma f(z)
%
%
%   Here are a list of common usual proximal operators available in the
%   UnLocBoX. We remember the reader that projections are particular cases of proximal
%   operators.
%
%  General Proximal operators
%    FAST_PROJ_B2       -  Fast projection on a B2-Ball
%    PROJ_B1            -  Projection on a B1-Ball
%    PROJ_B2            -  Projection on a B2-Ball
%    PROX_L1            -  Proximal operator of the L1 norm
%    PROX_L2            -  Proximal operator of the L2 norm
%    PROX_L2grad        -  Proximal operator of the L1 norm
%    PROX_L1inf         -  Proximal operator of the gradient of the L2 norm
%    PROX_L12           -  Proximal operator of the L12 norm
%    PROX_NUCLEARNORM   -  Proximal operator of the nuclear norm
%    PROX_SUMG          -  Proximal operator of a sum of function
%    PROX_TV            -  Proximal operator of the TV norm
%    PROX_TV3D          -  Proximal operator of the 3D TV norm
%
%  For help, bug reports, suggestions etc. please send email to
%  unlocbox-help@lists.sourceforge.net
%
%   Url: http://unlocbox.sourceforge.net/doc/prox/Contents.php

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
