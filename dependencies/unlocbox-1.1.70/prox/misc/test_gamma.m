function [gamma]=test_gamma(gamma,warning)
%TEST_GAMMA test if gamma is correct
%   Usage:  gamma=test_gamma(gamma)
%           test_gamma(gamma)
%           gamma=test_gamma(gamma,warning)
%           test_gamma(gamma,warning)
%
%   Input parameters:
%         gamma : number
%         warning: boolean
%   Output parameters:
%         gamma : number
%
%   This function test is gamma is stricly positiv
%   
%   If gamma is negativ, this function return an error. If gamma is zero
%   this function add eps to gamma and display a warning. This can be used
%   to set 0 weight to some objectiv function.
%
%   The warning parameter is a flag to display warning or not. Default=1;
%
%
%   Url: http://unlocbox.sourceforge.net/doc/prox/misc/test_gamma.php

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
% Date: Nov 2012
%

if nargin<2
   warning=1; 
end

if gamma<0
    error('gamma can not be negativ!');
elseif (gamma==0) && warning
    gamma=gamma+eps;
    fprintf(' WARNING!!! gamma is 0. We add eps to gamma to keep going...\n');
else
   gamma=gamma; 
end
    