function []=test_multigroup(x,g_d,g_t)
%TEST_MULTIGROUP test if the parameter g_d and g_t are correct
%   Usage:  test_multigroup(x,g_d,g_t)
%
%   Input parameters:
%         x       : vector
%         g_d     : numerics
%         g_t     : numerics
%   Output parameters:
%
%   This function 
%
%
%   Url: http://unlocbox.sourceforge.net/doc/prox/misc/test_multigroup.php

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


L=length(x); % lenght of the signal

if size(g_d,2)~=L
   fprintf(' WARNING!!! The number of collum of gd should be the size of the signal and it is not.\n'); 
end

if sum(g_t,2)~=size(g_d,2)
   error( ' The sum of the group size should be equal to the total number of element grouped. Check if g_d,g_t are row vectors!');
end


