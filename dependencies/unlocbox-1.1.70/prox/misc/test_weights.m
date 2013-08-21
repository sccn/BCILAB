function [weights]=test_weights(weights,warning)
%TEST_GAMMA test if the weights are corrects
%   Usage:  weights=test_weights(weights)
%           test_weights(weights)
%           weights=test_weights(weights,warning)
%           test_weights(weights,warning)
%
%   Input parameters:
%         weights : vector
%         warning: boolean
%   Output parameters:
%         weights : vector
%
%   This function test is the weights are stricly positivs
%   
%   For each element of the vector weights
%   If it is negativ, this function return an error. If it is zero
%   this function add eps to this weights and display a warning. This can be used
%   to set 0 weight to some objectiv function.
%
%   The warning parameter is a flag to display warning or not. Default=1;
%
%   Url: http://unlocbox.sourceforge.net/doc/prox/misc/test_weights.php

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

if sum(weights<0)
    error('gamma can not be negativ!');
elseif logical(sum(weights(:)==0)) && warning
    weights=weights+eps;
    fprintf(' WARNING!!! weights is 0. We add eps to weights to keep going...\n');
else
   weights=weights; 
end
    