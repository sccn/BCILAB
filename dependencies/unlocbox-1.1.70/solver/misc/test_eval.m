function [F] = test_eval(F)
%TEST_NORM Test if the norm is given with the actual way
%   Usage: F = test_eval(F);
%          test_eval(F);
%
%   Input parameters:
%         F         : Scuctures or array of structures
%
%   Output parameters:
%         F         : Scuctures or array of structures
%
%   This function test if the user still use the old .norm name instead of
%   .eval. In this case it will fix the bug by coppying .norm to
%   .eval and display a warning.
%
%   This function also check if the function contain a .norm or a
%   .eval. If not, it returns an error.
%
%   At convergence, the flag stop is set to one.
%
%
%   Url: http://unlocbox.sourceforge.net/doc/solver/misc/test_eval.php

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

%   Nathanael Perraudin
%   Date: 22 jan 2013

% Optional input arguments

% param
if nargin<1, F=struct; end

% number of function
m = length(F);
if m>1

    for ii=1:m
        [F2(ii)] = test_eval(F(ii)); %#ok<AGROW>
    end
    F=F2;
else
    flag=0;
    
    if ~isfield(F, 'eval')
        if ~isfield(F, 'norm')
            error('No eval function!');
        else
            flag=1;
            F.eval=F.norm;
        end
    end

    if flag
       fprintf('WARNING!!! You are using the old name "norm" to evalute a function. \n    Please upgrade to the new name: "eval". \n    This name will not be supported anymore in the future versions of the UNLocBoX.\n'); 
    end
end

