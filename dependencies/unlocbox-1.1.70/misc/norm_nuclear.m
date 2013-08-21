function n = norm_nuclear(x,shape)
%NORM_NUCLEAR - Nuclear norm of x
%   Usage: norm_nuclear(x,shape) 
%   
%   return the nuclear norm of x
%
%
%   The input arguments are:
%
%   - x : the matrix which we want the norm
%
%
%   Url: http://unlocbox.sourceforge.net/doc/misc/norm_nuclear.php

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
% Date: June 2012
%

if nargin > 1
    x = reshape(x,shape(1),shape(2),[]); end
if ndims(x) > 2 %#ok<ISMAT>
    n = 0;
    for d = 1:size(x,3)
        n = n + norm_nuclear(x(:,:,d)); end
else
    [~,Sigma,~] =  svd(x,'econ');
    n = sum(diag(Sigma));
end