function y = tv_norm(u,shape)
%TV_NORM 2 Dimentional TV norm
%   Usage:  y = tv_norm(x)
%
%   Input parameters:
%         x     : Input data (2 dimentional matrix)
%   Output parameters:
%         sol   : Norm
%
%   Compute the 2-dimentional TV norm of x
%
%   Url: http://unlocbox.sourceforge.net/doc/misc/tv_norm.php

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


if nargin > 1
    u = reshape(u,shape(1),shape(2),[]); end
if ndims(u) > 2 %#ok<ISMAT>
    y = 0;
    for d = 1:size(u,3)
        y = y + tv_norm(u(:,:,d)); end
else
    [dx, dy] = gradient_op(u);    
    temp = sqrt(abs(dx).^2 + abs(dy).^2);
    y = sum(temp(:));
end