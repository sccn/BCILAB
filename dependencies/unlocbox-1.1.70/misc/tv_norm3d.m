function y = tv_norm3d(u,shape)
%TV_NORM3D 3 Dimentional TV norm
%   Usage:  y = tv_norm(x)
%
%   Input parameters:
%         x     : Input data (3 dimentional matrix)
%   Output parameters:
%         sol   : Norm
%
%   Compute the 3-dimentional TV norm of x
%
%   Url: http://unlocbox.sourceforge.net/doc/misc/tv_norm3d.php

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
    u = reshape(u,shape(1),shape(2),shape(3),[]); end
if ndims(u) > 3
    y = 0;
    for d = 1:size(u,4)
        y = y + tv_norm3d(u(:,:,:,d)); end
else
    [dx, dy, dz] = gradient_op3d(u);    
    temp = sqrt(abs(dx).^2 + abs(dy).^2 + abs(dz).^2);
    y = sum(temp(:));
end
