function n = back_perm(k)
%BACK_PERM(k)
%   Usage: n = back_perm(k)
%
%   Input parameters:
%         k     : row vectors of permutation (Can be a matrix)
%   Output parameters:
%         n     : Permuted vector
%
%   Warning! k has to be a row vector.
%
%   Url: http://unlocbox.sourceforge.net/doc/misc/back_perm.php

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
% Date: October 2011
%
n=zeros(size(k));

for j=1:size(k,2)
for p=1:size(k,1)
    n(k(p),j)=p;
end
end
