function n = norm_sumg(x, G, weights)
% NORM_SUMG - Sum of norm
%
% n = norm_sumg(x, G, weights) gives back the sum of the norm x given in
% the structur array G
%
%
% The input arguments are:
%
%   - x : the vector which we want the norm
%
%   - G : The structur array of norm operator: the norm should be called "norm"
%           For instance G=[norm1 norm2 norm3]
%
%
%   - weights: to weight differently the different norm
%       (default = 1)
%
%
% Author: Gilles Puy, Nathanael Perraudin
% E-mail: gilles.puy@epfl.ch, nathanael.perraudin@epfl.ch
% Date: October 2011
%
%
%   Url: http://unlocbox.sourceforge.net/doc/misc/norm_sumg.php

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


% Optional input arguments
if nargin<2, G=struct([]); end
if nargin<3, weights=ones(size(G,2)); end


% Compute the norm
n=0;
k=1;

for i=G

    n=n+weights(k)*i.eval(x);
    k=k+1;
end

end
