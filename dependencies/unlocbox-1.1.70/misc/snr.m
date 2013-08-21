function snr = snr(map_init, map_recon)
%SNR Compute the SNR between two maps
%   Usage X = snr(MAP_INIT, MAP_NOISY) 
%   
%   computes the SNR between the maps MAP_INIT
%   and MAP_NOISY. The SNR is computed as:
%       10  log10( var(MAP_INIT) / var(MAP_INIT-MAP_NOISY) )
%   where var stands for the matlab built-in function that computes the
%   variance.
% 
%   See also: VAR
% 
%
%   Url: http://unlocbox.sourceforge.net/doc/misc/snr.php

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
 
% Author: Gilles Puy
% Date: 2009
% 

noise = map_init(:)-map_recon(:);
var_init = var(map_init(:));
var_den = var(noise(:));
snr = 10 * log10(var_init/var_den);

end