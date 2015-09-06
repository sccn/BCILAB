function dipfit = hlp_chanlocs2dipfit(varargin)
% create dipfit structure containing 3D channel locations
%
% Author: Tim Mullen 2012, SCCN/INC, UCSD.
% Email:  tim@sccn.ucsd.edu

% This function is part of the Source Information Flow Toolbox (SIFT)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

g = arg_define([0 Inf],varargin, ...
    arg_norep('chanlocs',mandatory,[],'Channel locations structure'), ...
    arg({'scalingFactor','ScalingFactor'},1000,[],'Scaling factor. Each [X,Y,Z] channel coordinate will be multiplied by this factor. This is useful for converting between unit types (e.g. meters -> mm)','type','denserealdouble') ... 
);

arg_toworkspace(g);

nbchan = length(chanlocs);

dipfit.hdmfile  = '';
dipfit.mrifile  = '';
dipfit.chanfile = '';
dipfit.chansel  = 1:nbchan;
dipfit.coordformat = 'Spherical';
dipfit.coord_transform = [0 0 0 0 0 0 1 1 1];

for i=1:nbchan
    dipfit.model(i).posxyz      = [chanlocs(i).X chanlocs(i).Y chanlocs(i).Z]*scalingFactor;
    dipfit.model(i).momxyz      = [0 0 0];
    dipfit.model(i).rv          = 0;
    dipfit.model(i).select      = 1;
    dipfit.model(i).diffmap     = [];
    dipfit.model(i).sourcepot   = [];
    dipfit.model(i).datapot     = [];
end

