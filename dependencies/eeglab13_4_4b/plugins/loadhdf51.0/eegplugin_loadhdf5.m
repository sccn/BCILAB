% eegplugin_loadhdf5() - EEGLAB plugin for importing *.hdf5 files in the EEGlab Matlab toolbox
%
% Usage:
%   >> eegplugin_biosig(fig, trystrs, catchstrs);
%
% Inputs:
%   fig        - [integer]  EEGLAB figure
%   trystrs    - [struct] "try" strings for menu callbacks.
%   catchstrs  - [struct] "catch" strings for menu callbacks. 
%
% Create a plugin:
%   For more information on how to create an EEGLAB plugin see the
%   help message of eegplugin_loadhdf5() or visit http://www.sccn.ucsd.edu/eeglab/contrib.html
%
% Author: Simon Lind Kappel, Aarhus University, 2015

% Copyright (C) 2015 Simon Lind Kappel, Aarhus University, s@lkappel.dk
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
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

% $Log: eegplugin_loadhdf5.m,v $

function vers = eegplugin_loadhdf5(fig, trystrs, catchstrs)

vers = 'loadhdf5 v1.0';

if nargin < 3
    error('eegplugin_loadhdf5 requires 3 arguments');
end;

% add readhdf5 folder to path
% -----------------------
if ~exist('eegplugin_loadhdf5')
    p = which('eegplugin_loadhdf5');
    p = p(1:findstr(p,'eegplugin_loadhdf5.m')-1);
    addpath(p);
end;

% find import data menu
% ---------------------
menu = findobj(fig, 'tag', 'import data');

% menu callbacks
% --------------
comhdf5 = [ trystrs.no_check '[EEGTMP LASTCOM] = pop_loadhdf5;' catchstrs.new_non_empty ];

% create menus
% -------------------------
% add new submenu
uimenu( menu, 'Label', 'From .hdf5 recorded using g.recorder', 'CallBack', comhdf5,'separator', 'on');