% eegplugin_xdfimport() - EEGLAB plugin for importing XDF data files.
%                         With this menu it is possible to import a raw XDF file (*.xdf) or
%                         a compressed file (*.xdfz).
%
% Usage:
%   >> eegplugin_xdfimport(menu);
%   >> eegplugin_xdfimport(menu, trystrs, catchstrs);
%
% Inputs:
%   menu       - [float]  EEGLAB menu handle
%   trystrs    - [struct] "try" strings for menu callbacks. See notes on EEGLab plugins.
%                (http://www.sccn.ucsd.edu/eeglab/contrib.html)
%   catchstrs  - [struct] "catch" strings for menu callbacks. See notes on EEGLab plugins.
%                (http://www.sccn.ucsd.edu/eeglab/contrib.html)
%
%
% Notes:
%   This plugins consist of the following Matlab files:
%   pop_loadxdf.m           eeg_load_xdf.m
%   load_xdf.m
%
% Authors:
% Christian Kothe, Swartz Center for Computational Neuroscience UCSD, 7 May 2012
%
% See also: eeglab(), pop_loadxdf(), eeg_load_xdf(), load_xdf()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2012 Christian Kothe, ckothe@ucsd.edu
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


function vers = eegplugin_xdfimport(fig, trystrs, catchstrs)

    vers = 'xdfimport1.13';

    if nargin < 3
        error('eegplugin_xdfimport requires 3 arguments');
    end;

    % add folder to path
    % ------------------
    if ~exist('pop_loadxdf','file')
        p = which('eegplugin_xdfimport.m');
        p = p(1:findstr(p,'eegplugin_xdfimport.m')-1);
        addpath( p );
    end;

    % find import data menu
    % ---------------------
    menu = findobj(fig, 'tag', 'import data');

    % menu callbacks
    % --------------
    comcnt = [ trystrs.no_check '[EEG LASTCOM] = pop_loadxdf;' catchstrs.new_non_empty ];
    uimenu( menu, 'label', 'From .XDF or .XDFZ file', 'callback', comcnt, 'separator', 'on');
