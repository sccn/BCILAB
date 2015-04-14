% eegplugin_cogniscan() - Cogniscan plugin
%
% Usage:
%   >> eegplugin_cogniscan(fig, trystrs, catchstrs);
%
% Inputs:
%   fig        - [integer] eeglab figure.
%   trystrs    - [struct] "try" strings for menu callbacks.
%   catchstrs  - [struct] "catch" strings for menu callbacks.
%
% Authors: Adam Gerson (adg71@columbia.edu, 2004),
%          with Lucas Parra (parra@ccny.cuny.edu, 2004)

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2004 Adam Gerson and Lucas Parra
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

function vers = eegplugin_cogniscan( fig, try_strings, catch_strings); 

vers='cogniscan1.1';

if nargin < 3
    error('eegplugin_cogniscan requires 3 arguments');
end;

% add cogniscan folder to path
% -----------------------
if ~exist('pop_readcog')
    p = which('eegplugin_cogniscan');
    p = p(1:findstr(p,'eegplugin_cogniscan.m')-1);
    addpath([ p vers ] );
end;


% create menu
menu = findobj(fig, 'tag', 'import data');

% menu callback commands
% ----------------------

cmd = [ try_strings.no_check '[EEG LASTCOM] = pop_readcog;'     catch_strings.new_and_hist ];

% add new submenu
uimenu( menu, 'label', 'From CogniScan .BIN file', 'callback', cmd, 'separator', 'on');
