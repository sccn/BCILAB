% eegplugin_amica() - AMICA plugin 
%
% Usage:
%   >> eegplugin_amica(fig, trystrs, catchstrs);
%
% Inputs:
%   fig        - [integer] eeglab figure.
%   trystrs    - [struct] "try" strings for menu callbacks.
%   catchstrs  - [struct] "catch" strings for menu callbacks.
%
% Notes:
%   To create a new plugin, simply create a file beginning with "eegplugin_"
%   and place it in your eeglab folder. It will then be automatically 
%   detected by eeglab. See also this source code internal comments.
%   For eeglab to return errors and add the function's results to 
%   the eeglab history, menu callback must be nested into "try" and 
%   a "catch" strings. For more information on how to create eeglab 
%   plugins, see http://www.sccn.ucsd.edu/eeglab/contrib.html
%
% Author: Jason Palmer, SCCN, INC, UCSD
%
% See also: eeglab()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2008, Jason Palmer
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
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1.07  USA

% $Log: eegplugin_amica.m,v $

function vers = eegplugin_amica(fig, trystrs, catchstrs)
    
    vers = 'amica1.0';
    if nargin < 3
        error('eegplugin_amica requires 3 arguments');
    end;
    
    % add amica folder to path
    % -----------------------
    if ~exist('runamica')
        p = which('eegplugin_amica');
        p = p(1:findstr(p,'eegplugin_amica.m')-1);
        addpath(p);
    end;

    % find tools menu
    % ---------------
    menu = findobj(fig, 'tag', 'tools'); 
    % tag can be 
    % 'import data'  -> File > import data menu
    % 'import epoch' -> File > import epoch menu
    % 'import event' -> File > import event menu
    % 'export'       -> File > export
    % 'tools'        -> tools menu
    % 'plot'         -> plot menu
    
    % menu callback commands
    % ----------------------
    comrun    = 'pop_runamica(EEG);'; 
    comload   = [  trystrs.no_check '[EEG LASTCOM] = pop_loadmodout(EEG);' catchstrs.store_and_hist ];
    comchangeweights = [trystrs.no_check '[EEG LASTCOM] = pop_changeweights(EEG);' catchstrs.store_and_hist];
    comselectmodel = [ trystrs.check_data '[EEG LASTCOM] = pop_selectmodel(EEG);' catchstrs.new_and_hist ];
    complotmodelprob = [trystrs.no_check '[EEG LASTCOM] = pop_modprobplot(EEG);' catchstrs.store_and_hist];
    comtopohistplot = [ trystrs.check_data  'LASTCOM = pop_topohistplot(EEG, 0);'  catchstrs.add_to_hist];
    commoderp = [trystrs.check_data 'LASTCOM = pop_moderp(EEG,1);' catchstrs.add_to_hist];
    combackgeegplot = [trystrs.check_data  'LASTCOM = pop_ozeegplot(EEG,1,1,1);' catchstrs.add_to_hist];
    combackgcompplot = [trystrs.check_data  'LASTCOM = pop_ozeegplot(EEG,0,1,1);' catchstrs.add_to_hist];
    % create menus
    % ------------
    submenu = uimenu( menu, 'Label', 'AMICA', 'separator', 'on');
    uimenu( submenu, 'Label', 'Run AMICA'  , 'CallBack', comrun);
    uimenu( submenu, 'Label', 'Load AMICA components (beta)'  , 'CallBack', comload);
    uimenu( submenu, 'Label', 'Change current ICA weights (beta)'  , 'CallBack', comchangeweights);
    uimenu( submenu, 'Label', 'Plot Model Probabilities (beta)'  , 'CallBack', complotmodelprob);
    uimenu( submenu, 'Label', 'Plot 2-D AMICA comp. maps (beta)'  , 'CallBack', comtopohistplot);
    uimenu( submenu, 'Label', 'Scroll data + active AMICA model (beta)'  , 'CallBack', combackgeegplot);
    uimenu( submenu, 'Label', 'Scroll activations + active AMICA model (beta)'  , 'CallBack', combackgcompplot);
    uimenu( submenu, 'Label', 'Plot event-related model act.'  , 'CallBack', commoderp);
    uimenu( submenu, 'Label', 'Select data using model probability (beta)'  , 'CallBack', comselectmodel);
    
    