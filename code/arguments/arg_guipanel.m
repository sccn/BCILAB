function result = arg_guipanel(varargin)
% Create a uipanel that displays an argument property inspector for a Function.
% Result = arg_guipanel(Options ...)
% Result = arg_guipanel(Parent, Options ...)
%
% The handle supports the method .GetPropertySpecification(), by means of which the edited argument
% specification can be retrieved. This result can be turned into a valid Function argument using
% arg_tovals(). Additional Parameters may be passed to the Function, in order to override some of
% its defaults.
%
% In:
%   Parent : optional parent widget
%
%   Options : name-value pairs; possible names are:
%              'Function' : the function for which to display arguments
%
%              'Parameters' : cell array of parameters to the function
%
%              'Position' : position of the panel within the parent widget
%
%              'PanelOnly' : if true, generate only a uipanel that can be embedded in a dialog;
%                            otherwise generate a figure and wait until it is closed.
%                            (default: true)
%
% Out:
%   Result : * if PanelOnly, this is the handle to the panel; supports .GetPropertySpecification() 
%              to obain the edited specification (when done)
%            * otherwise this is the PropertySpecification of the function at the time when the 
%              figure is closed
%
% Examples:
%   % get a uipanel that allows to edit parameters to a function
%   f = figure;
%   h = arg_guipanel(f,'Function',@myfunction);
%
%   % get a uipanel that allows to edit parameters to a function, and customize initial settings
%   f = figure;
%   h = arg_guipanel(f,'Function',@myfunction,'Parameters',{3,21,'test'});
%
%   % get a uipanel that allows to edit parameters to a function, and put it in a specific position
%   f = figure;
%   h = arg_guipanel(f,'Function',@myfunction,'Position',[0 0 100 200]);
%
% See also:
%   arg_guidialog, arg_guidialog_ex, arg_define
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-10-24

% Copyright (C) Christian Kothe, SCCN, 2010, christian@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License as published by the Free Software Foundation; either version 2 of the
% License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
% even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License along with this program; if not,
% write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
% USA

global tracking;
if ~isfield(tracking,'gui') || ~isfield(tracking.gui,'show_guru')
    tracking.gui.show_guru = true; end

% separate the parent, if specified
if isscalar(varargin{1}) && ishandle(varargin{1})
    parent = varargin{1};
    varargin = varargin(2:end);
else
    mp = get(0,'MonitorPositions')';
    parent = figure('MenuBar','none','Position',[mp(3)/2-200,mp(4)/2-200,400,400]);
end

% get the options
opts = hlp_varargin2struct(varargin, {'Function','function','func'},mandatory, {'Parameters','parameters','params'},{}, {'Position','position','pos'},[0 0 1 1],{'PanelOnly','panel_only'},true);

% obtain the argument specification for the function
spec = arg_report('rich', opts.Function, opts.Parameters);
% ... and turn it into an array of PropertyGridField's
properties = PropertyGridField.GenerateFrom(spec);

% instantiate the grid
args = hlp_struct2varargin(opts,'suppress',{'Function','Parameters','PanelOnly'});
hpanel = PropertyGrid(parent,args{:},'Properties',properties,'ShowExpert',tracking.gui.show_guru);

if ~opts.PanelOnly
    % in this case a figure is generated and we wait until the figure is closed
    set(hpanel.Parent,'CloseRequestFcn',@(hfig,v)extract_properties(hfig));
    uiwait(hpanel.Parent);
else
    result = hpanel;
end

function extract_properties(hfig)
    result = arg_tovals(hpanel.GetPropertySpecification);
    delete(hfig);
end

end
