function varargout = ssg_featureselect(varargin)
% SSG_FEATURESELECT M-file for ssg_featureselect.fig
%      SSG_FEATURESELECT(axes_handle) creates a new SSG_FEATURESELECT.  Not
%      a user GUI; intended for internal use by SS GUI applications.

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ssg_featureselect_OpeningFcn, ...
                   'gui_OutputFcn',  @ssg_featureselect_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before ssg_featureselect is made visible.
function ssg_featureselect_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;


% Argument error checking
if (length(varargin) < 2)
    error('SSG:invalid_init_arg', 'Incorrect initalization syntax.');
end
handles.axes_handle = varargin{1};
handles.switchaxis = varargin{2};

if (~ishandle(handles.axes_handle) || ~strcmp(get(handles.axes_handle, 'Type'), 'axes') || ...
    isempty(guidata(handles.axes_handle)))
    error('SSG:invalid_init_arg', 'First initialization argument must be an axis handle created by an SSG function.');
end
if (~ischar(handles.switchaxis) || ~any(strcmpi({'x', 'y', 'z'}, handles.switchaxis)))
    error('SSG:invalid_init_arg', 'Second initialization argument must be one of ''x'', ''y'', or ''z''.');
end
ssg = guidata(handles.axes_handle);

choices = get(handles.popup_datatype, 'String');
value = find(strcmpi(choices, ssg.([handles.switchaxis 'choice'])));
if (isempty(value))
    error('SSG:invalid_default_choice', ['Invalid entry in initial ' handles.switchaxis 'choice.']);
end

% Fill in controls to reflect exisiting data choice.
set(handles.popup_datatype, 'Value', value);
set(handles.edit_param1, 'String', ssg.([handles.switchaxis 'param1']));
handles.oldparam = ssg.([handles.switchaxis 'param1']);
handles.maxparam = size(ssg.ss_object.waveforms, 2);
set(handles.figure_ssg_featureselect, 'Name', [upper(handles.switchaxis) '-Axis Selection']);

% Move this window to the current location of the corresponding label for convenience.
set(handles.figure_ssg_featureselect, 'Units', 'Pixel');
newpos = labelposition_screen_coords(handles);
currpos = get(handles.figure_ssg_featureselect, 'Position');
set(handles.figure_ssg_featureselect, 'Position', [(newpos-currpos(3:4)*0.75) currpos(3:4)]);
set(handles.figure_ssg_featureselect, 'Color', get(handles.axes_handle, [handles.switchaxis 'Color']));

% Update guidata structures
ssg.([handles.switchaxis 'control']) = handles.figure_ssg_featureselect;
guidata(handles.axes_handle, ssg);
guidata(hObject, handles);


% Outputs from this function are returned to the command line.
function varargout = ssg_featureselect_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;


%%%%%%%%%%%%%%%%%%%%%%%  GUI INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%
% Executes during object creation, after setting all properties.
%  Most of this is just for visual consistency.

function edit_param1_CreateFcn(hObject, eventdata, handles)
if ispc,  set(hObject,'BackgroundColor','white');
else      set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function popup_datatype_CreateFcn(hObject, eventdata, handles)
if ispc,  set(hObject,'BackgroundColor','white');
else      set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function slider_aspect_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor', [.9 .9 .9]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%% CALLBACKS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function popup_datatype_Callback(hObject, eventdata, handles)
% Make GUI state consistent before calling for data update.
choices = get(hObject, 'String');
choice = choices{get(hObject, 'Value')};
paramhandles = [handles.edit_param1, handles.button_next, handles.button_last];
switch (choice),
    case {'Signal', 'PC'}
        set(paramhandles, 'Enable', 'on');
    case {'Event Time', 'Height', 'Width', 'ISI Preceding', 'Total Energy', 'Cluster #'}
        set(paramhandles, 'Enable', 'off');
end
alter_data(handles);


function edit_param1_Callback(hObject, eventdata, handles)
% Bounds checking on the parameter edit box.
param1 = str2double(get(handles.edit_param1, 'String'));
if (isnan(param1) || (param1 ~= round(param1)))
    set(handles.edit_param1, 'String', handles.oldparam);
    beep;
else
    handles = param1_bounds_check(handles, param1);
    alter_data(handles);
end


% These next two functions just advance or review the param1 value
% in integer steps; just finger-candy to make it easier to step
% through a range of values.
function button_next_Callback(hObject, eventdata, handles)
param1 = str2double(get(handles.edit_param1, 'String'));
handles = param1_bounds_check(handles, param1 + 1);
alter_data(handles);

function button_last_Callback(hObject, eventdata, handles)
param1 = str2double(get(handles.edit_param1, 'String'));
handles = param1_bounds_check(handles, param1 - 1);
alter_data(handles);


function figure_delete_Callback(hObject, eventdata, handles)
% We need to let the parent axis know when this GUI is deleted.
ssg = guidata(handles.axes_handle);
ssg.([handles.switchaxis 'control']) = [];
guidata(handles.axes_handle, ssg);

%%%%%%%%%%%%%%%%%%%%%%%%% HELPER FUNCTIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pos = labelposition_screen_coords(handles)
% Returns [left, lower] pixel coordinates where the GUI can be placed
% (in screen coordinates) so that it appears somewhere in the vicinity
% of the axes (i tried to have them pop up on top of the axis labels, but
% there's something wrong with the units of the label's 'Extent' property
% in MatlabR13).
axeshndl = handles.axes_handle;
figurehndl = get(axeshndl, 'Parent');

oldunits = get([figurehndl, axeshndl], 'Units');
set([figurehndl, axeshndl], 'Units', 'Pixel');

axes_relative = get(axeshndl, 'Position');
figure_pos = get(figurehndl, 'Position');

set(figurehndl, 'Units', oldunits{1});
set(axeshndl, 'Units', oldunits{2});

% Offsets to get near default axis orientation.  Won't do the right thing if the
% axes have been rotated . . . its an imperfect world.
offset = axes_relative(3:4);
switch(handles.switchaxis),  
    case 'x',
        offset = offset .* [0.75 0];
    case 'y',
        offset = offset .* [-0.10 0.10];
    case 'z',
        offset = offset .* [0 0.75];
end
pos = figure_pos + axes_relative;
pos = pos(1:2) + offset;


function handles = param1_bounds_check(handles, param1)
% Checks bounds on 'param1', clipping to [1,maxparam] if necessary, then
% sets the param1 edit box to reflect the clipped value.
% Note that this function also sets the 'handles.oldparam' to the newly
% validated string but does NOT save the modified 'handles'.  Make sure to
% accept 'handles' as an output and to save if you want to remember the
% old value (useful for recovering from user typos).
if (param1 < 1)
    param1 = '1';
elseif (param1 > handles.maxparam)
    param1 = num2str(handles.maxparam);
else
    param1 = num2str(param1);
end
set(handles.edit_param1, 'String', param1);
handles.oldparam = param1;

    
function alter_data(handles)
%%% NEEDS COMMENTS!
choices = get(handles.popup_datatype, 'String');
choice = choices{get(handles.popup_datatype, 'Value')};
param1 = str2double(get(handles.edit_param1, 'String'));
ssg = guidata(handles.axes_handle);
axname = handles.switchaxis;  % shorthand

% First get the relevant data values ...
switch (choice), 
    case 'Signal', 
        vector = ssg.ss_object.waveforms(:,param1);
      
    case 'PC', % (for principal components, we might need to compute the first time)
        if (~isfield(ssg.ss_object, 'pca'))
            parentfig = get(handles.axes_handle, 'Parent');  oldptr = get(parentfig, 'Pointer');
            set([handles.figure_ssg_featureselect, parentfig], 'Pointer', 'watch');

			[pca.scores,pca.u,pca.s,pca.v] = pcasvd(ssg.ss_object.waveforms);
			ssg.ss_object.pca = pca;
			
            set(handles.figure_ssg_featureselect, 'Pointer', 'arrow');
            set(parentfig, 'Pointer', oldptr);
        end
        vector = ssg.ss_object.pca.scores(:,param1);
        
	case 'Cluster #',
        vector = ssg.ss_assigns;
		
    case 'Event Time',
        vector = ssg.ss_object.spiketimes;
        
    case 'ISI Preceding',
        vector = [NaN; diff(ssg.ss_object.spiketimes)];
        
    case 'Total Energy',
        vector = sum(ssg.ss_object.waveforms.^2, 2);
        
    case {'Height', 'Width'},  % again, we might need to compute these the first time they're used
        if (~isfield(ssg.ss_object, 'heightwidth'))
            [hw.height, hw.width] = thresholded_peaks(ssg.ss_object);
            ssg.ss_object.heightwidth = hw;
        end
        vector = ssg.ss_object.heightwidth.(lower(choice));
        
end 
ssg.([axname 'choice']) = choice;
ssg.([axname 'param1']) = num2str(param1);

% ... then feed them to the appropriate group of data.
for gp = 1:length(ssg.group_handles)
    set(ssg.group_handles(gp), [axname 'Data'], vector(ssg.group_indices{gp}));
end

% Reset the camera in case the axis scales have changed dramatically.
haxes = handles.axes_handle;
axis(haxes, 'tight');
if ~isempty(get(ssg.group_handles(1), 'ZData'))
    campos(haxes,'auto'); camtarget(haxes,'auto');  view(haxes,3);
    camva(haxes,'auto');  camup(haxes,'auto');      camva(haxes,camva(haxes)); 
end
daspect(haxes,'auto');   daspect(haxes,daspect(haxes));

% Change the axis labels to reflect the new reality.
label = ssg.([axname 'choice']);
if (strcmp(get(handles.edit_param1, 'Enable'), 'on'))
    label = [label ssg.([axname 'param1'])];
end
set(get(haxes, [axname 'Label']), 'String', label);

% Store guidata's and raise this axis control to the top.
guidata(handles.axes_handle, ssg);
guidata(handles.figure_ssg_featureselect, handles);
figure(handles.figure_ssg_featureselect);
