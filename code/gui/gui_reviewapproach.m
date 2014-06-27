function varargout = gui_reviewapproach(varargin)
% bring up a modal configuration panel for the given approach
% [Approach,Action] = gui_configapproach(Approach)
%
% In:
%   Approach : an approach; struct with fields 'paradigm' and 'parameters' (and optionally 'description' and 'name')
%              or cell array {paradigm, parameter1, parameter2, ...}
%
%   DoSave: Whether to bring up a save approach gui, after clicking okay (default: false)
%
%
% Out:
%   Result : a (re-)configured version of the Approach, or the unmodified input Approach (though possibly reformatted) if the user pressed 'Cancel'
%   Action : either 'OK' or 'Cancel'
%
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2010-10-25

% Last Modified by GUIDE v2.5 14-Aug-2013 15:38:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_reviewapproach_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_reviewapproach_OutputFcn, ...
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


% --- Executes just before gui_reviewapproach is made visible.
function gui_reviewapproach_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for gui_reviewapproach
handles.output = 'Cancel';
% get the approach
if length(varargin) < 1 || isempty(varargin{1})
    handles.approach = gui_chooseapproach();
else
    handles.approach = varargin{1};
end
if length(varargin) < 2
    handles.dosave = false;
else
    handles.dosave = varargin{2};
end
handles.origapproach = handles.approach;

if isempty(handles.approach)
    % nothing to do here...
    guidata(hObject, handles);
    close(handles.figure1);
else
    fprintf('Generating approach GUI...');
    calibrate_func = handles.approach.paradigm;
    if ischar(calibrate_func)
        instance = eval(calibrate_func); %#ok<NASGU>
        calibrate_func = eval('@instance.calibrate');
    end
    % create a new sub-panel...
    handles.hProperties = arg_guipanel(handles.uipanel2, 'Function',calibrate_func,'Parameters',handles.approach.parameters);
    % Update handles structure
    guidata(hObject, handles);
    fprintf('done.\n');
end

% --- Outputs from this function are returned to the command line.
function varargout = gui_reviewapproach_OutputFcn(hObject, eventdata, handles) 
try
    % hack to get the guipanel to render on Macs...
    if ismac
        p=get(handles.uipanel2,'Position');
        set(handles.uipanel2,'Position',p.*[ 1 1 0.99 1]);
        set(handles.uipanel2,'Position',p.*[ 1 1 1/0.99 1]);
    end
    uiwait(handles.figure1);
    handles = guidata(hObject);
    if strcmp(handles.output,'OK')
        % get the edited specification
        spec = handles.hProperties.GetPropertySpecification();
        % for nodes where the children have been removed (unchecked subtoggle),
        % reinstate an arg_selection argument..
        spec = reinstate_empty_children(spec);
        % turn it into a parameter structure and store that in the approach
        handles.approach.parameters = {arg_tovals(spec)};
        % get rid of hidden references...
        handles.approach = utl_prune_handles(handles.approach);
        % Get default command line output from handles structure
        varargout = {handles.approach,handles.output};
    else
        varargout = {handles.origapproach,handles.output};
        handles.dosave = false;
    end
    ... and close the figure!
    close(handles.figure1);
    if handles.dosave
        gui_saveapproach(varargout{1}); end
catch
    % someone closed the figure...
    varargout = {[],'Cancel'};
end


% --- Reinstate empty children of unchecked subtoggles
function x = reinstate_empty_children(x)
empty_children = cellfun('isempty',{x.children});
fix_pos = empty_children & ~cellfun('isempty',{x.alternatives});
[x(fix_pos).children] = celldeal(cellfun(@(v)cached_argument('arg_selection',v),{x(fix_pos).value},'UniformOutput',false));
[x(~empty_children).children] = celldeal(cellfun(@reinstate_empty_children,{x(~empty_children).children},'UniformOutput',false));


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
handles.output = 'OK';
guidata(hObject,handles);
uiresume(handles.figure1);

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
handles.output = 'Cancel';
guidata(hObject,handles);
uiresume(handles.figure1);

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
calibrate_func = handles.approach.paradigm;
if ischar(calibrate_func)
    instance = eval(calibrate_func); %#ok<NASGU>
    calibrate_func = eval('@instance.calibrate');
end
env_doc(calibrate_func);


% --- Executes during object creation, after setting all properties.
function uipanel2_CreateFcn(hObject, eventdata, handles)


function figure1_KeyPressFcn(hObject, eventdata, handles)
if strcmp(eventdata.Key,'return')
    pushbutton1_Callback(handles.pushbutton1, eventdata, handles); end
if strcmp(eventdata.Key,'escape')
    pushbutton2_Callback(handles.pushbutton2, eventdata, handles); end


function show_guru1_Callback(hObject, eventdata, handles)
set(handles.hProperties,'ShowExpert',get(hObject,'Value')==get(hObject,'Max'));

function show_guru1_CreateFcn(hObject, eventdata, handles)
global tracking;
set(hObject,'Value',double(tracking.gui.show_guru));