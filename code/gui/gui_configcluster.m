function varargout = gui_configcluster(varargin)
% GUI_CONFIGCLUSTER MATLAB code for gui_configcluster.fig
%      GUI_CONFIGCLUSTER, by itself, creates a new GUI_CONFIGCLUSTER or raises the existing
%      singleton*.
%
%      H = GUI_CONFIGCLUSTER returns the handle to a new GUI_CONFIGCLUSTER or the handle to
%      the existing singleton*.
%
%      GUI_CONFIGCLUSTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_CONFIGCLUSTER.M with the given input arguments.
%
%      GUI_CONFIGCLUSTER('Property','Value',...) creates a new GUI_CONFIGCLUSTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_configcluster_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_configcluster_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_configcluster

% Last Modified by GUIDE v2.5 12-Apr-2012 22:13:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_configcluster_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_configcluster_OutputFcn, ...
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


% --- Executes just before gui_configcluster is made visible.
function gui_configcluster_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_configcluster (see VARARGIN)
global tracking;

% get the current cache setup from the config file
try
    handles.parallel = hlp_config(tracking.configscript,'get','parallel');
    handles.acquire_options = hlp_config(tracking.configscript,'get','acquire_options');
catch
    handles.parallel = {};
    handles.acquire_options = {};
end
% turn into a struct
handles.parallel = hlp_varargin2struct(handles.parallel, ...
    'engine','local', ...
    'pool',{'localhost:23547','localhost:23548','localhost:23549','localhost:23550','localhost:23551','localhost:23552','localhost:23553','localhost:23554'}, ...
    'policy','par_reschedule_policy');

engines = {'local', 'BLS', 'ParallelComputingToolbox', 'Reference'};
set(handles.popupmenu1,'Value',find(strcmp(engines,handles.parallel.engine)));
set(handles.edit1,'String',hlp_tostring(handles.parallel.pool));
tmp = hlp_tostring(handles.acquire_options);
set(handles.edit2,'String',tmp(2:end-1));


% Choose default command line output for gui_configcache
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);
% UIWAIT makes gui_configcache wait for user response (see UIRESUME)
uiwait(handles.figure1);



% --- Outputs from this function are returned to the command line.
function varargout = gui_configcluster_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    varargout{1} = handles.output;
    close(handles.figure1);
catch
    varargout = {};
end


function pushbutton1_Callback(hObject, eventdata, handles)
global tracking;
% get acquire options
handles.acquire_options = eval(['{' get(handles.edit2,'String') '}']);
% get engine
engines = {'local', 'BLS', 'ParallelComputingToolbox', 'Reference'};
handles.parallel.engine = engines{get(handles.popupmenu1,'Value')};
% get workers
handles.parallel.pool = eval(get(handles.edit1,'String'));
% prune redundant fields
if strcmp(handles.parallel.policy,'par_reschedule_policy')
    handles.parallel = rmfield(handles.parallel,'policy'); end
% turn back into cell-string array
handles.parallel = hlp_struct2varargin(handles.parallel);
% and assign
if utl_update_config('set','parallel',hlp_tostring(handles.parallel),'acquire_options',hlp_tostring(handles.acquire_options))
    tracking.acquire_options = handles.acquire_options;
    tracking.parallel = hlp_varargin2struct(handles.parallel);
    uiresume(handles.figure1);
end

function popupmenu1_Callback(hObject, eventdata, handles)
engines = {'local', 'BLS', 'ParallelComputingToolbox', 'Reference'};
handles.parallel.engine = engines{get(hObject,'Value')};
guidata(hObject,handles);

function edit1_Callback(hObject, eventdata, handles)
gui_check_cellstr(hObject,hlp_tostring(handles.parallel.pool));

function checkbox1_Callback(hObject, eventdata, handles)


function pushbutton3_Callback(hObject, eventdata, handles)
env_doc env_startup

function pushbutton2_Callback(hObject, eventdata, handles)
uiresume(handles.figure1);

function popupmenu1_CreateFcn(hObject, eventdata, handles)
function listbox1_Callback(hObject, eventdata, handles)
function listbox1_CreateFcn(hObject, eventdata, handles)
function edit1_CreateFcn(hObject, eventdata, handles)
function edit2_Callback(hObject, eventdata, handles)
function edit2_CreateFcn(hObject, eventdata, handles)


function pushbutton4_Callback(hObject, eventdata, handles)
% try to reuse the current arguments as a starting point...
acquireargs = get_cell_args(handles.edit2,[]);
% make a dialog for this
acquireargs = arg_guidialog(@par_getworkers_ssh,'Parameters',acquireargs,'Title','Cluster acquistion options','Invoke',false);
% and assign to the edit field, if not empty
if ~isempty(acquireargs)
    tmp = hlp_tostring(hlp_struct2varargin(acquireargs));
    set(handles.edit2,'String',tmp(2:end-1)); end

function args = get_cell_args(handle,name)
% get a cell argument pack
args = get(handle,'String');
if ~isempty(args)
    if args(1) ~= '{' || args(end) ~= '}'
        args = ['{' args '}']; end
    try
        args = evalin('base',args);
    catch e
        if isempty(name)
            % revert empty cell
            args = {};
        else
            % return [] to indicate formal failure
            warndlg2(sprintf('The %s cannot be evaluated; error message: %s',name,e.message),'Mal-formed load arguments.'); 
            args = [];
        end
    end
else
    args = {};
end
