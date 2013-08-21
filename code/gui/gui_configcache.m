function varargout = gui_configcache(varargin)
% GUI_CONFIGCACHE MATLAB code for gui_configcache.fig
%      GUI_CONFIGCACHE, by itself, creates a new GUI_CONFIGCACHE or raises the existing
%      singleton*.
%
%      H = GUI_CONFIGCACHE returns the handle to a new GUI_CONFIGCACHE or the handle to
%      the existing singleton*.
%
%      GUI_CONFIGCACHE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_CONFIGCACHE.M with the given input arguments.
%
%      GUI_CONFIGCACHE('Property','Value',...) creates a new GUI_CONFIGCACHE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_configcache_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_configcache_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_configcache

% Last Modified by GUIDE v2.5 09-Dec-2011 23:02:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @gui_configcache_OpeningFcn, ...
    'gui_OutputFcn',  @gui_configcache_OutputFcn, ...
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


% --- Executes just before gui_configcache is made visible.
function gui_configcache_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_configcache (see VARARGIN)
global tracking;

% get the current cache setup from the config file
try
    handles.cache = hlp_config(tracking.configscript,'get','cache');
catch
    handles.cache = {};
end
% make sure that the first entry is canonicalized
if isempty(handles.cache)
    handles.cache = {{}}; end
if ~iscell(handles.cache{1})
    handles.cache = {handles.cache}; end
handles.cache{1} = hlp_varargin2struct(handles.cache{1},'dir','','tag','location_1','time',30,'free',0.1);
set(handles.edit2,'String',handles.cache{1}.dir);
set(handles.edit3,'String',num2str(handles.cache{1}.time));
set(handles.edit4,'String',num2str(handles.cache{1}.free));

% get the mem capacity from the config file
try
    handles.memcap = hlp_config(tracking.configscript,'get','mem_capacity');
catch
    handles.memcap = 2;
end
set(handles.edit5,'String',num2str(handles.memcap));

% Choose default command line output for gui_configcache
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_configcache wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_configcache_OutputFcn(hObject, eventdata, handles)
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





function pushbutton6_Callback(hObject, eventdata, handles)
global tracking;
memcap = str2num(get(handles.edit5,'String'));
% take original cache setting
cache = handles.cache;
% override with new state
cache{1}.dir = get(handles.edit2,'String');
cache{1}.time = str2num(get(handles.edit3,'String'));
cache{1}.free = str2num(get(handles.edit4,'String'));
if isempty(cache{1}.dir)
    % clear 1st location if dir cleared
    cache{1} = {};
else
    % turn 1st location back to cell array of NVPs
    cache{1} = hlp_struct2varargin(cache{1});
end
% reduce to a simpler format if there is just 1 entry
if length(cache) == 1
    cache = cache{1}; end

if strcmp('OK',questdlg2('Changes will only take effect after you restart BCILAB.','Notice','OK','Cancel','OK'))
    if utl_update_config('set', ...
        'cache',hlp_tostring(cache), ...
        'mem_capacity',hlp_tostring(memcap))
        uiresume(handles.figure1);
    end
end

function pushbutton7_Callback(hObject, eventdata, handles)
uiresume(handles.figure1);


function pushbutton3_Callback(hObject, eventdata, handles)
set(handles.edit2,'String','');

function pushbutton8_Callback(hObject, eventdata, handles)
env_doc env_startup

function pushbutton1_Callback(hObject, eventdata, handles)
newdir = uigetdir(get(handles.edit2,'String'),'Select storage directory');
if ischar(newdir) && ~isempty(newdir)
    set(handles.edit2,'String',newdir); end

function edit3_Callback(hObject, eventdata, handles)
gui_check_scalar(hObject,handles.cache{1}.time);

function edit4_Callback(hObject, eventdata, handles)
gui_check_scalar(hObject,handles.cache{1}.free);

function edit5_Callback(hObject, eventdata, handles)
gui_check_scalar(hObject,handles.memcap);

function edit2_Callback(hObject, eventdata, handles)
gui_check_dirname(hObject,handles.cache{1}.dir);

function edit2_CreateFcn(hObject, eventdata, handles)
function edit3_CreateFcn(hObject, eventdata, handles)
function edit4_CreateFcn(hObject, eventdata, handles)
function edit5_CreateFcn(hObject, eventdata, handles)
