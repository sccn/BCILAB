function varargout = gui_configpaths(varargin)
% GUI_CONFIGPATHS M-file for gui_configpaths.fig
%      GUI_CONFIGPATHS, by itself, creates a new GUI_CONFIGPATHS or raises the existing
%      singleton*.
%
%      H = GUI_CONFIGPATHS returns the handle to a new GUI_CONFIGPATHS or the handle to
%      the existing singleton*.
%
%      GUI_CONFIGPATHS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_CONFIGPATHS.M with the given input arguments.
%
%      GUI_CONFIGPATHS('Property','Value',...) creates a new GUI_CONFIGPATHS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_configpaths_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_configpaths_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_configpaths

% Last Modified by GUIDE v2.5 19-Nov-2010 01:57:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_configpaths_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_configpaths_OutputFcn, ...
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


% --- Executes just before gui_configpaths is made visible.
function gui_configpaths_OpeningFcn(hObject, eventdata, handles, varargin)
global tracking;
set(handles.edit1,'String',hlp_tostring(tracking.paths.data_paths));
set(handles.edit2,'String',tracking.paths.store_path);
set(handles.edit3,'String',tracking.paths.temp_path);

handles.output = hObject;
guidata(hObject, handles);
uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = gui_configpaths_OutputFcn(hObject, eventdata, handles) 
try
    varargout{1} = handles.output;
    close(handles.figure1);
catch
    varargout = {};
end

function edit1_Callback(hObject, eventdata, handles)
global tracking;
gui_check_dirnames(hObject,hlp_tostring(tracking.paths.data_paths));

function pushbutton1_Callback(hObject, eventdata, handles)
dir = evalin('base',get(handles.edit1,'String'));
newdir = uigetdir(dir{1},'Select data directory');
if ischar(newdir) && ~isempty(newdir)
    set(handles.edit1,'String',hlp_tostring({newdir})); end

function edit2_Callback(hObject, eventdata, handles)
global tracking;
gui_check_dirname(hObject,tracking.paths.store_path);

function pushbutton2_Callback(hObject, eventdata, handles)
newdir = uigetdir(get(handles.edit2,'String'),'Select results directory');
if ischar(newdir) && ~isempty(newdir)
    set(handles.edit2,'String',newdir); end

function edit3_Callback(hObject, eventdata, handles)
global tracking;
gui_check_dirname(hObject,tracking.paths.temp_path);

function pushbutton3_Callback(hObject, eventdata, handles)
newdir = uigetdir(get(handles.edit3,'String'),'Select temp directory');
if ischar(newdir) && ~isempty(newdir)
    set(handles.edit3,'String',newdir); end

% OK pressed
function pushbutton4_Callback(hObject, eventdata, handles)
global tracking;
% set new state
tracking.paths.data_paths = evalin('base',get(handles.edit1,'String'));
tracking.paths.store_path = get(handles.edit2,'String');
tracking.paths.temp_path = get(handles.edit3,'String');
% if we're still using the default temp path, just leave it [] in the config file
temp_path = tracking.paths.temp_path;
if strcmp(temp_path,[tracking.paths.bcilab_path '-temp'])
    temp_path = []; end
% make changes persistent
if utl_update_config('set', ...
        'data',hlp_tostring(tracking.paths.data_paths), ...
        'store',hlp_tostring(tracking.paths.store_path), ...
        'temp',hlp_tostring(temp_path))
    uiresume(handles.figure1);
end

function pushbutton6_Callback(hObject, eventdata, handles)
uiresume(handles.figure1);

function pushbutton7_Callback(hObject, eventdata, handles)
% help button
env_doc env_startup

function checkbox1_Callback(hObject, eventdata, handles)

function edit3_CreateFcn(hObject, eventdata, handles)

function edit2_CreateFcn(hObject, eventdata, handles)

function edit1_CreateFcn(hObject, eventdata, handles)
