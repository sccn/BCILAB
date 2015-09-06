function varargout = gui_onl_controlpanel(varargin)
% GUI_ONL_CONTROLPANEL MATLAB code for gui_onl_controlpanel.fig
%      GUI_ONL_CONTROLPANEL, by itself, creates a new GUI_ONL_CONTROLPANEL or raises the existing
%      singleton*.
%
%      H = GUI_ONL_CONTROLPANEL returns the handle to a new GUI_ONL_CONTROLPANEL or the handle to
%      the existing singleton*.
%
%      GUI_ONL_CONTROLPANEL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_ONL_CONTROLPANEL.M with the given input arguments.
%
%      GUI_ONL_CONTROLPANEL('Property','Value',...) creates a new GUI_ONL_CONTROLPANEL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_onl_controlpanel_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_onl_controlpanel_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_onl_controlpanel

% Last Modified by GUIDE v2.5 12-Jan-2014 00:29:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_onl_controlpanel_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_onl_controlpanel_OutputFcn, ...
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


% --- Executes just before gui_onl_controlpanel is made visible.
function gui_onl_controlpanel_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_onl_controlpanel (see VARARGIN)

% Choose default command line output for gui_onl_controlpanel
handles.output = hObject;

% intialize state
handles.state.startPipeline = false;
handles.state.pausePipeline = false;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_onl_controlpanel wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_onl_controlpanel_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function mnu_file_Callback(hObject, eventdata, handles)
% hObject    handle to mnu_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnu_bcilab_Callback(hObject, eventdata, handles)
% hObject    handle to mnu_bcilab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnu_sift_Callback(hObject, eventdata, handles)
% hObject    handle to mnu_sift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnu_output_Callback(hObject, eventdata, handles)
% hObject    handle to mnu_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnu_outCfg_Callback(hObject, eventdata, handles)
% hObject    handle to mnu_outCfg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnu_siftPipCfg_Callback(hObject, eventdata, handles)
% hObject    handle to mnu_siftPipCfg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnu_fltPipCfg_Callback(hObject, eventdata, handles)
% hObject    handle to mnu_fltPipCfg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnu_load_Callback(hObject, eventdata, handles)
% hObject    handle to mnu_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% load pipelines
% load the configs
[fname fpath] = uigetfile('*.mat','Load Config File');
if ~fname
    return;
end
tmp = load(fullfile(fpath,fname));
if isfield(tmp,'opts')
    handles.opts = tmp.opts;
else
    handles.opts = tmp;
end

% --------------------------------------------------------------------
function mnu_save_Callback(hObject, eventdata, handles)
% hObject    handle to mnu_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function ui_startPipeline_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to ui_startPipeline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~handles.state.startPipeline
    handles.state.startPipeline = true;
    set(handles.lblStatus,'String','Pipeline running');
end

% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function ui_startPipeline_OnCallback(hObject, eventdata, handles)
% hObject    handle to ui_startPipeline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.state.startPipeline = true;
set(handles.lblStatus,'String','Pipeline running');

% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function ui_startPipeline_OffCallback(hObject, eventdata, handles)
% hObject    handle to ui_startPipeline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.state.startPipeline = false;
set(handles.lblStatus,'String','Pipeline stopped');

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function ui_pausePipeline_OnCallback(hObject, eventdata, handles)
% hObject    handle to ui_pausePipeline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~handles.state.startPipeline
    set(hObject,'State','off');
    return;
end

handles.state.pausePipeline = true;
set(handles.lblStatus,'String','Pipeline paused');

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function ui_pausePipeline_OffCallback(hObject, eventdata, handles)
% hObject    handle to ui_pausePipeline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~handles.state.startPipeline
    set(hObject,'State','off');
    return;
end

handles.state.pausePipeline = false;
set(handles.lblStatus,'String','Pipeline running');

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function ui_streamViewer_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to ui_streamViewer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    gui_vis_filtered;
catch err
    errordlg2(err.message,'StreamViewer Error');
end


% --------------------------------------------------------------------
function mnu_Preferences_Callback(hObject, eventdata, handles)
% hObject    handle to mnu_Preferences (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function ui_senddata_OnCallback(hObject, eventdata, handles)
% hObject    handle to ui_senddata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function ui_senddata_OffCallback(hObject, eventdata, handles)
% hObject    handle to ui_senddata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
