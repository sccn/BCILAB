function varargout = gui_producevisualization(varargin)
% GUI_PRODUCEVISUALIZATION MATLAB code for gui_producevisualization.fig
%      GUI_PRODUCEVISUALIZATION, by itself, creates a new GUI_PRODUCEVISUALIZATION or raises the existing
%      singleton*.
%
%      H = GUI_PRODUCEVISUALIZATION returns the handle to a new GUI_PRODUCEVISUALIZATION or the handle to
%      the existing singleton*.
%
%      GUI_PRODUCEVISUALIZATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_PRODUCEVISUALIZATION.M with the given input arguments.
%
%      GUI_PRODUCEVISUALIZATION('Property','Value',...) creates a new GUI_PRODUCEVISUALIZATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_producevisualization_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_producevisualization_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_producevisualization

% Last Modified by GUIDE v2.5 17-Nov-2010 20:40:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_producevisualization_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_producevisualization_OutputFcn, ...
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


% --- Executes just before gui_producevisualization is made visible.
function gui_producevisualization_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_producevisualization (see VARARGIN)

% populate model selector listbox with admissible models
vars = evalin('base','whos');
vars = {vars(strcmp({vars.class},'struct')).name};
model_vars = {};
model_times = [];
% for all struct vars...
for v=vars
    % select those that are results
    var = evalin('base',v{1});
    if isfield(var,'tracking') && ~isfield(var,'pipelines') && isfield(var.tracking,'prediction_function')
        model_vars{end+1} = v{1}; 
        if isfield(var,'timestamp')
            model_times(end+1) = var.timestamp;
        else
            model_times(end+1) = 0;
        end
    end
end

if isempty(model_vars)
    % disable the OK button
    set(handles.popupmenu1,'String','(no model available)');
    set(handles.pushbutton1,'Enable','off');
else    
    % add to pulldown menu...
    set(handles.popupmenu1,'String',model_vars);
end
guidata(hObject, handles);

% populate stream selector listbox with admissible streams
vars = evalin('base','whos');
vars = {vars(strcmp({vars.class},'struct')).name};
stream_vars = {};
% for all struct vars...
for v=vars
    % select those that are results
    var = evalin('base',v{1});
    if isfield(var,{'buffer','smax'})
        stream_vars{end+1} = v{1}; end
end

if isempty(stream_vars)
    % disable the OK button
    set(handles.popupmenu2,'String','(no stream available)');
    set(handles.pushbutton1,'Enable','off');
else    
    % add to pulldown menu...
    set(handles.popupmenu2,'String',stream_vars);
end
guidata(hObject, handles);


% Choose default command line output for gui_producevisualization
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_producevisualization wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_producevisualization_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
try
    varargout{1} = handles.output;
    close(handles.figure1);
catch
    varargout = {};
end

function pushbutton1_Callback(hObject, eventdata, handles)
% ok pressed
% get model
choices = get(handles.popupmenu1,'String');
selection = get(handles.popupmenu1,'Value');
model = evalin('base',strtok(choices{selection}));
% get stream
choices = get(handles.popupmenu2,'String');
selection = get(handles.popupmenu2,'Value');
stream = strtok(choices{selection});
% get format
choices = get(handles.popupmenu3,'String');
selection = get(handles.popupmenu3,'Value');
fmt = choices{selection};
% start visualizer
run_producevisualization('pred_model',model,'pred_name',get(handles.edit2,'String'),'in_stream',stream, ...
    'vis_func',get(handles.edit3,'String'), 'freq',str2num(get(handles.edit1,'String')), 'out_form',fmt, ...
    'create_fig',true,'start_delay',1);
% resume gui
uiresume(handles.figure1);

function pushbutton5_Callback(hObject, eventdata, handles)
% help pressed
disp('coming soon...');

function pushbutton4_Callback(hObject, eventdata, handles)
% cancel pressed
uiresume(handles.figure1);

function edit1_Callback(hObject, eventdata, handles)
gui_check_scalar(hObject,'10');


function edit2_Callback(hObject, eventdata, handles)
gui_check_varname(hObject,'lastpredictor');

% --- auto-generated junk ---

function popupmenu1_Callback(hObject, eventdata, handles)
function popupmenu1_CreateFcn(hObject, eventdata, handles)
function popupmenu2_Callback(hObject, eventdata, handles)
function popupmenu2_CreateFcn(hObject, eventdata, handles)
function edit1_CreateFcn(hObject, eventdata, handles)
function popupmenu3_Callback(hObject, eventdata, handles)
function popupmenu3_CreateFcn(hObject, eventdata, handles)
function edit2_CreateFcn(hObject, eventdata, handles)
function edit3_Callback(hObject, eventdata, handles)
function edit3_CreateFcn(hObject, eventdata, handles)
