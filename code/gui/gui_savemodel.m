function varargout = gui_savemodel(varargin)
% GUI_SAVEMODEL M-file for gui_savemodel.fig
%      GUI_SAVEMODEL, by itself, creates a new GUI_SAVEMODEL or raises the existing
%      singleton*.
%
%      H = GUI_SAVEMODEL returns the handle to a new GUI_SAVEMODEL or the handle to
%      the existing singleton*.
%
%      GUI_SAVEMODEL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_SAVEMODEL.M with the given input arguments.
%
%      GUI_SAVEMODEL('Property','Value',...) creates a new GUI_SAVEMODEL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_savemodel_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_savemodel_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_savemodel

% Last Modified by GUIDE v2.5 12-Apr-2012 21:06:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_savemodel_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_savemodel_OutputFcn, ...
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


% --- Executes just before gui_savemodel is made visible.
function gui_savemodel_OpeningFcn(hObject, eventdata, handles, varargin)
% populate listbox with all result data structures...
vars = evalin('base','whos');
vars = {vars(strcmp({vars.class},'struct')).name};
model_vars = {};
model_times = [];
% for all struct vars...
for v=vars
    % select those that are results
    var = evalin('base',v{1});
    if isfield(var,'tracking') && isfield(var.tracking,'prediction_function')        
        model_vars{end+1} = v{1}; 
        if isfield(var,'timestamp')
            model_times(end+1) = var.timestamp;
        else
            model_times(end+1) = 0;
        end
    end
end

handles.selection = [];
guidata(hObject, handles);

if isempty(model_vars)
    % disable the OK button
    set(handles.popupmenu1,'String','(no models available)');
    set(handles.pushbutton1,'Enable','off');
    uiwait(handles.figure1);
elseif length(model_vars) == 1
    % save it right away
    close(handles.figure1);
    if ~isempty(model_vars{1})
        do_save(model_vars{1}); end
else
    % add to pulldown menu...
    set(handles.popupmenu1,'String',model_vars);
    set(handles.popupmenu1,'Value',argmax(model_times));
    
    % Update handles structure
    guidata(hObject, handles);
    
    % UIWAIT makes gui_savemodel wait for user response (see UIRESUME)
    uiwait(handles.figure1);
end



% --- Outputs from this function are returned to the command line.
function varargout = gui_savemodel_OutputFcn(hObject, eventdata, handles) 
try
    selection = handles.selection;
    close(handles.figure1);
    if ~isempty(selection)
        do_save(selection); end
catch
end
varargout = {};


% OK button
function pushbutton1_Callback(hObject, eventdata, handles)
% get the selection...
contents = get(handles.popupmenu1,'String');
selection = get(handles.popupmenu1,'Value');
handles.selection = contents{selection};
guidata(hObject,handles);
uiresume(handles.figure1);

% Cancel button
function pushbutton2_Callback(hObject, eventdata, handles)
uiresume(handles.figure1);


% Perform the actual saving after a model has been selected
function do_save(selection)
io_mkdirs('home:/.bcilab/models/');
[FileName,PathName] = uiputfile('*.mdl','Save model',env_translatepath('home:/.bcilab/models/untitled.mdl'));
if FileName
    % assign as struct field
    tmp.(selection) = evalin('base',selection);
    % and then save it from there...
    io_save([PathName,FileName],'-struct tmp',selection);
end


function figure1_KeyPressFcn(hObject, eventdata, handles)
if strcmp(eventdata.Key,'return')
    pushbutton1_Callback(handles.pushbutton1, eventdata, handles); end
if strcmp(eventdata.Key,'escape')
    pushbutton2_Callback(handles.pushbutton2, eventdata, handles); end

function popupmenu1_KeyPressFcn(hObject, eventdata, handles)
if strcmp(eventdata.Key,'return')
    pushbutton1_Callback(handles.pushbutton1, eventdata, handles); end
if strcmp(eventdata.Key,'escape')
    pushbutton2_Callback(handles.pushbutton2, eventdata, handles); end

% --- auto-generated junk 
function popupmenu1_Callback(hObject, eventdata, handles)
function popupmenu1_CreateFcn(hObject, eventdata, handles)
