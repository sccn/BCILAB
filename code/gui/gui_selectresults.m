function varargout = gui_selectresults(varargin)
% GUI_SELECTRESULTS M-file for gui_selectresults.fig
%      GUI_SELECTRESULTS, by itself, creates a new GUI_SELECTRESULTS or raises the existing
%      singleton*.
%
%      H = GUI_SELECTRESULTS returns the handle to a new GUI_SELECTRESULTS or the handle to
%      the existing singleton*.
%
%      GUI_SELECTRESULTS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_SELECTRESULTS.M with the given input arguments.
%
%      GUI_SELECTRESULTS('Property','Value',...) creates a new GUI_SELECTRESULTS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_selectresults_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_selectresults_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_selectresults

% Last Modified by GUIDE v2.5 10-Dec-2011 18:43:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_selectresults_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_selectresults_OutputFcn, ...
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


% --- Executes just before gui_selectresults is made visible.
function gui_selectresults_OpeningFcn(hObject, eventdata, handles, varargin)
% populate listbox with all result data structures...
vars = evalin('base','whos');
vars = {vars(strcmp({vars.class},'struct')).name};
result_vars = {};
result_times = [];
% for all struct vars...
for v=vars
    % select those that are results
    data = evalin('base',v{1});
    if isfield(data,'is_result')
        result_vars{end+1} = v{1}; 
        if isfield(data,'timestamp')
            result_times(end+1) = data.timestamp;
        else
            result_times(end+1) = 0;
        end
    end
end


if isempty(result_vars)
    % disable the OK button
    set(handles.popupmenu1,'String','(no results to choose)');
    set(handles.pushbutton1,'Enable','off');
else
    % add to pulldown menu...
    set(handles.popupmenu1,'String',result_vars);
    % by default select the newest result
    set(handles.popupmenu1,'Value',argmax(result_times));
end

% Choose default command line output for gui_selectresults
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_selectresults wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_selectresults_OutputFcn(hObject, eventdata, handles) 
try
    varargout{1} = handles.output;
    close(handles.figure1);
catch
end



% OK button
function pushbutton1_Callback(hObject, eventdata, handles)
% get the selection...
contents = get(handles.popupmenu1,'String');
selection = get(handles.popupmenu1,'Value');
gui_reviewresults(evalin('base',contents{selection}));
uiresume(handles.figure1);

% Cancel button
function pushbutton2_Callback(hObject, eventdata, handles)
uiresume(handles.figure1);

function figure1_KeyPressFcn(hObject, eventdata, handles)
if strcmp(eventdata.Key,'return')
    pushbutton1_Callback(handles.pushbutton1, eventdata, handles); end
if strcmp(eventdata.Key,'escape')
    pushbutton2_Callback(handles.pushbutton2, eventdata, handles); end


% --- auto-generated junk 
function popupmenu1_Callback(hObject, eventdata, handles)
function popupmenu1_CreateFcn(hObject, eventdata, handles)
