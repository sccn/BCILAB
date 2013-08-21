function varargout = gui_chooseapproach(varargin)
% GUI_CHOOSEAPPROACH MATLAB code for gui_chooseapproach.fig
%      GUI_CHOOSEAPPROACH, by itself, creates a new GUI_CHOOSEAPPROACH or raises the existing
%      singleton*.
%
%      H = GUI_CHOOSEAPPROACH returns the handle to a new GUI_CHOOSEAPPROACH or the handle to
%      the existing singleton*.
%
%      GUI_CHOOSEAPPROACH('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_CHOOSEAPPROACH.M with the given input arguments.
%
%      GUI_CHOOSEAPPROACH('Property','Value',...) creates a new GUI_CHOOSEAPPROACH or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_chooseapproach_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_chooseapproach_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_chooseapproach

% Last Modified by GUIDE v2.5 10-Dec-2011 18:39:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_chooseapproach_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_chooseapproach_OutputFcn, ...
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



function gui_chooseapproach_OpeningFcn(hObject, eventdata, handles, varargin)

% retrieve the known workspace approaches
approaches = {};
names = {};
times = [];
vars = evalin('base','whos');
for v=find(strcmp({vars.class},'struct'))
    var = evalin('base',vars(v).name);
    if isfield(var,{'paradigm','parameters'})
        % list it
        approaches{end+1} = var;
        if isfield(var,'timestamp')
            times(end+1) = var.timestamp;
        else
            times(end+1) = 0;
        end
        if isfield(var,'name')
            names{end+1} = [vars(v).name  ' ("' var.name '")']; 
        else
            names{end+1} = [vars(v).name  ' (based on ' char(var.paradigm) ')']; 
        end
    end
end
handles.approaches = approaches;
handles.approachnames = names;
guidata(hObject,handles);

if ~isempty(handles.approaches)
    % populate the popup menu with approaches...
    set(handles.popupmenu1,'String',handles.approachnames);

    % preferably select the newest one...
    selected_id = argmax(times);
    if isempty(selected_id)
        selected_id = 1; end
    set(handles.popupmenu1,'Value',selected_id);
    guidata(hObject, handles);
    
    if length(handles.approaches) > 1
        % UIWAIT makes gui_chooseapproach wait for user response (see UIRESUME)
        uiwait(handles.figure1);
    else
        handles.output = handles.approaches{1};
        guidata(hObject, handles);
    end
else
    % error dialog, if there is no approach
    errordlg2('You first need to create an approach before you can operate on it.');
    % return nothing, if there is no approach
    handles.output = [];
    % Update handles structure
    guidata(hObject, handles);
end

function pushbutton1_Callback(hObject, eventdata, handles)
% retrieve the selected approach, and return
handles.output = handles.approaches{get(handles.popupmenu1,'Value')};
guidata(hObject, handles);
uiresume(handles.figure1);


function varargout = gui_chooseapproach_OutputFcn(hObject, eventdata, handles) 
try
    varargout{1} = handles.output;
    close(handles.figure1);
catch
    varargout{1} = [];
end


% --- auto-generated stufff.... 

function popupmenu1_Callback(hObject, eventdata, handles)



function popupmenu1_CreateFcn(hObject, eventdata, handles)


function figure1_KeyPressFcn(hObject, eventdata, handles)
if strcmp(eventdata.Key,'return')
    pushbutton1_Callback(handles.pushbutton1, eventdata, handles); end
