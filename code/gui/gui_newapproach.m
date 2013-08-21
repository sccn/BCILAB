function varargout = gui_newapproach(varargin)
% Open a dialog to select a new BCI approach
% Approach = gui_newapproch()
%
% Out:
%   Approach : the newly created approach, or [] if the user cancelled the process
%
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2010-10-25

% Last Modified by GUIDE v2.5 10-Dec-2011 18:46:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_newapproach_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_newapproach_OutputFcn, ...
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


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% retrieve and attach the known approaches...
[handles.approaches,flat,handles.indexable_approaches] = gui_listapproaches;
guidata(hObject,handles);


% --- Executes just before gui_newapproach is made visible.
function gui_newapproach_OpeningFcn(hObject, eventdata, handles, varargin)
% invoke the popupmenu at least once
popupmenu1_Callback(handles.popupmenu1,{},guidata(hObject));
% UIWAIT makes gui_newapproach wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_newapproach_OutputFcn(hObject, eventdata, handles) 
try
    varargout = handles.output;
    full_edit = get(handles.checkbox1,'Value') == get(handles.checkbox1,'Max');
    % ... and close the figure
    close(handles.figure1);
catch
    varargout = {[],'Cancel'};
end

if strcmp(varargout{2},'OK')
    if full_edit
        % full review
        [varargout{1:2}] = gui_reviewapproach(varargout{1},true);
    else
        % quick edit
        [varargout{1:2}] = gui_configapproach(varargout{1},true);
    end
else
    varargout = {[],'Cancel'};
end

% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
idx = gui_select_bracketed(hObject);
% store the current approach in handles
handles.approach = handles.indexable_approaches{idx};
% update the description
set(handles.edit1,'String',handles.approach.description);
% write back
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% populate with pretty-printed list
set(hObject,'String',gui_print_grouplist(handles.approaches,@(x)x.name));


function pushbutton1_Callback(hObject, eventdata, handles)
handles.output = {handles.approach,'OK'};
guidata(hObject,handles);
uiresume(handles.figure1);

function pushbutton3_Callback(hObject, eventdata, handles)
handles.output = {[],'Cancel'};
guidata(hObject,handles);
uiresume(handles.figure1);

function pushbutton4_Callback(hObject, eventdata, handles)
env_doc code/paradigms


function edit1_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)

function checkbox1_Callback(hObject, eventdata, handles)


function pushbutton1_KeyPressFcn(hObject, eventdata, handles)


function figure1_KeyPressFcn(hObject, eventdata, handles)
if strcmp(eventdata.Key,'return')
    pushbutton1_Callback(handles.pushbutton1, eventdata, handles); end
if strcmp(eventdata.Key,'escape')
    pushbutton3_Callback(handles.pushbutton3, eventdata, handles); end


function popupmenu1_KeyPressFcn(hObject, eventdata, handles)
if strcmp(eventdata.Key,'return')
    pushbutton1_Callback(handles.pushbutton1, eventdata, handles); end
if strcmp(eventdata.Key,'escape')
    pushbutton3_Callback(handles.pushbutton3, eventdata, handles); end
