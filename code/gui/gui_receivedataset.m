function varargout = gui_receivedataset(varargin)
% GUI_RECEIVEDATASET MATLAB code for gui_receivedataset.fig
%      GUI_RECEIVEDATASET, by itself, creates a new GUI_RECEIVEDATASET or raises the existing
%      singleton*.
%
%      H = GUI_RECEIVEDATASET returns the handle to a new GUI_RECEIVEDATASET or the handle to
%      the existing singleton*.
%
%      GUI_RECEIVEDATASET('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_RECEIVEDATASET.M with the given input arguments.
%
%      GUI_RECEIVEDATASET('Property','Value',...) creates a new GUI_RECEIVEDATASET or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_receivedataset_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_receivedataset_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_receivedataset

% Last Modified by GUIDE v2.5 19-Nov-2010 21:07:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_receivedataset_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_receivedataset_OutputFcn, ...
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


% --- Executes just before gui_receivedataset is made visible.
function gui_receivedataset_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_receivedataset (see VARARGIN)


% update the popup contents (and datasets variables)
update_datasets(hObject);

% select dataset
handles = guidata(hObject);
if isempty([handles.datasets{2:2:end}])
    % no dataset specified: create a new one
    dataset = gui_loadset;
    % cancelled: cancel this dialog, too
    if isempty(dataset)
        return; end
    % (re-)create the dataset list
    update_datasets(hObject);
end

% preferably select the lastdata in the pulldown menu
handles = guidata(hObject);
dataids = strtok(cellstr(get(handles.popupmenu1,'String')));
selected_id = find(strcmp(dataids,'lastdata'),1);
if isempty(selected_id)
    selected_id = 1; end
set(handles.popupmenu1,'Value',selected_id);
guidata(hObject, handles);
% invoke the popup menu to propagate selection
popupmenu1_Callback(handles.popupmenu1,{},guidata(hObject));


% Choose default command line output for gui_receivedataset
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_receivedataset wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_receivedataset_OutputFcn(hObject, eventdata, handles) 
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


function update_datasets(hObject)
handles = guidata(hObject);
% retrieve and attach the known datasets...
[handles.datasets,printed,handles.indexable_datasets] = gui_listdatasets;
guidata(hObject,handles);
set(handles.popupmenu1,'String',printed);

function edit1_Callback(hObject, eventdata, handles)
gui_check_scalar(hObject,'25');

function edit2_Callback(hObject, eventdata, handles)
gui_check_varname(hObject,'laststream');

function pushbutton1_Callback(hObject, eventdata, handles)
% get dataset name
choices = get(handles.popupmenu1,'String');
selection = get(handles.popupmenu1,'Value');
dataset = strtok(choices{selection});
% start process
run_receivedataset(get(handles.edit2,'String'),dataset,str2num(get(handles.edit1,'String')));
% resume GUI
uiresume(handles.figure1);


function pushbutton2_Callback(hObject, eventdata, handles)
% cancel pressed
uiresume(handles.figure1);

function pushbutton3_Callback(hObject, eventdata, handles)
% help pressed
disp('coming soon...');

% --- auto-generated junk ---

function popupmenu1_Callback(hObject, eventdata, handles)
function popupmenu1_CreateFcn(hObject, eventdata, handles)
function edit1_CreateFcn(hObject, eventdata, handles)
function edit2_CreateFcn(hObject, eventdata, handles)
