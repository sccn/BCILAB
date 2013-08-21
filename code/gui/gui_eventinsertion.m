function varargout = gui_eventinsertion(varargin)
% GUI_EVENTINSERTION MATLAB code for gui_eventinsertion.fig
%      GUI_EVENTINSERTION, by itself, creates a new GUI_EVENTINSERTION or raises the existing
%      singleton*.
%
%      H = GUI_EVENTINSERTION returns the handle to a new GUI_EVENTINSERTION or the handle to
%      the existing singleton*.
%
%      GUI_EVENTINSERTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_EVENTINSERTION.M with the given input arguments.
%
%      GUI_EVENTINSERTION('Property','Value',...) creates a new GUI_EVENTINSERTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_eventinsertion_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_eventinsertion_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_eventinsertion

% Last Modified by GUIDE v2.5 28-Oct-2010 04:42:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_eventinsertion_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_eventinsertion_OutputFcn, ...
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


% --- Executes just before gui_eventinsertion is made visible.
function gui_eventinsertion_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_eventinsertion (see VARARGIN)

% Choose default command line output for gui_eventinsertion
handles.output = {};

% Update handles structure
guidata(hObject, handles);

% invoke the callback to init the dialog
edit1_Callback(handles.edit1,{},guidata(hObject));

% UIWAIT makes gui_eventinsertion wait for user response (see UIRESUME)
%uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_eventinsertion_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% hack to get the guipanel to render on Macs...
if ismac
    p=get(handles.uipanel1,'Position'); 
    set(handles.uipanel1,'Position',p.*[ 1 1 0.99 1]); 
    set(handles.uipanel1,'Position',p.*[ 1 1 1/0.99 1]);
end
uiwait(handles.figure1);
handles = guidata(hObject);

try
    % Get default command line output from handles structure
    varargout{1} = handles.output;
    close(handles.figure1);
catch
    varargout{1} = {};
end


% --- set number of insertion operations!

function edit1_Callback(hObject, eventdata, handles)
gui_check_scalar(hObject,'1');
num = max(1,round(str2num(get(hObject,'String'))));
set(hObject,'String',num2str(num));

if ~isfield(handles,'tabpanel1') || length(get(handles.tabpanel1,'Title')) ~= num
    % generate titles
    title = {};
    for k=1:num
        title{end+1} = ['Step ' num2str(k)]; end
    % save old parameters
    if ~isfield(handles,'opgrids')
        handles.opgrids = {}; end    
    parameters = cell(1,num);
    for h=1:min(num,length(handles.opgrids))
        parameters{h} = arg_tovals(handles.opgrids{h}.GetPropertySpecification(),false); end
    % delete old propertygrids
    for h=1:length(handles.opgrids)
        delete(handles.opgrids{h}); end    
    % (re-) create tabpanel
    if isfield(handles,'tabpanel1')
        item = num;
        delete(handles.tabpanel1); 
    else
        item = 1;
    end
    handles.tabpanel1 = uitabpanel('active',item,'parent',handles.uipanel2,'title',title,'TitleBackgroundColor',[0.93,0.91,0.85]*0.66,'PanelBackgroundColor',[0.93,0.91,0.85],'frameBackgroundColor',[0.66,0.76,1],'TitleForegroundColor',[0 0 0]);    
    %set(handles.tabpanel1,'SelectedItem',item);
    subhandles = getappdata(handles.tabpanel1,'panels');
    % create new propertygrids
    handles.opgrids = {};
    for h=1:length(subhandles)
        handles.opgrids{h} = arg_guipanel(subhandles(h), 'Function',@set_insert_markers, 'Parameters',parameters(h)); end
end
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get the parameter specification for all steps
handles.output = {};
for h=1:length(handles.opgrids)
    handles.output{h} = arg_tovals(handles.opgrids{h}.GetPropertySpecification(),false); end
guidata(hObject,handles);
uiresume(handles.figure1);

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.output = {};
guidata(hObject,handles);
uiresume(handles.figure1);

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
env_doc set_insert_markers

% --- Executes during object creation, after setting all properties.
function uipanel1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
