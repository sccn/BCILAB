function varargout = lfgui(varargin)
% LFGUI M-file for lfgui.fig
%      LFGUI, by itself, creates a new LFGUI or raises the existing
%      singleton*.
%
%      H = LFGUI returns the handle to a new LFGUI or the handle to
%      the existing singleton*.
%
%      LFGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LFGUI.M with the given input arguments.
%
%      LFGUI('Property','Value',...) creates a new LFGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before lfgui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to lfgui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help lfgui

% Last Modified by GUIDE v2.5 27-Dec-2005 18:56:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @lfgui_OpeningFcn, ...
                   'gui_OutputFcn',  @lfgui_OutputFcn, ...
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


% --- Executes just before lfgui is made visible.
function lfgui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to lfgui (see VARARGIN)

fit = locfit(varargin{:});
lfplot(fit);
handles.lfargs = varargin;

% Choose default command line output for lfgui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes lfgui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = lfgui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

n = get(hObject,'Value');
n0 = get(hObject,'Min');
n1 = get(hObject,'Max');
nn = 0.1+(n-n0)/(n1-n0);
fit = locfit(handles.lfargs{:},'nn',nn);
lfplot(fit);

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


