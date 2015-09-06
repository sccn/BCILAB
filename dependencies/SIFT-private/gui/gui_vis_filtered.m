function varargout = gui_vis_filtered(varargin)
% GUI_VIS_FILTERED MATLAB code for gui_vis_filtered.fig
%      GUI_VIS_FILTERED, by itself, creates a new GUI_VIS_FILTERED or raises the existing
%      singleton*.
%
%      H = GUI_VIS_FILTERED returns the handle to a new GUI_VIS_FILTERED or the handle to
%      the existing singleton*.
%
%      GUI_VIS_FILTERED('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_VIS_FILTERED.M with the given input arguments.
%
%      GUI_VIS_FILTERED('Property','Value',...) creates a new GUI_VIS_FILTERED or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_vis_filtered_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_vis_filtered_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_vis_filtered

% Last Modified by GUIDE v2.5 11-Dec-2013 18:55:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_vis_filtered_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_vis_filtered_OutputFcn, ...
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


% --- Executes just before gui_vis_filtered is made visible.
function gui_vis_filtered_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_vis_filtered (see VARARGIN)

% Choose default command line output for gui_vis_filtered
handles.output = hObject;

% render the PropertyGrids in the correct panels
handles.propgrid = arg_guipanel(...
                        handles.pnlVisFiltered,     ...
                        'Function',@vis_filtered,   ...
                        'Parameters',{});
            
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_vis_filtered wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = gui_vis_filtered_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in cmdVisStream.
function cmdVisStream_Callback(hObject, eventdata, handles)
% hObject    handle to cmdVisStream (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
vis_filtered(handles.propgrid.GetPropertyValues);
