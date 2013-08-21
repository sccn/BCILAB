function varargout = gui_saveapproach(varargin)
% Save the given approach in the workspace and optionally on disk
% gui_saveapproach(approach)
%
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2010-10-28

% Last Modified by GUIDE v2.5 28-Mar-2012 23:28:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_saveapproach_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_saveapproach_OutputFcn, ...
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


% --- Executes just before gui_saveapproach is made visible.
function gui_saveapproach_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_saveapproach (see VARARGIN)
%            the first argument is the approach data structure...

if isempty(varargin)
    varargin{1} = gui_chooseapproach(); 
    if isempty(varargin{1})
        try
            close(handles.figure1);
        catch,end
        return; 
    end
end
    
handles.approach = varargin{1};
if ~isfield(handles.approach,'name') || isempty(handles.approach.name)
    % propose a new name
    handles.approach.name = 'BCI Approach';
else
    % create a new name based on an existing names
    xp = hlp_split(handles.approach.name,' ');
    try
        if length(xp) > 1
            n = sscanf(xp{end},'%f');
        else
            n = [];
        end
    catch
        n = [];
    end
    if ~isempty(n) && ~isnan(n)
        handles.approach.name = [sprintf('%s ',xp{1:end-1}) num2str(n+1)];
    else
        handles.approach.name = [handles.approach.name ' 1'];
    end
end
if ~isfield(handles.approach,'description')
    handles.approach.description = ''; end

% Update handles structure
guidata(hObject, handles);

% update the two edit fields
set(handles.edit2,'String',handles.approach.description);
set(handles.edit3,'String',handles.approach.name);

% UIWAIT makes gui_saveapproach wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_saveapproach_OutputFcn(hObject, eventdata, handles) 
try
    % Get default command line output from handles structure
    varargout{1} = handles.approach;
    % close the figure
    close(handles.figure1);
catch
    varargout{1} = [];
end


function edit1_Callback(hObject, eventdata, handles)
gui_check_varname(hObject,'lastapproach');


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% save in workspace...
handles.approach.timestamp = now;
assignin('base',get(handles.edit1,'String'),handles.approach);
uiresume(handles.figure1);


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% Save on disk: create a file chooser dialog
handles.approach.timestamp = now;
io_mkdirs('home:/.bcilab/approaches/');
[FileName,PathName] = uiputfile('*.apr','Save approach',env_translatepath('home:/.bcilab/approaches/untitled.apr'));
if FileName
    % ... and save the file
    identifier = get(handles.edit1,'String');
    % assign as struct field
    tmp.(identifier) = handles.approach;
    % and then save it from there...
    io_save([PathName,FileName],'-struct tmp',identifier); 
end


function edit2_Callback(hObject, eventdata, handles)
handles.approach.description = get(hObject,'String');
guidata(hObject,handles);

function edit3_Callback(hObject, eventdata, handles)
handles.approach.name = get(hObject,'String');
guidata(hObject,handles);

function figure1_KeyPressFcn(hObject, eventdata, handles)
if strcmp(eventdata.Key,'return')
    pushbutton1_Callback(handles.pushbutton1, eventdata, handles); end

% --- auto-generated trash ---

% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
function edit2_CreateFcn(hObject, eventdata, handles)
function edit1_CreateFcn(hObject, eventdata, handles)
function edit1_KeyPressFcn(hObject, eventdata, handles)


function pushbutton2_ButtonDownFcn(hObject, eventdata, handles)
