function varargout = gui_loadset(varargin)
% load a dataset via the GUI
% [Dataset,Action] = gui_loadset()
% [Dataset,Action] = gui_loadset('Initialpath',Initialpath)
%
% Out:
%   Dataset     : the dataset, an expression (or [], if the user pressed Cancel)
%   Action      : either 'OK' or 'Cancel'
%
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2010-10-27

% Last Modified by GUIDE v2.5 05-Jan-2012 07:14:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_loadset_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_loadset_OutputFcn, ...
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


% --- Executes just before gui_loadset is made visible.
function gui_loadset_OpeningFcn(hObject, eventdata, handles, varargin)
if length(varargin) < 2
    initialpath = env_translatepath('data:/');
else
    initialpath = env_translatepath(varargin{2});
end

% run file dialog...
[handles.filename,handles.filepath] = uigetfile({'*.*', 'any supported file'; '*.set', 'EEGLAB dataset'; '*.vhdr', 'BrainProducts Vision Recorder file'; ...
    '*.bdf', 'BioSEMI BDF file'; '*.cnt', 'Neuroscan CNT / ANT EEProbe CNT file'; '*.avr', 'ANT EEProbe average / Megis/BESA file'; '*.mat', 'MATLAB file'; ...
    '*.raw', 'Electrical Geodesics / ERPSS / Yokogawa raw data file '; '*.sna', 'Snapmaster SNA file'; '*.eeg', 'Neuroscan EEG / ANT EEProbe / BrainVision data file'; ...
    '*.edf', 'European data format file'; '*.ds', 'CTF folder'; '*.rdf', 'ERPSS data file'; '*.asc', 'INStep ASC file'; '*.m4d', '4D pdf file'; ...
    '*.smr', 'Cambridge Electronic Design file'; '*.egis', 'Electrical Geodesics file'; '*.ave', 'Yokogawa / Electrical Geodesics file'; ...
    '*.gave', 'Electrical Geodesics file'; '*.ses', 'Electrical Geodesics file '; '*.swf', 'Megis/BESA file'; '*.avg', 'NeuroScan file'; '*.nxe', 'Nexstim file'; ...
    '*.dat', 'BrainVision data file'; '*.seg', 'BrainVision data file'; '*.res4', 'VSM MedTech file'; '*.meg4', 'VSM MedTech file'; ...
    '*.m4d', 'Neuromag Elektra / 4D Neuroimaging file'; '*.pdf', 'Neuromag Elektra / 4D Neuroimaging file'; '*.xyz', 'Neuromag Elektra / 4D Neuroimaging file'},'Select dataset(s) to load',initialpath);

% Choose default command line output for gui_loadset
handles.output = {[],'Cancel'};
% Update handles structure
guidata(hObject, handles);

if ~isnumeric(handles.filename)
    % UIWAIT makes gui_loadset wait for user response (see UIRESUME)
    uiwait(handles.figure1);
else
    % otherwise we directly close the figure again...
    close(handles.figure1);
end


% --- Outputs from this function are returned to the command line.
function varargout = gui_loadset_OutputFcn(hObject, eventdata, handles) 
try 
    varargout = handles.output;
    close(handles.figure1);
catch
    varargout = {[],'Cancel'};
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% construct optional arguments
args = {};
channels = str2num(get(handles.edit2,'String'));
if ~isempty(channels)
    args = [args {'channels' channels}]; end
timerange = str2num(get(handles.edit5,'String'));
if ~isempty(timerange)
    args = [args {'timerange' timerange}]; end
samplerange = str2num(get(handles.edit4,'String'));
if ~isempty(samplerange)
    args = [args {'samplerange' samplerange}]; end
types = get(handles.edit8,'String');
if ~isempty(types)
    % neither quoted nor a cell array --> put in quotes
    if types(1) ~= '''' && types(1) ~= '{'        
        types = ['''' types '''']; end
    % evaluate & make a cell if necessary
    types = eval(types);
    if ~iscell(types)
        types = {types}; end
    args = [args {'types' types}]; 
end
miscopts = get(handles.edit6,'String');
if ~isempty(miscopts)
    args = [args evalin('base',['{' miscopts '}'])]; end
if isfield(handles,'markerchannel')
    args = [args {'markerchannel', {handles.markerchannel}}]; end
if ~isfield(handles,'insert_params')
    handles.insert_params = {}; end

% construct the loader expression
if iscell(handles.filename)
    for c=1:length(handles.filename)
        dataset{c} = io_loadset([handles.filepath,handles.filename{c}],args{:}); 
        for p=1:length(handles.insert_params)
            dataset{c} = set_insert_markers(dataset{c},handles.insert_params{p}); end
    
    end
    dataset = set_concat(dataset{:});
else
    dataset = io_loadset([handles.filepath,handles.filename],args{:});
    for p=1:length(handles.insert_params)
        dataset = set_insert_markers(dataset,handles.insert_params{p}); end
end

% optionally check if the dataset is loadable
if get(handles.checkbox1,'Value') == get(handles.checkbox1,'Max')
    tmp = exp_eval(dataset); end
        
% assign in workspace
assignin('base',get(handles.edit1,'String'),dataset);
handles.output = {dataset,'OK'};
guidata(hObject,handles);
uiresume(handles.figure1);


function pushbutton2_Callback(hObject, eventdata, handles)
handles.output = {[],'Cancel'};
guidata(hObject,handles);
uiresume(handles.figure1);

function pushbutton3_Callback(hObject, eventdata, handles)
handles.insert_params = gui_eventinsertion();
guidata(hObject,handles);

function pushbutton6_Callback(hObject, eventdata, handles)
tmp = arg_guidialog(@set_infer_markers,'Invoke',false);
if ~isempty(tmp)
    handles.markerchannel = tmp; end
guidata(hObject,handles);

function edit1_Callback(hObject, eventdata, handles)
gui_check_varname(hObject,'lastdata');

function edit2_Callback(hObject, eventdata, handles)
gui_check_vector(hObject);

function edit3_Callback(hObject, eventdata, handles)
gui_check_scalar(hObject);

function edit4_Callback(hObject, eventdata, handles)
gui_check_range(hObject);

function edit5_Callback(hObject, eventdata, handles)
gui_check_range(hObject);

function edit8_Callback(hObject, eventdata, handles)

function pushbutton5_Callback(hObject, eventdata, handles)
% help button
env_doc io_loadset

function figure1_KeyPressFcn(hObject, eventdata, handles)
if strcmp(eventdata.Key,'return')
    pushbutton1_Callback(handles.pushbutton1, eventdata, handles); end
if strcmp(eventdata.Key,'escape')
    pushbutton2_Callback(handles.pushbutton2, eventdata, handles); end


% --- auto-generated trash ---

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
function edit2_CreateFcn(hObject, eventdata, handles)
function edit3_CreateFcn(hObject, eventdata, handles)
function edit4_CreateFcn(hObject, eventdata, handles)
function edit5_CreateFcn(hObject, eventdata, handles)
function edit6_Callback(hObject, eventdata, handles)
function edit6_CreateFcn(hObject, eventdata, handles)
function checkbox1_Callback(hObject, eventdata, handles)
function edit8_CreateFcn(hObject, eventdata, handles)
function figure1_CreateFcn(hObject, eventdata, handles)

