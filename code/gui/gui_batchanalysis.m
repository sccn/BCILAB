function varargout = gui_batchanalysis(varargin)
% GUI_BATCHANALYSIS MATLAB code for gui_batchanalysis.fig
%      GUI_BATCHANALYSIS, by itself, creates a new GUI_BATCHANALYSIS or raises the existing
%      singleton*.
%
%      H = GUI_BATCHANALYSIS returns the handle to a new GUI_BATCHANALYSIS or the handle to
%      the existing singleton*.
%
%      GUI_BATCHANALYSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_BATCHANALYSIS.M with the given input arguments.
%
%      GUI_BATCHANALYSIS('Property','Value',...) creates a new GUI_BATCHANALYSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_batchanalysis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_batchanalysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_batchanalysis

% Last Modified by GUIDE v2.5 11-Dec-2011 22:25:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_batchanalysis_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_batchanalysis_OutputFcn, ...
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


% --- Executes just before gui_batchanalysis is made visible.
function gui_batchanalysis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_batchanalysis (see VARARGIN)

% populate some of the fields with defaults of bci_batchtrain
defaults = arg_report('vals',@bci_batchtrain);
if ~isempty(defaults.loadargs)
    set(handles.edit7,'String',hlp_tostring(defaults.loadargs)); end
if ~isempty(defaults.trainargs)
    set(handles.edit8,'String',hlp_tostring(defaults.trainargs)); end
if ~isempty(defaults.predictargs)
    set(handles.edit9,'String',hlp_tostring(defaults.predictargs)); end
if ~isempty(defaults.saveargs)
    set(handles.edit10,'String',hlp_tostring(defaults.saveargs)); end
set(handles.edit11,'String',defaults.storepatt);
set(handles.edit12,'String',defaults.resultpatt);

handles.output = [];



% Update handles structure
guidata(hObject, handles);




% --- Outputs from this function are returned to the command line.
function varargout = gui_batchanalysis_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global tracking;

% turn all large edit panes to multi-line editors
for obj = {handles.edit1,handles.edit20,handles.edit19}
    jScrollPane = findjobj(obj{1});
    jViewPort = jScrollPane.getViewport;
    jEditbox = jViewPort.getComponent(0);
    jEditbox.setWrapping(false);
    set(jScrollPane,'HorizontalScrollBarPolicy',30);
end

% get all paradigm names
para_files = [];
para_paths = {'functions:/paradigms', 'functions:/paradigms/in_development', 'home:/.bcilab/code/paradigms'};
if ~isempty(tracking.paths.private_path)
    para_paths = [para_paths {'private:/code/paradigms','private:/code/paradigms/in_development'}]; end
for p = para_paths
    para_files = [para_files; dir([env_translatepath(p{1}) filesep 'Paradigm*.m'])]; end
para_acronyms = cellfun(@(s)s(9:end-2),{para_files.name},'UniformOutput',false);

% populate the list of approaches with available in the workspace
approaches = {};
vars = evalin('base','whos');
for v=find(strcmp({vars.class},'struct'))
    var = evalin('base',vars(v).name);
    if isfield(var,{'paradigm','parameters'})
        approaches{end+1} = vars(v).name; end
end
for v=find(strcmp({vars.class},'cell'))
    var = evalin('base',vars(v).name);
    if length(var) >= 1 && ischar(var{1}) && any(strcmp(var{1},para_acronyms))
        approaches{end+1} = vars(v).name; end
end
set(handles.edit20,'String',char(approaches));

traindata = {};
vars = evalin('base','whos');
for v=find(strcmp({vars.class},'struct'))
    var = evalin('base',vars(v).name);
    if all(isfield(var,{'head','parts'})) && strcmp(char(var.head),'io_loadset')
        traindata{end+1} = vars(v).name; end
end
set(handles.edit1,'String',char(traindata));

% also set up the targetmarkers
global tracking;
try
    set(handles.edit6,'String',tracking.gui.last_markers);
catch
end

% wait for further input
uiwait(handles.figure1);

% Get default command line output from handles structure
try
    varargout{1} = handles.output;
    close(handles.figure1);
catch
    varargout = {};
end


function pushbutton1_Callback(hObject, eventdata, handles)

% --- process approaches, training sets, and test sets ---

% get the approaches, traindata and testdata cell arrays
approaches = cellstr(get(handles.edit20,'String'));
traindata = cellstr(get(handles.edit1,'String'));
testdata = cellstr(get(handles.edit19,'String'));
% prune empty entries
approaches = approaches(~cellfun('isempty',approaches));
trainsets = traindata(~cellfun('isempty',traindata));
testsets = testdata(~cellfun('isempty',testdata));

problematic = {};
for k=1:length(approaches)
    if ~isvarname(approaches{k})
        problematic{end+1} = approaches{k};
    else
        try
            if ~evalin('base',sprintf('exist(''%s'',''var'')',approaches{k}))
                problematic{end+1} = approaches{k}; end
        catch
            problematic{end+1} = approaches{k};
        end
    end
end
if ~isempty(problematic)
    warndlg2(sprintf(['The following approaches do not exist in the workspace:\n' hlp_tostring(problematic)]),'Non-existing approaches.'); 
    return
end

if isempty(approaches)
    warndlg2('You need to specify at least one approach. By default, any approach that you have created (using "New Approach" from the menu) will be listed here.','No approach specified');
    return;
end
if isempty(trainsets)
    warndlg2('You need to specify at least one training set. By default, any set that you have loaded (using "Load recording(s)" from the menu) will be listed here.','No training data specified');
    return;
end

% resolve all data references properly
trainsets = resolve_data(trainsets,'training sets');
if ~iscell(trainsets)
    return; end
testsets = resolve_data(testsets,'test sets');
if ~iscell(testsets)
    return; end

% --- process custom pipeline settings ---

% get target markers
targetmarkers = evalin('base',get(handles.edit6,'String'));
% get load arguments
loadargs = get_cell_args(handles.edit7,'load arguments');
if ~iscell(loadargs)
    return; end
% get train arguments
trainargs = get_cell_args(handles.edit8,'train arguments');
if ~iscell(trainargs)
    return; end
% get predict arguments
predictargs = get_cell_args(handles.edit9,'predict arguments');
if ~iscell(predictargs)
    return; end
% get save arguments
saveargs = get_cell_args(handles.edit10,'save arguments');
if ~iscell(saveargs)
    return; end

% --- process computation settings ---

% get reuse flag
reuseexisting = get(handles.checkbox1,'Value') == get(handles.checkbox1,'Max');

% get cluster resources flag
if get(handles.checkbox2,'Value') == get(handles.checkbox2,'Max')
    clusterresources = 'global';
else
    clusterresources = 'local';
end

% get worker pool
workerpool = get(handles.edit18,'String');
if strcmp(workerpool,'(use global settings)')
    workerpool = 'global';
else
    try
        workerpool = evalin('base',workerpool);
        assert(iscellstr(workerpool));
    catch
        warndlg2('If not using global settings, the worker pool must be a cell array of ''host:port'' items.'); 
        return
    end
end


% --- process output formatting ---

% get study tag
studytag = get(handles.edit21,'String');
if ~isvarname(studytag)
    disp('It is recommended that the study identifier is a valid MATLAB variable name.'); end % perhaps for future compatibility...

% get the storage pattern
storagepatt = get(handles.edit11,'String');
resultpatt = get(handles.edit12,'String');

% invoke batch analysis
try
    disp('beginning computation...');
    results = bci_batchtrain('Datasets',trainsets,'Approaches',approaches,'PredictSets',testsets,'TargetMarkers',targetmarkers, ....
        'LoadArguments',loadargs,'TrainArguments',trainargs,'PredictArguments',predictargs,'SaveArguments',saveargs, ...
        'StudyTag',studytag,'StoragePattern',storagepatt,'ResultPattern',resultpatt,'ReuseExisting',reuseexisting, ...
        'ClusterResources',clusterresources,'ClusterPool',workerpool);
catch e
    env_handleerror(e);
    warndlg2(['Error while starting the batch analysis: ' e.message],'Error'); 
    return;
end
    
global tracking;
tracking.gui.last_markers = get(handles.edit6,'String');

% assign results into workspace variable
assignin('base',studytag,results);

% all went well - resume!
handles.output = results;
uiresume(handles.figure1);



function pushbutton2_Callback(hObject, eventdata, handles)
uiresume(handles.figure1);

function pushbutton3_Callback(hObject, eventdata, handles)
env_doc bci_batchtrain



    
function edit6_Callback(hObject, eventdata, handles)
global tracking;
try
    fallback = tracking.gui.last_markers;
catch
    fallback = '{''1'',''2''}';
end
gui_check_cell(hObject,fallback);


function pushbutton4_Callback(hObject, eventdata, handles)
% try to reuse the current arguments as a starting point...
loadargs = get_cell_args(handles.edit7,[]);
% make a dialog for this
loadargs = arg_guidialog(@io_loadset,'Parameters',loadargs,'Title','Configure load arguments','Invoke',false);
% and assign to the edit field, if not empty
if ~isempty(loadargs)
    set(handles.edit7,'String',hlp_tostring(hlp_struct2varargin(loadargs))); end


function pushbutton5_Callback(hObject, eventdata, handles)
% try to reuse the current arguments as a starting point...
trainargs = get_cell_args(handles.edit8,[]);
% make a dialog for this
trainargs = arg_guidialog(@bci_train,'Parameters',trainargs,'Title','Configure train arguments','Invoke',false);
% and assign to the edit field, if not empty
if ~isempty(trainargs)
    set(handles.edit8,'String',hlp_tostring(hlp_struct2varargin(trainargs))); end


function pushbutton6_Callback(hObject, eventdata, handles)
% try to reuse the current arguments as a starting point...
predictargs = get_cell_args(handles.edit9,[]);
% make a dialog for this
predictargs = arg_guidialog(@bci_predict,'Parameters',predictargs,'Title','Configure prediction arguments','Invoke',false);
% and assign to the edit field, if not empty
if ~isempty(predictargs)
    set(handles.edit9,'String',hlp_tostring(hlp_struct2varargin(predictargs))); end

function edit21_Callback(hObject, eventdata, handles)
gui_check_varname(hObject,'mystudy');

% --- utility functions ---


function data = resolve_data(data,name)
% resolve all references to data in the given cell array of strings
% allowed formats are:
% - varname or quoted string, or unquoted file name
% - cell array of varnames or quoted strings (but not of file names)
% for each cell...
try
    for c=1:length(data)
        data{c} = strtrim(data{c});
        % check if this is an unquoted file name reference
        if (data{c}(1) ~= '{' || ~data{c}(end) ~= '}') && (~data{c}(1) ~= '''' || ~data{c}(end) ~= '''')
            if exist(data{c},'file')
                % it's a file name: wrap in quotes
                data{c} = ['''' data{c} ''''];
            elseif (any(data{c}=='*') || any(data{c}=='?')) && ~isempty(rdir(data{c}))
                % it's apparently a file name pattern: wrap in quotes
                data{c} = ['''' data{c} ''''];
            elseif isvarname(data{c}) && evalin('base',sprintf('exist(''%s'',''var'')',data{c}))
                % assume it's a varname
            else
                error(['The follwoing unqouted data set refers neither to a file on disk nor a variable in the workspace: ' data{c}]);
            end
        end
        
        % try to evaluate it
        try
            data{c} = evalin('base',data{c});
        catch e
            env_handleerror(e);
            error(['The following data set reference cannot be resolved: ' data{c}]);
        end
    end
catch e
    warndlg2(sprintf('At least one of the %s is malformed:\n%s',name,e.message),'Mal-formed training set.');
    data = [];
end

function args = get_cell_args(handle,name)
% get a cell argument pack
args = get(handle,'String');
if ~isempty(args)
    if args(1) ~= '{' || args(end) ~= '}'
        args = ['{' args '}']; end
    try
        args = evalin('base',args);
    catch e
        if isempty(name)
            % revert empty cell
            args = {};
        else
            % return [] to indicate formal failure
            warndlg2(sprintf('The %s cannot be evaluated; error message: %s',name,e.message),'Mal-formed load arguments.'); 
            args = [];
        end
    end
else
    args = {};
end

% --- auto-generated stuff ---

function edit1_Callback(hObject, eventdata, handles)
function edit1_CreateFcn(hObject, eventdata, handles)
function edit2_Callback(hObject, eventdata, handles)
function edit2_CreateFcn(hObject, eventdata, handles)
function edit4_Callback(hObject, eventdata, handles)
function edit4_CreateFcn(hObject, eventdata, handles)
function edit20_Callback(hObject, eventdata, handles)
function edit20_CreateFcn(hObject, eventdata, handles)
function edit19_Callback(hObject, eventdata, handles)
function edit19_CreateFcn(hObject, eventdata, handles)
function checkbox1_Callback(hObject, eventdata, handles)
function checkbox2_Callback(hObject, eventdata, handles)
function edit18_Callback(hObject, eventdata, handles)
function edit18_CreateFcn(hObject, eventdata, handles)
function edit12_Callback(hObject, eventdata, handles)
function edit12_CreateFcn(hObject, eventdata, handles)
function edit11_Callback(hObject, eventdata, handles)
function edit11_CreateFcn(hObject, eventdata, handles)
function edit7_Callback(hObject, eventdata, handles)
function edit7_CreateFcn(hObject, eventdata, handles)
function edit8_Callback(hObject, eventdata, handles)
function edit8_CreateFcn(hObject, eventdata, handles)
function edit9_Callback(hObject, eventdata, handles)
function edit9_CreateFcn(hObject, eventdata, handles)
function edit10_Callback(hObject, eventdata, handles)
function edit10_CreateFcn(hObject, eventdata, handles)
function edit6_CreateFcn(hObject, eventdata, handles)
function edit21_CreateFcn(hObject, eventdata, handles)
