function varargout = gui_calibratemodel(varargin)
% Calibrate a model
% [Model, Stats] = gui_calibratemodel()
%
% Out:
%   Model : the calibrated model; [] if the process was cancelled
%   
%   Stats : calibration statistics, including performance estimate if selected; [] if the process was cancelled
%
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2010-10-25

% Last Modified by GUIDE v2.5 02-Aug-2012 16:32:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_calibratemodel_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_calibratemodel_OutputFcn, ...
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



% --- Executes just before gui_calibratemodel is made visible.
function gui_calibratemodel_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = {[],[]};
guidata(hObject, handles);

% update the popup contents (and approaches / datasets variables)
update_approaches(hObject);
update_datasets(hObject);

% select approach
handles = guidata(hObject);
if isempty(handles.approaches)
    warning off MATLAB:hg:uicontrol:ParameterValuesMustBeValid
    % no approach specified: create a new one
    approach = gui_newapproach;
    % cancelled: cancel this dialog, too
    if isempty(approach)
        return; end
    % otherwise (re-)create the approach list
    update_approaches(hObject);
end

% preferably select the lastapproach in the pulldown menus
handles = guidata(hObject);
guidata(hObject, handles);

% invoke the popup menu to propagate selection
popupmenu1_Callback(handles.popupmenu1,{},guidata(hObject));


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
dataids = strtok(cellstr(get(handles.popupmenu5,'String')));
selected_id = find(strcmp(dataids,'lastdata'),1);
if isempty(selected_id)
    selected_id = 1; end
set(handles.popupmenu5,'Value',selected_id);
guidata(hObject, handles);
% invoke the popup menu to propagate selection
popupmenu5_Callback(handles.popupmenu5,{},guidata(hObject));
% invoke the loss selection menu, too
popupmenu2_Callback(handles.popupmenu2,{},guidata(hObject));

% UIWAIT makes gui_calibratemodel wait for user response (see UIRESUME)
uiwait(handles.figure1);



% --- Outputs from this function are returned to the command line.
function varargout = gui_calibratemodel_OutputFcn(hObject, eventdata, handles) 
try
    varargout = handles.output;
    try 
        close(handles.figure1);
    catch,end    
catch
    % cancelled
    varargout = {[],[]};
end


% --- OK button 

function pushbutton1_Callback(hObject, eventdata, handles)
% start the computation
global tracking;
tracking.gui.last_markers = get(handles.edit8,'String');
tracking.gui.last_eventfield = get(handles.targetfieldlist, 'String');

% get opt scheme
opt_foldness = str2num(get(handles.edit1,'String'));
opt_margin = str2num(get(handles.edit2,'String'));
opt_scheme = {'chron' opt_foldness opt_margin};
% get eval scheme
if get(handles.checkbox1,'Value') == get(handles.checkbox1,'Max')
    eval_foldness = str2num(get(handles.edit3,'String'));
    eval_margin = str2num(get(handles.edit4,'String'));
    eval_scheme = {'chron' eval_foldness eval_margin};
else
    eval_scheme = 0;
end
% crossval compute engine
if get(handles.checkbox2,'Value') == get(handles.checkbox2,'Max')
    engine_cv = 'global';
else
    engine_cv = 'local';
end
tracking.gui.cluster_checked = get(handles.checkbox2,'Value');
% get nodepool
nodepool = get(handles.edit5,'String');
if isempty(nodepool) || strcmp(nodepool,'(use current config)')
    nodepool = 'global';
else
    try
        nodepool = eval(nodepool);
    catch
        disp('Error interpreting Node pool entry; falling back to current config');
        nodepool = 'global';
    end
end
% collect computation-relevant data
parameters = {'Data',evalin('base',strtok(handles.dataset)), 'Approach',handles.approach, 'EvaluationMetric', handles.metric, 'OptimizationScheme',opt_scheme, ...
    'EvaluationScheme',eval_scheme, 'CrossvalidationResources',engine_cv, 'ResourcePool',nodepool,'TargetMarkers',eval(get(handles.edit8,'String')), ...
    'EventField', eval(get(handles.targetfieldlist,'String'))};
% now start the computation
fprintf('\nbeginning new computation...\n');
try
    [loss,mdl,stats] = bci_train(parameters{:}); %#ok<ASGLU>
catch e
    fprintf('\nComputation failed; error trace:\n');
    hlp_handleerror(e);
    errordlg2([hlp_handleerror(e,0,false) sprintf('\n') '(See the MATLAB command window for a detailed stack trace.)'],'Error during model calibration');
    return;
end
    
assignin('base',get(handles.edit6,'String'),mdl);
assignin('base',get(handles.edit7,'String'),stats);
handles.output = {mdl, stats};
guidata(hObject,handles);
uiresume(handles.figure1);
% if cross-validation was checked, display statistics
if ~isequal(eval_scheme,0)
    gui_reviewresults(stats); end


% --- approaches & datasets popups

function popupmenu1_Callback(hObject, eventdata, handles)
% select an approach
handles.approach = handles.approaches{get(hObject,'Value')};
guidata(hObject,handles);

function update_approaches(hObject)
handles = guidata(hObject);
% retrieve and attach the known workspace approaches
approaches = {};
names = {};
approach_times = [];
vars = evalin('base','whos');
for v=find(strcmp({vars.class},'struct'))
    var = evalin('base',vars(v).name);
    if isfield(var,{'paradigm','parameters'})
        % list it
        approaches{end+1} = var;
        if isfield(var,'timestamp')
            approach_times(end+1) = var.timestamp;
        else
            approach_times(end+1) = 0;
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
handles.approach_times = approach_times;
guidata(hObject,handles);
set(handles.popupmenu1,'String',handles.approachnames);
set(handles.popupmenu1,'Value',argmax(approach_times));


function popupmenu5_Callback(hObject, eventdata, handles)
% make sure that we don't select [...] entries
idx = gui_select_bracketed(hObject);
% store the current approach in handles
handles.dataset = handles.indexable_datasets{idx};
guidata(hObject,handles);

function update_datasets(hObject)
handles = guidata(hObject);
% retrieve and attach the known datasets...
[handles.datasets,printed,handles.indexable_datasets] = gui_listdatasets;
guidata(hObject,handles);
set(handles.popupmenu5,'String',printed);


% --- loss popup 

function popupmenu2_Callback(hObject, eventdata, handles)
shortlosses = {'auto', 'kld', 'nll', 'mcr', 'mae', 'mse', 'smse', 'max', 'rms', 'bias', 'medse', 'auc', 'cond_entropy', 'cross_entropy', 'f_measure'};
handles.metric = shortlosses{get(hObject,'Value')};
guidata(hObject,handles);

function popupmenu2_CreateFcn(hObject, eventdata, handles)
set(hObject,'String',{'Automatically chosen', 'Kullback-Leibler divergence', 'Negative log-likelihood (NLL)', 'Mis-classification rate (MCR)', 'Mean absolute error (MAE)', ...
       'Mean square error (MSE)', 'Standardized mean square error (SMSE)', 'Maximum absolute error', 'Root mean square error (RMSE)', 'Directed mean bias', 'Median square error', 'Area under ROC (AUC)', ...
       'Conditional Entropy', 'Cross-entropy', 'negative F-Score'});

   
% --- edit fields and checkboxes...

function edit5_Callback(hObject, eventdata, handles)
% Node pool edit
str = get(hObject,'String');
if ~isempty(str) && ~strcmp('(use current config)',str)
    try
        assert(iscellstr(eval(str)),'Error!');
    catch
        errordlg2(sprintf('Must be a cell array of strings - example:\n{''computer1:25624'',''computer2:13456'',''computer3:34553''}'));
        set(hObject,'String','(use current config)');
    end
end

function edit5_CreateFcn(hObject, eventdata, handles)

function edit1_Callback(hObject, eventdata, handles)
gui_check_scalar(hObject,'5');

function edit2_Callback(hObject, eventdata, handles)
gui_check_scalar(hObject,'5');

function edit3_Callback(hObject, eventdata, handles)
gui_check_scalar(hObject,'10');

function edit4_Callback(hObject, eventdata, handles)
gui_check_scalar(hObject,'5');

function edit6_Callback(hObject, eventdata, handles)
gui_check_varname(hObject,'lastmodel');

function edit7_Callback(hObject, eventdata, handles)
gui_check_varname(hObject,'laststats');

function pushbutton3_Callback(hObject, eventdata, handles)
% help button
env_doc bci_train

function pushbutton4_Callback(hObject, eventdata, handles)
% "Inspect data..." button
try
    data = evalin('base',strtok(handles.dataset));
catch e
    env_handleerror(e);
    errordlg2(sprintf([e.message '\n(See the MATLAB command window for a detailed stack trace.)']),'Cannot find the data in the MATLAB workspace.');
    return
end
try
    data = exp_eval_optimized(data);
catch e
    env_handleerror(e);
    errordlg2(sprintf([e.message '\n(See the MATLAB command window for a detailed stack trace.)']),'Cannot open the data file.');
    return
end
try
    vis_artifacts(data,data,'ShowEventLegend',true);
catch e
    hlp_handleerror(e);
    errordlg2(sprintf([e.message '\n(See the MATLAB command window for a detailed stack trace.)']),'Cannot open the data viewer.');
end
    

function pushbutton2_Callback(hObject, eventdata, handles)
% cancel button
uiresume(handles.figure1);

function figure1_KeyPressFcn(hObject, eventdata, handles)
if strcmp(eventdata.Key,'return')
    pushbutton1_Callback(handles.pushbutton1, eventdata, handles); end
if strcmp(eventdata.Key,'escape')
    pushbutton2_Callback(handles.pushbutton2, eventdata, handles); end

function edit8_Callback(hObject, eventdata, handles)
global tracking;
try
    fallback = tracking.gui.last_markers;
catch
    fallback = '{''1'',''2''}';
end
gui_check_cell(hObject,fallback);

function targetfieldlist_Callback(hObject, eventdata, handles)
global tracking;
try
    fallback = tracking.gui.last_eventfield;
catch
    fallback = '''type''';
end
gui_check_cellstr_or_str(hObject,fallback);




function edit8_KeyPressFcn(hObject, eventdata, handles)

% --- auto-generated stuff...

function figure1_CreateFcn(hObject, eventdata, handles)

function popupmenu1_CreateFcn(hObject, eventdata, handles)

function checkbox1_Callback(hObject, eventdata, handles)

function checkbox2_Callback(hObject, eventdata, handles)

function edit3_CreateFcn(hObject, eventdata, handles)

function edit4_CreateFcn(hObject, eventdata, handles)

function edit1_CreateFcn(hObject, eventdata, handles)

function edit2_CreateFcn(hObject, eventdata, handles)

function popupmenu5_CreateFcn(hObject, eventdata, handles)

function edit6_CreateFcn(hObject, eventdata, handles)

function edit7_CreateFcn(hObject, eventdata, handles)



function edit8_CreateFcn(hObject, eventdata, handles)
global tracking;
try
    set(hObject,'String',tracking.gui.last_markers);
catch
end


function checkbox2_CreateFcn(hObject, eventdata, handles)
global tracking;
try
    set(hObject,'Value',tracking.gui.cluster_checked);
catch
end

% --- Executes during object creation, after setting all properties.
function targetfieldlist_CreateFcn(hObject, eventdata, handles)
global tracking;
try
    set(hObject,'String',tracking.gui.last_eventfield);
catch
end
