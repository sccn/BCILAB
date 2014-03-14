function varargout = gui_applymodel(varargin)
% GUI_APPLYMODEL M-file for gui_applymodel.fig
%      GUI_APPLYMODEL, by itself, creates a new GUI_APPLYMODEL or raises the existing
%      singleton*.
%
%      H = GUI_APPLYMODEL returns the handle to a new GUI_APPLYMODEL or the handle to
%      the existing singleton*.
%
%      GUI_APPLYMODEL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_APPLYMODEL.M with the given input arguments.
%
%      GUI_APPLYMODEL('Property','Value',...) creates a new GUI_APPLYMODEL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_applymodel_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_applymodel_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_applymodel

% Last Modified by GUIDE v2.5 12-Dec-2011 03:25:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @gui_applymodel_OpeningFcn, ...
    'gui_OutputFcn',  @gui_applymodel_OutputFcn, ...
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


% --- Executes just before gui_applymodel is made visible.
function gui_applymodel_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = {[],[]};
guidata(hObject, handles);


% populate model selector listbox with admissible models
vars = evalin('base','whos');
vars = {vars(strcmp({vars.class},'struct')).name};
model_vars = {};
model_times = [];
% for all struct vars...
for v=vars
    % select those that are results
    var = evalin('base',v{1});
    if isfield(var,'tracking') && isfield(var.tracking,'prediction_function')        
        model_vars{end+1} = v{1}; 
        if isfield(var,'timestamp')
            model_times(end+1) = var.timestamp;
        else
            model_times(end+1) = 0;
        end
    end
end

if isempty(model_vars)
    % disable the OK button
    set(handles.popupmenu2,'String','(no model available)');
    set(handles.pushbutton1,'Enable','off');
else    
    % add to pulldown menu...
    set(handles.popupmenu2,'String',model_vars);
    % and select the most recent one
    set(handles.popupmenu2,'Value',argmax(model_times));
end
guidata(hObject, handles);

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
% invoke the loss selection menu, too
popupmenu3_Callback(handles.popupmenu3,{},guidata(hObject));

% UIWAIT makes gui_calibratemodel wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- OK button
function pushbutton1_Callback(hObject, eventdata, handles)
    
% get dataset
dataset = evalin('base',strtok(handles.dataset));
% chosen model
models = get(handles.popupmenu2,'String');
model = evalin('base',models{get(handles.popupmenu2,'Value')});
% get loss metric
metric = handles.metric;

% now start the computation
disp('beginning evaluation...');
[pred,loss,stats,targ] = bci_predict('model',model,'data',dataset,'metric',metric); %#ok<NASGU>

% write back the results
assignin('base',get(handles.edit1,'String'),stats);
handles.output = {stats};
guidata(hObject,handles);
uiresume(handles.figure1);

% display statistics
gui_reviewresults(stats);



% --- Outputs from this function are returned to the command line.
function varargout = gui_applymodel_OutputFcn(hObject, eventdata, handles)
try
    varargout = handles.output;
    try 
        close(handles.figure1);
    catch,end    
catch
    % cancelled
    varargout = {[],[]};
end



% --- cancel button
function pushbutton2_Callback(hObject, eventdata, handles)
uiresume(handles.figure1);


% --- dataset popup menu

function popupmenu1_Callback(hObject, eventdata, handles)
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
set(handles.popupmenu1,'String',printed);


% --- loss selector menu

function popupmenu3_Callback(hObject, eventdata, handles)
shortlosses = {'auto', 'kld', 'nll', 'mcr', 'mae', 'mse', 'smse', 'max', 'rms', 'bias', 'medse', 'auc', 'cond_entropy', 'cross_entropy', 'f_measure'};
handles.metric = shortlosses{get(hObject,'Value')};
guidata(hObject,handles);

function popupmenu3_CreateFcn(hObject, eventdata, handles)
set(hObject,'String',{'Automatically chosen', 'Kullback-Leibler divergence', 'Negative log-likelihood (NLL)', 'Mis-classification rate (MCR)', 'Mean absolute error (MAE)', ...
    'Mean square error (MSE)', 'Standardized mean square error (SMSE)', 'Maximum absolute error', 'Root mean square error (RMSE)', 'Directed mean bias', 'Median square error', 'Area under ROC (AUC)', ...
    'Conditional Entropy', 'Cross-entropy', 'negative F-Score'});


% --- results selector edit

function edit1_Callback(hObject, eventdata, handles)
gui_check_varname(hObject,'lastresults');

% help button
function pushbutton3_Callback(hObject, eventdata, handles)
env_doc bci_predict


% --- key events ---

function figure1_KeyPressFcn(hObject, eventdata, handles)
if strcmp(eventdata.Key,'return')
    pushbutton1_Callback(handles.pushbutton1, eventdata, handles); end
if strcmp(eventdata.Key,'escape')
    pushbutton2_Callback(handles.pushbutton2, eventdata, handles); end

% ---- auto-generated trash ----

function edit1_CreateFcn(hObject, eventdata, handles)
function popupmenu1_CreateFcn(hObject, eventdata, handles)
function popupmenu2_CreateFcn(hObject, eventdata, handles)
function popupmenu2_Callback(hObject, eventdata, handles)
