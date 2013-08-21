function varargout = gui_preprocessdata(varargin)
% GUI_PREPROCESSDATA M-file for gui_preprocessdata.fig
%      GUI_PREPROCESSDATA, by itself, creates a new GUI_PREPROCESSDATA or raises the existing
%      singleton*.
%
%      H = GUI_PREPROCESSDATA returns the handle to a new GUI_PREPROCESSDATA or the handle to
%      the existing singleton*.
%
%      GUI_PREPROCESSDATA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_PREPROCESSDATA.M with the given input arguments.
%
%      GUI_PREPROCESSDATA('Property','Value',...) creates a new GUI_PREPROCESSDATA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_preprocessdata_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_preprocessdata_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_preprocessdata

% Last Modified by GUIDE v2.5 09-Nov-2010 16:05:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_preprocessdata_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_preprocessdata_OutputFcn, ...
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


% --- Executes just before gui_preprocessdata is made visible.
function gui_preprocessdata_OpeningFcn(hObject, eventdata, handles, varargin)
% populate listbox with all applicable models...
vars = evalin('base','whos');
vars = {vars(strcmp({vars.class},'struct')).name};
model_vars = {};
% for all struct vars...
for v=vars
    % select those that are results
    var = evalin('base',v{1});
    if isfield(var,'tracking') && isfield(var.tracking,'prediction_function')
        model_vars{end+1} = v{1}; end
end

if isempty(model_vars)
    % disable the OK button
    set(handles.popupmenu1,'String','(no models available)');
    set(handles.pushbutton1,'Enable','off');
else
    % add to pulldown menu...
    set(handles.popupmenu1,'String',model_vars);
end

guidata(hObject, handles);

% update the popup contents (and approaches / datasets variables)
update_datasets(hObject);



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
dataids = strtok(cellstr(get(handles.popupmenu2,'String')));
selected_id = find(strcmp(dataids,'lastdata'),1);
if isempty(selected_id)
    selected_id = 1; end
set(handles.popupmenu2,'Value',selected_id);
guidata(hObject, handles);
% invoke the popup menu to propagate selection
popupmenu2_Callback(handles.popupmenu2,{},guidata(hObject));





% Choose default command line output for gui_preprocessdata
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_preprocessdata wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_preprocessdata_OutputFcn(hObject, eventdata, handles) 
try
    varargout{1} = handles.output;
    close(handles.figure1);
catch
end



% OK button
function pushbutton1_Callback(hObject, eventdata, handles)
% get the selection...
mdcontents = get(handles.popupmenu1,'String');
mdselection = get(handles.popupmenu1,'Value');
dscontents = get(handles.popupmenu2,'String');
dsselection = get(handles.popupmenu2,'Value');
uiresume(handles.figure1);
EEG = bci_preproc(evalin('base',strtok(dscontents{dsselection})),evalin('base',mdcontents{mdselection}));
% write back...
assignin('base',get(handles.edit2,'String'),EEG);
% re-load if desired...
if get(handles.checkbox1,'Value') == get(handles.checkbox1,'Max')
    eeglab redraw; end
    
% Cancel button
function pushbutton2_Callback(hObject, eventdata, handles)
uiresume(handles.figure1);


% --- dataset chooser
function update_datasets(hObject)
handles = guidata(hObject);
% retrieve and attach the known datasets...
[handles.datasets,printed,handles.indexable_datasets] = gui_listdatasets;
guidata(hObject,handles);
set(handles.popupmenu2,'String',printed);


function popupmenu2_Callback(hObject, eventdata, handles)
% make sure that we don't select [...] entries
idx = gui_select_bracketed(hObject);
% store the current approach in handles
handles.dataset = handles.indexable_datasets{idx};
guidata(hObject,handles);


function edit2_Callback(hObject, eventdata, handles)
gui_check_varname(hObject,'EEG');
if strcmp(get(hObject,'String'),{'EEG','ALLEEG'})
   set(handles.checkbox1,'Enable','on');
else
   set(handles.checkbox1,'Enable','off');    
end
   
% --- auto-generated junk 
function popupmenu1_Callback(hObject, eventdata, handles)
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function popupmenu2_CreateFcn(hObject, eventdata, handles)

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function checkbox1_Callback(hObject, eventdata, handles)


