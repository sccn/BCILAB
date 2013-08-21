function varargout = gui_reviewresults(varargin)
% GUI_REVIEWRESULTS M-file for gui_reviewresults.fig
%      GUI_REVIEWRESULTS, by itself, creates a new GUI_REVIEWRESULTS or raises the existing
%      singleton*.
%
%      H = GUI_REVIEWRESULTS returns the handle to a new GUI_REVIEWRESULTS or the handle to
%      the existing singleton*.
%
%      GUI_REVIEWRESULTS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_REVIEWRESULTS.M with the given input arguments.
%
%      GUI_REVIEWRESULTS('Property','Value',...) creates a new GUI_REVIEWRESULTS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_reviewresults_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_reviewresults_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_reviewresults

% Last Modified by GUIDE v2.5 13-Apr-2012 20:57:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @gui_reviewresults_OpeningFcn, ...
    'gui_OutputFcn',  @gui_reviewresults_OutputFcn, ...
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


% --- Executes just before gui_reviewresults is made visible.
function gui_reviewresults_OpeningFcn(hObject, eventdata, handles, varargin)
% get data and options
data = varargin{1};
% set up default options for utl_summarize_structarray
opts = hlp_varargin2struct(varargin(2:end),'leftindent',[0 15],'rightindent',[0 15], ...
    'rewrite_fields',{'kld','KL divergence','nll','Negative log-likelihood', ...
                      'mcr','Error rate', 'mae','Mean abs error', 'mse','Mean square error',...
                      'max','Max absolut error','rms','Root mean square error','bias','Mean bias',...
                      'medse','Median square error','auc','Neg. area under ROC', ...
                      'cond_entropy','Conditional entropy','cross_entropy','Cross-Entropy', ...
                      'f_measure','Negative F-score','classes','Target classes','class_ratio','Class ratio', ...
                      'per_class','Per-class error','TP','True positive rate','TN','True negative rate',...
                      'FP','False positive rate','FN','False negative rate','measure','Measure',...
                      'targ','Targets','pred','Predictions','time','Runtime','target','Targets',...
                      'prediction','Predictions','expression','Reproducing expression','value','Value'});
                  
% pick out the correct part of the array
if isscalar(data)
    sub_arrays = {};
    % check if some of the sub-fields are a struct arrays...
    for f=fieldnames(data)'
        field = data.(f{1});
        if isstruct(field) && ~isscalar(field)
            sub_arrays{end+1} = field; end
    end
    if length(sub_arrays) == 1
        % exactly one sub-array: take this as the data to summarize...
        data = sub_arrays{1};
    end
end

if isfield(data,'timestamp')
    data = rmfield(data,'timestamp'); end

% set the summary text
[text,handles.data] = utl_summarize_structarray(data,opts);
set(handles.edit1,'String',[{''} text]);
set(handles.edit1,'Position',[0 0 1 1],'Units','normalized');

% create uitable and popuate with data
conts = handles.data.conts;
fullconts = handles.data.fullconts;
fnames = handles.data.fnames;
% no column is editable
columneditable = false(1,length(fnames));
% by default all columns are of type 'char'
columnformat = repmat({'char'},1,length(fnames));
for c=1:length(fullconts)
    col = fullconts{c};
    % check if all fields of a column are numeric or logical
    if all(cellfun(@isnumeric,conts{c}))
        columnformat{c} = 'numeric';
        for r=1:length(col)
            if isempty(col{r}) || ~isnumeric(col{r})
                col{r} = NaN; end
        end
    elseif all(cellfun(@islogical,conts{c}))
        columnformat{c} = 'logical';        
    else
        % otherwise we convert to char what is not yet char
        for r=find(~cellfun('isclass',col,'char'))
            col{r} = hlp_tostring(col{r}); end
    end
    % ... and put it into the table data
    tabledata(:,c) = col;
end
handles.tabledata = tabledata;
handles.table1 = uitable('Data',tabledata,'ColumnName',fnames,'Parent',handles.uipanel2, ...
    'Units','normalized','Position',[0 0 1 1],'ColumnEditable',columneditable,'ColumnFormat',columnformat);
set(handles.table1,'Units','pixels')
wh = get(handles.table1,'Position');
colwidth = (wh(3)-35) / length(fnames);
set(handles.table1,'Units','normalized')
set(handles.table1,'ColumnWidth',repmat({round(colwidth)},1,length(fnames)));


% Choose default command line output for gui_reviewresults
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);
% UIWAIT makes gui_reviewresults wait for user response (see UIRESUME)
uiwait(handles.figure1);



% --- Outputs from this function are returned to the command line.
function varargout = gui_reviewresults_OutputFcn(hObject, eventdata, handles)
try
    varargout{1} = handles.output;
    close(handles.figure1);
catch
end

% help button
function pushbutton5_Callback(hObject, eventdata, handles)
disp('coming soon...');

% ok button
function pushbutton1_Callback(hObject, eventdata, handles)
uiresume(handles.figure1);

function figure1_KeyPressFcn(hObject, eventdata, handles)
if strcmp(eventdata.Key,'return') || strcmp(eventdata.Key,'escape') 
    pushbutton1_Callback(handles.pushbutton1, eventdata, handles); end


% Explore button
function pushbutton4_Callback(hObject, eventdata, handles)
function pushbutton2_Callback(hObject, eventdata, handles)
function pushbutton3_Callback(hObject, eventdata, handles)
function figure1_CreateFcn(hObject, eventdata, handles)
function uipanel1_CreateFcn(hObject, eventdata, handles)
function edit1_Callback(hObject, eventdata, handles)
function edit1_CreateFcn(hObject, eventdata, handles)


function pushbutton6_Callback(hObject, eventdata, handles)
io_mkdirs('home:/.bcilab/results/');
[FileName,PathName] = uiputfile('*.csv','Save table',env_translatepath('home:/.bcilab/results/untitled.csv'));
if FileName
    csvwrite([PathName,FileName],handles.tabledata);  end


function pushbutton6_ButtonDownFcn(hObject, eventdata, handles)
