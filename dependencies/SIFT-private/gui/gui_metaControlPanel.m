function varargout = gui_metaControlPanel(varargin)
%
% GUI_METACONTROLPANEL M-file for gui_metaControlPanel.fig
%      GUI_METACONTROLPANEL, by itself, creates a new GUI_METACONTROLPANEL or raises the existing
%      singleton*.
%
%      H = GUI_METACONTROLPANEL returns the handle to a new GUI_METACONTROLPANEL or the handle to
%      the existing singleton*.
%
%      GUI_METACONTROLPANEL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_METACONTROLPANEL.M with the given input arguments.
%
%      GUI_METACONTROLPANEL('Property','Value',...) creates a new GUI_METACONTROLPANEL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_metaControlPanel_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_metaControlPanel_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_metaControlPanel

% Last Modified by GUIDE v2.5 11-Dec-2012 21:36:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_metaControlPanel_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_metaControlPanel_OutputFcn, ...
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


% --- Executes just before gui_metaControlPanel is made visible.
function gui_metaControlPanel_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_metaControlPanel (see VARARGIN)

% Choose default command line output for gui_metaControlPanel
handles.output = hObject;
handles.gcf    = hObject;
set(hObject,'Name','Meta Control Panel');
% set default termination behavior
handles.ExitButtonClicked = 'Cancel';

% extract some data from command-line input
if isempty(varargin)
    error('First input to gui_metaControlPanel must be a valid bcilab dataset');
end
if length(varargin)<2
    error('Second input to gui_metaControlPanel must be an opts structure');
end

% set default EEGLAB background and text colors
%-----------------------------------------------
try, icadefs;
catch,
	GUIBACKCOLOR        =  [.8 .8 .8];     
	GUIPOPBUTTONCOLOR   = [.8 .8 .8];    
	GUITEXTCOLOR        = [0 0 0];
end;

allhandlers = [hObject handles.pnlPropertyGridMisc];

hh = findobj(allhandlers,'style', 'text');
%set(hh, 'BackgroundColor', get(hObject, 'color'), 'horizontalalignment', 'left');
set(hh, 'Backgroundcolor', GUIBACKCOLOR);
set(hh, 'foregroundcolor', GUITEXTCOLOR);
set(hObject, 'color',GUIBACKCOLOR );
% set(hh, 'horizontalalignment', g.horizontalalignment);

hh = findobj(allhandlers, 'style', 'edit');
set(hh, 'BackgroundColor', [1 1 1]); %, 'horizontalalignment', 'right');

hh =findobj(allhandlers, 'style', 'pushbutton');
if ~strcmpi(computer, 'MAC') && ~strcmpi(computer, 'MACI') % this puts the wrong background on macs
    set(hh, 'backgroundcolor', GUIPOPBUTTONCOLOR);
    set(hh, 'foregroundcolor', GUITEXTCOLOR);
end;
hh =findobj(allhandlers, 'style', 'popupmenu');
set(hh, 'backgroundcolor', GUIPOPBUTTONCOLOR);
set(hh, 'foregroundcolor', GUITEXTCOLOR);
hh =findobj(allhandlers, 'style', 'checkbox');
set(hh, 'backgroundcolor', GUIBACKCOLOR);
set(hh, 'foregroundcolor', GUITEXTCOLOR);
hh =findobj(allhandlers, 'style', 'listbox');
set(hh, 'backgroundcolor', GUIPOPBUTTONCOLOR);
set(hh, 'foregroundcolor', GUITEXTCOLOR);
hh =findobj(allhandlers, 'style', 'radio');
set(hh, 'foregroundcolor', GUITEXTCOLOR);
set(hh, 'backgroundcolor', GUIPOPBUTTONCOLOR);
set(hObject, 'visible', 'on');

% set(handles.pnlPropertyGridFltPip,'backgroundcolor', GUIBACKCOLOR);
% set(handles.pnlPropertyGridFltPip,'foregroundcolor', GUITEXTCOLOR);
% 
% set(handles.pnlPropertyGridSiftPip,'backgroundcolor', GUIBACKCOLOR);
% set(handles.pnlPropertyGridSiftPip,'foregroundcolor', GUITEXTCOLOR);

set(handles.pnlPropertyGridMisc,'backgroundcolor', GUIBACKCOLOR);
set(handles.pnlPropertyGridMisc,'foregroundcolor', GUITEXTCOLOR);

%-----------------------------------------------

% store opts structure
handles.opts    = varargin{2};
handles.dataset = varargin{1};

handles = redrawPropertyGrids(hObject,handles);

% save figure handle
handles.figureHandle = hObject;

% Update handles structure
guidata(hObject, handles);


% Wait for user to click OK, Cancel or close figure
% uiwait(handles.gui_MetaControlPanel);


% UIWAIT makes gui_metaControlPanel wait for user response (see UIRESUME)
% uiwait(handles.gui_MetaControlPanel);


% --- User-defined function to (re)draw Property Grids
function handles = redrawPropertyGrids(hObject,handles)

try
    % delete existing property grids
    delete([handles.SiftPipPropertyGridHandle.Control, ...
            handles.FltPipPropertyGridHandle.Control, ....
            handles.MiscPropertyGridHandle.Control]);
catch err
end
 
handles = arg_setdirect(handles,false);

% render the PropertyGrids in the correct panels
handles.SiftPipPropertyGridHandle = arg_guipanel( ...
                handles.pnlPropertyGridSiftPip, ...
                'Function',@onl_siftpipeline, ...
                'Parameters',{'EEG',handles.dataset,handles.opts.siftPipCfg});
    
handles.FltPipPropertyGridHandle = arg_guipanel( ...
                 handles.pnlPropertyGridFltPip, ...
                'Function',@flt_pipeline, ...
                'Parameters',{'signal',handles.dataset,handles.opts.fltPipCfg});

handles.MiscPropertyGridHandle = arg_guipanel(...
                 handles.pnlPropertyGridMisc, ...
                 'Function',@onl_MetaCtrlPnlMiscOpts,'Parameters',{handles.opts.miscOptCfg});
             

% initialize some control contents from opts input
% try set(handles.txtWindowLength,'String',num2str(handles.opts.winLenSec));
% catch; end
% try set(handles.chkDispBenchmarking,'Value',handles.opts.doBenchmark);
% catch; end
% try set(handles.chkDispBrainMovie,'Value',handles.opts.doBrainMovie);
% catch; end
% try set(handles.chkDispSpectrum,'Value',handles.opts.dispSpectrum);
% catch; end
% try set(handles.chkAdaptBMLimits,'Value',handles.opts.adaptLimits);
% catch; end
% try set(handles.doSIFT,'Value',handles.opts.doSIFT);
% catch; end
drawnow

% --- Get contents of Property Grids and store in subfields of handles.opts
function handles = getPropertyGridContents(handles)
% handles  handles structure with opts subfield 

% get the property specifications
handles.opts.fltPipCfg  = arg_tovals(handles.FltPipPropertyGridHandle.GetPropertySpecification,true);
handles.opts.siftPipCfg  = arg_tovals(handles.SiftPipPropertyGridHandle.GetPropertySpecification,true);
handles.opts.miscOptCfg = arg_tovals(handles.MiscPropertyGridHandle.GetPropertySpecification,false);


% --- Outputs from this function are returned to the command line.
function varargout = gui_metaControlPanel_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

varargout = {hObject};


% --- Executes on button press in cmdCancel.
function cmdCancel_Callback(hObject, eventdata, handles)

if isempty(handles)
    % user closed the figure
    varargout = {[] hObject};
% elseif strcmpi(handles.ExitButtonClicked,'OK')
%     % user clicked OK
%     % get PropertySpecification
%     varargout = {handles.PropertyGridHandle handles.output};
else
    % user clicked cancel
    varargout = {[] handles.output};
end

handles.opts.exitPipeline = true;
assignin('base','opts',handles.opts);
guidata(hObject,handles);

try, 
    close(handles.figureHandle);
catch
end



% --- Executes on button press in cmdUpdate.
function cmdUpdate_Callback(hObject, eventdata, handles)

% handles.ExitButtonClicked ='OK';

% get the contents of the prop grids
% results will be stored in subfields of handles.opts
handles = getPropertyGridContents(handles);

% save pipeline configurations in base workspace and set new data flags
assignin('base','opts',handles.opts);
pause(0.1);
assignin('base','newPipeline',true);
set(handles.gcf,'UserData','initialized');
% evalin('base','opts.fltpipeline  = evalin(''caller'',''fltcfg'')');
% evalin('base','opts.siftpipeline = evalin(''caller'',''siftcfg'')');

% handles.opts.winLenSec = str2double(get(handles.txtWindowLength,'String'));

guidata(hObject,handles);
% guidata(hObject,handles);
% uiresume(handles.gui_MetaControlPanel);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADDITIONAL CALLBACKS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% % --- Executes during object creation, after setting all properties.
% function slideCurTime_CreateFcn(hObject, eventdata, handles)
% 
% 
% % Hint: slider controls usually have a light gray background.
% if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor',[.9 .9 .9]);
% end


% --- Executes during object deletion, before destroying properties.
function cmdUpdate_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to cmdUpdate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% 
% % --- Executes on button press in cmdHelp.
% function cmdHelp_Callback(hObject, eventdata, handles)
% % hObject    handle to cmdHelp (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% warndlg2('Coming soon!');


% --- Executes on button press in cmdPause.
function cmdPause_Callback(hObject, eventdata, handles)
% hObject    handle to cmdPause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% pause the pipeline
if strcmpi(get(hObject,'String'),'pause');
    handles.opts.holdPipeline = true;
    set(hObject,'String','Unpause');
else
    handles.opts.holdPipeline = false;
    set(hObject,'String','Pause');
end

assignin('base','opts',handles.opts);
guidata(hObject,handles);


% --- Executes during object deletion, before destroying properties.
function chkDispBenchmarking_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to chkDispBenchmarking (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in cmdVisStream.
function cmdVisStream_Callback(hObject, eventdata, handles)
% hObject    handle to cmdVisStream (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% open stream visualizer
vis_filtered;



% --- Executes when pnlPropertyGridFltPip is resized.
function pnlPropertyGridFltPip_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to pnlPropertyGridFltPip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when pnlPropertyGridSiftPip is resized.
function pnlPropertyGridSiftPip_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to pnlPropertyGridSiftPip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when pnlPropertyGridMisc is resized.
function pnlPropertyGridMisc_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to pnlPropertyGridMisc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in cmdSaveConfig.
function cmdSaveConfig_Callback(hObject, eventdata, handles)
% hObject    handle to cmdSaveConfig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get current contents of property grid
handles = getPropertyGridContents(handles);
opts    = handles.opts;

% select the path/file for saving configs
[fname fpath] = uiputfile('*.mat','Save Config File');
if ~fname
    return;
end
% save the opts structure
save(fullfile(fpath,fname),'opts');


% --------------------------------------------------------------------
function mnu_File_Callback(hObject, eventdata, handles)
% hObject    handle to mnu_File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Load a configuration file and update Property Grids ------------
function mnu_LoadCfg_Callback(hObject, eventdata, handles)
% hObject    handle to mnu_LoadCfg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% load the configs
[fname fpath] = uigetfile('*.mat','Load Config File');
if ~fname
    return;
end
tmp = load(fullfile(fpath,fname));
if isfield(tmp,'opts')
    handles.opts = tmp.opts;
else
    handles.opts = tmp;
end

% redraw the property grids
handles = redrawPropertyGrids(hObject,handles);
guidata(hObject,handles);

% --------------------------------------------------------------------
function mnu_SaveCfg_Callback(hObject, eventdata, handles)
% hObject    handle to mnu_SaveCfg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cmdSaveConfig_Callback(hObject,eventdata,handles);


% --- Executes on key press with focus on gui_MetaControlPanel or any of its controls.
function gui_MetaControlPanel_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to gui_MetaControlPanel (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
