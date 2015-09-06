function varargout = gui_causalBrainMovie3D_online(varargin)
%
% GUI_CAUSALBRAINMOVIE3D_ONLINE M-file for gui_causalBrainMovie3D_online.fig
%      GUI_CAUSALBRAINMOVIE3D_ONLINE, by itself, creates a new GUI_CAUSALBRAINMOVIE3D_ONLINE or raises the existing
%      singleton*.
%
%      H = GUI_CAUSALBRAINMOVIE3D_ONLINE returns the handle to a new GUI_CAUSALBRAINMOVIE3D_ONLINE or the handle to
%      the existing singleton*.
%
%      GUI_CAUSALBRAINMOVIE3D_ONLINE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_CAUSALBRAINMOVIE3D_ONLINE.M with the given input arguments.
%
%      GUI_CAUSALBRAINMOVIE3D_ONLINE('Property','Value',...) creates a new GUI_CAUSALBRAINMOVIE3D_ONLINE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_causalBrainMovie3D_online_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_causalBrainMovie3D_online_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_causalBrainMovie3D_online

% Last Modified by GUIDE v2.5 19-Feb-2013 00:02:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_causalBrainMovie3D_online_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_causalBrainMovie3D_online_OutputFcn, ...
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


% --- Executes just before gui_causalBrainMovie3D_online is made visible.
function gui_causalBrainMovie3D_online_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_causalBrainMovie3D_online (see VARARGIN)

% Choose default command line output for gui_causalBrainMovie3D_online
handles.output = hObject;
handles.gcf    = hObject;
set(hObject,'Name','BrainMovie3D Control Panel');
% set default termination behavior
handles.ExitButtonClicked = 'Cancel';

% extract some data from command-line input
if isempty(varargin)
    error('You must pass ALLEEG and Conn to gui_CausalBrainMovie3D');
end

% Extract input parameters/data and store
handles.ud.ALLEEG  = varargin{1};  
handles.ud.Conn    = varargin{2};  
varargin([1 2]) = [];

% set default EEGLAB background and text colors
%-----------------------------------------------
try, icadefs;
catch,
	GUIBACKCOLOR        =  [.8 .8 .8];     
	GUIPOPBUTTONCOLOR   = [.8 .8 .8];    
	GUITEXTCOLOR        = [0 0 0];
end;

allhandlers = hObject;

hh = findobj(allhandlers,'style', 'text');
%set(hh, 'BackgroundColor', get(hObject, 'color'), 'horizontalalignment', 'left');
set(hh, 'Backgroundcolor', GUIBACKCOLOR);
set(hh, 'foregroundcolor', GUITEXTCOLOR);
set(hObject, 'color',GUIBACKCOLOR );
% set(hh, 'horizontalalignment', g.horizontalalignment);

hh = findobj(allhandlers, 'style', 'edit');
set(hh, 'BackgroundColor', [1 1 1]); %, 'horizontalalignment', 'right');

hh =findobj(allhandlers, 'parent', hObject, 'style', 'pushbutton');
if ~strcmpi(computer, 'MAC') && ~strcmpi(computer, 'MACI') % this puts the wrong background on macs
    set(hh, 'backgroundcolor', GUIPOPBUTTONCOLOR);
    set(hh, 'foregroundcolor', GUITEXTCOLOR);
end;
hh =findobj(allhandlers, 'parent', hObject, 'style', 'popupmenu');
set(hh, 'backgroundcolor', GUIPOPBUTTONCOLOR);
set(hh, 'foregroundcolor', GUITEXTCOLOR);
hh =findobj(allhandlers, 'parent', hObject, 'style', 'checkbox');
set(hh, 'backgroundcolor', GUIBACKCOLOR);
set(hh, 'foregroundcolor', GUITEXTCOLOR);
hh =findobj(allhandlers, 'parent', hObject, 'style', 'listbox');
set(hh, 'backgroundcolor', GUIPOPBUTTONCOLOR);
set(hh, 'foregroundcolor', GUITEXTCOLOR);
hh =findobj(allhandlers, 'parent', hObject, 'style', 'radio');
set(hh, 'foregroundcolor', GUITEXTCOLOR);
set(hh, 'backgroundcolor', GUIPOPBUTTONCOLOR);
set(hObject, 'visible', 'on');

set(handles.pnlPropertyGrid,'backgroundcolor', GUIBACKCOLOR);
set(handles.pnlPropertyGrid,'foregroundcolor', GUITEXTCOLOR);

%-----------------------------------------------

drawnow

% render the PropertyGrid in the correct panel
handles.PropertyGridHandle = arg_guipanel( ...
                 handles.pnlPropertyGrid, ...
                'Function',@vis_causalBrainMovie3D, ...
                'params',{'ALLEEG',handles.ud.ALLEEG, 'Conn',handles.ud.Conn, varargin{:}});

handles.BrainMovieFigure = [];
handles.figureHandle = hObject;

% Update handles structure
guidata(hObject, handles);

% Wait for user to click OK, Cancel or close figure
% uiwait(handles.gui_BrainMovie3D);


% UIWAIT makes gui_causalBrainMovie3D_online wait for user response (see UIRESUME)
% uiwait(handles.gui_BrainMovie3D);





% --- Outputs from this function are returned to the command line.
function varargout = gui_causalBrainMovie3D_online_OutputFcn(hObject, eventdata, handles) 
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

evalin('base','opts.miscOptCfg.doSIFT.dispBrainMovie.arg_selection=false;');
% evalin('base','if ishandle(figHandles.BMDisplay), close(figHandles.BMDisplay); end');

try
    close(handles.figureHandle);
catch
end



% --- Executes on button press in cmdOK.
function cmdOK_Callback(hObject, eventdata, handles)

% handles.ExitButtonClicked ='OK';
% get the property specification
pg = handles.PropertyGridHandle;
ps = pg.GetPropertySpecification;
cfg = arg_tovals(ps,false);
assignin('base','BMCFG',cfg);
assignin('base','BM_CFG_CHANGED',true);
set(handles.gcf,'UserData','newBrainmovie');
% guidata(hObject,handles);
% uiresume(handles.gui_BrainMovie3D);




% % --- Executes on slider movement.
% function slideCurTime_Callback(hObject, eventdata, handles)
% 
% % get the PropertyGrid object
% % ud = get(handles.gui_BrainMovie3D,'UserData');
% % if isempty(ud)
% %     error('Could not obtain handle to PropertyGrid! Did you pass it to the UserData of the figure?');
% % end
% 
% % get the property specification
% pg = handles.PropertyGridHandle;
% ps = pg.GetPropertySpecification;
% cfg = arg_tovals(ps,false);
% 
% if isempty(cfg)
%     error('could not retrieve property specification from propertygrid');
% end
% 
% % get the index of the frame the user wants to render
% range = get(hObject,'max')-get(hObject,'min');
% curframe = round(get(hObject,'Value'));
% 
% % update the curTime textbox
% set(handles.txtCurTime,'string',handles.ud.Conn(1).erWinCenterTimes(curframe));
% 
% minTime = str2num(get(handles.txtMinTime,'String'));
% maxTime = str2num(get(handles.txtMaxTime,'String'));
% 
% figh = handles.BrainMovieFigure; %findobj('tag','BrainMovieFigure');
% try figure(figh)
% catch
%     figh = [];
% end
% curpos = get(findobj(figh,'tag','brain'),'CameraPosition');
% % execute the brainmovie function to render this frame
% vis_causalBrainMovie3D('ALLEEG',handles.ud.ALLEEG,'Conn',handles.ud.Conn,cfg,'MovieTimeRange',[minTime maxTime],'BrainMovieOptions',{cfg.BMopts,'visible','off','FigureHandle',figh,'LatenciesToRender',[],'FramesToRender',curframe,'outputFormat',{'framefolder','','moviename',''}});
% figh = gcf;
% set(findobj(figh,'tag','brain'),'CameraPosition',curpos);
% set(gcf,'visible','on')
% handles.BrainMovieFigure = gcf;
% guidata(hObject,handles);



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
function cmdOK_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to cmdOK (see GCBO)
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


% --- Executes on key press with focus on gui_BrainMovie3D or any of its controls.
function gui_BrainMovie3D_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to gui_BrainMovie3D (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

% hObject    handle to gui_BrainMovie3D (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

% if strcmpi(eventdata.Key,'s') ...
%         && (strcmpi(eventdata.Modifier,'control') ...
%         || strcmpi(eventdata.Modifier,'command'))
%     
%     
%     % get current contents of property grid
%     pg    = handles.PropertyGridHandle;
%     ps    = pg.GetPropertySpecification;
%     BMCFG = arg_tovals(ps,false);
% 
%     % select the path/file for saving configs
%     [fname fpath] = uiputfile('*.mat','Save BrainMovie Config File');
%     if ~fname
%         return;
%     end
%     % save the opts structure
%     save(fullfile(fpath,fname),'BMCFG');
% 
%     
% %     varname = inputdlg2({'Enter name of variable in workspace to store configuration'}, ...
% %                         'Save GUI configuration to workspace',1,{''});
% %      
% %     if isempty(varname)
% %         return;
% %     end
% %     
% %     % get the property specification
% %     pg = handles.PropertyGridHandle;
% %     ps = pg.GetPropertySpecification;
% %     cfg = arg_tovals(ps,false);
% %     
% %     % save configuration structure to workspace
% %     assignin('base',varname{1},cfg);
% end
% 
% % load configuration
% if strcmpi(eventdata.Key,'l') ...
%         && (strcmpi(eventdata.Modifier,'control') ...
%         || strcmpi(eventdata.Modifier,'command'))
%     
%     % load the configs
%     [fname fpath] = uigetfile('*.mat','Load BrainMovie Config File');
%     if ~fname
%         return;
%     end
%     tmp   = load(fullfile(fpath,fname));
%     fn    = fieldnames(tmp);
%     BMCFG = tmp.(fn{1});
% 
%     % pause the brainmovie rendering
%     evalin('base','opts.miscOptCfg.doSIFT.dispBrainMovie.arg_selection=false;');
%     
%     % close existing brainmovie figure 
%     evalin('base','close(figHandles.BMDisplay);  figHandles.BMDisplay = [];');
%     % update configs in base workspace
%     assignin('base','BMCFG',BMCFG);
%     assignin('base','BM_CFG_CHANGED',true);
%     assignin('base','newBrainmovie',true);
%     % close the BrainMovie control panel (a new one will be opened)
%     close(hObject);
%     % resume the brainmovie rendering
%     evalin('base','opts.miscOptCfg.doSIFT.dispBrainMovie.arg_selection=true;');
%     
% %     % pause the brainmovie rendering
% %     evalin('base','opts.miscOptCfg.doSIFT.dispBrainMovie.arg_selection=false;');
% %     % redraw the property grids
% %     handles = redrawPropertyGrid(hObject,handles,BMCFG);
% %     % close existing brainmovie figure 
% %     evalin('base','close(figHandles.BMDisplay);  figHandles.BMDisplay = [];');
% %     % update configs in base workspace
% %     assignin('base','BMCFG',BMCFG);
% %     assignin('base','BM_CFG_CHANGED',true);
% %     % resume the brainmovie rendering
% %     evalin('base','opts.miscOptCfg.doSIFT.dispBrainMovie.arg_selection=true;');
% 
%     guidata(hObject,handles);
% end


% --- User-defined function to (re)draw Property Grids
function handles = redrawPropertyGrid(hObject,handles,BMCFG)

try
    % delete existing property grids
    delete(handles.PropertyGridHandle.Control);
catch err
end
 
% render the PropertyGrid in the correct panel
handles.PropertyGridHandle = arg_guipanel( ...
                 handles.pnlPropertyGrid, ...
                'Function',@vis_causalBrainMovie3D, ...
                'params',{'ALLEEG',handles.ud.ALLEEG, 'Conn',handles.ud.Conn, BMCFG});


             
             
             
% --- Executes on button press in cmdPause.
function cmdPause_Callback(hObject, eventdata, handles)
% hObject    handle to cmdPause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% pause the brainmovie
% pause the brainmovie
if strcmpi(get(hObject,'String'),'pause');
    doBrainMovie = false;
    set(hObject,'String','Unpause');
else
    doBrainMovie = true;
    set(hObject,'String','Pause');
end

evalin('base',sprintf('opts.miscOptCfg.doSIFT.dispBrainMovie.arg_selection=%s;',fastif(doBrainMovie,'true','false')));
guidata(hObject,handles);


% --------------------------------------------------------------------
function mnuLoad_Callback(hObject, eventdata, handles)
% hObject    handle to mnuLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% load the configs
[fname fpath] = uigetfile('*.mat','Load BrainMovie Config File');
if ~fname
    return;
end
tmp   = load(fullfile(fpath,fname));
fn    = fieldnames(tmp);
BMCFG = tmp.(fn{1});

if ~isstruct(BMCFG) || ~isfield(BMCFG,'BMopts')
    fprintf('Invalid BrainMovie Config\n');
    return; 
end

% pause the brainmovie rendering
evalin('base','opts.miscOptCfg.doSIFT.dispBrainMovie.arg_selection=false;');

% close existing brainmovie figure 
evalin('base','close(figHandles.BMDisplay);  figHandles.BMDisplay = [];');
% update configs in base workspace
assignin('base','BMCFG',BMCFG);
assignin('base','BM_CFG_CHANGED',true);
assignin('base','newBrainmovie',true);
% close the BrainMovie control panel (a new one will be opened)
close(handles.figureHandle);
% resume the brainmovie rendering
evalin('base','opts.miscOptCfg.doSIFT.dispBrainMovie.arg_selection=true;');

%     % pause the brainmovie rendering
%     evalin('base','opts.miscOptCfg.doSIFT.dispBrainMovie.arg_selection=false;');
%     % redraw the property grids
%     handles = redrawPropertyGrid(hObject,handles,BMCFG);
%     % close existing brainmovie figure 
%     evalin('base','close(figHandles.BMDisplay);  figHandles.BMDisplay = [];');
%     % update configs in base workspace
%     assignin('base','BMCFG',BMCFG);
%     assignin('base','BM_CFG_CHANGED',true);
%     % resume the brainmovie rendering
%     evalin('base','opts.miscOptCfg.doSIFT.dispBrainMovie.arg_selection=true;');

guidata(handles.figureHandle,handles);
    

% --------------------------------------------------------------------
function mnuSave_Callback(hObject, eventdata, handles)
% hObject    handle to mnuSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get current contents of property grid
pg    = handles.PropertyGridHandle;
ps    = pg.GetPropertySpecification;
BMCFG = arg_tovals(ps,false);

% select the path/file for saving configs
[fname fpath] = uiputfile('*.mat','Save BrainMovie Config File');
if ~fname
    return;
end
% save the opts structure
save(fullfile(fpath,fname),'BMCFG');


%     varname = inputdlg2({'Enter name of variable in workspace to store configuration'}, ...
%                         'Save GUI configuration to workspace',1,{''});
%      
%     if isempty(varname)
%         return;
%     end
%     
%     % get the property specification
%     pg = handles.PropertyGridHandle;
%     ps = pg.GetPropertySpecification;
%     cfg = arg_tovals(ps,false);
%     
%     % save configuration structure to workspace
%     assignin('base',varname{1},cfg);


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
