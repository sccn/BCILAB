function env_showmenu(varargin)
% Links the BCILAB menu into another menu, or creates a new root menu if necessary.
% env_showmenu(Options...)
%
% In:
%   Options... : optional name-value pairs; names are:
%                 'parent': parent menu to link into
%
%                 'shortcuts': whether to enable keyboard shortcuts
%
%                 'forcenew': whether to force creation of a new menu
%
% Example:
%   % bring up the BCILAB main menu if it had been closed
%   env_showmenu;
%
% See also:
%   env_startup
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-10-29

% parse options...
hlp_varargin2struct(varargin,'parent',[], 'shortcuts',true,'forcenew',false);

% check if we're an EEGLAB plugin
folders = hlp_split(fileparts(mfilename('fullpath')),filesep);
within_eeglab = length(folders) >= 5 && strcmp(folders{end-3},'plugins') && ~isempty(strfind(folders{end-4},'eeglab'));

% highlight the menu if already there
f = findobj('Tag','bcilab_menu');
if ~isempty(f) && ~forcenew    
    close(get(f,'Parent')); end

if isempty(parent) %#ok<NODEF>
    if within_eeglab && ~forcenew
        % try to link into the EEGLAB main menu
        try
            toolsmenu = findobj(0,'tag','tools');
            if ~isempty(toolsmenu)
                parent = uimenu(toolsmenu, 'Label','BCILAB');
                set(toolsmenu,'Enable','on');
            end
        catch
            disp('Unable to link BCILAB menu into EEGLAB menu.');
        end
    end
    if isempty(parent)
        % create new root menu, if no parent
        from_left = 100;
        from_top = 150;
        width = 500;
        height = 1;
        % determine position on primary monitor
        import java.awt.GraphicsEnvironment
        ge = GraphicsEnvironment.getLocalGraphicsEnvironment();
        gd = ge.getDefaultScreenDevice();
        scrheight = gd.getDisplayMode().getHeight();
        pos = [from_left, scrheight-from_top, width, height];        
        % create figure
        figtitle = ['BCILAB ' env_version ' (on ' hlp_hostname ')'];
        parent = figure('DockControls','off','NumberTitle','off','Name',figtitle,'Resize','off','MenuBar','none','Position',pos,'Tag','bcilab_toolwnd');
    end
end

% Data Source menu
source = uimenu(parent, 'Label','Data Source','Tag','bcilab_menu');
uimenu(source,'Label','Load recording(s)...','Accelerator',char(shortcuts*'l'),'Callback','gui_loadset');
wspace = uimenu(source,'Label','Workspace','Separator','on');
uimenu(wspace,'Label','Load...','Callback','io_loadworkspace');
uimenu(wspace,'Label','Save...','Callback','io_saveworkspace');
uimenu(wspace,'Label','Clear...','Callback','clear');
uimenu(source,'Label','Run script...','Separator','on','Callback',@invoke_script);
if isdeployed
    uimenu(source,'Label','Quit','Separator','on','Callback','exit'); end

% Offline Analysis menu
offline = uimenu(parent, 'Label','Offline Analysis');
uimenu(offline,'Label','New approach...','Accelerator',char(shortcuts*'n'),'Callback','gui_newapproach');
uimenu(offline,'Label','Modify approach...','Accelerator',char(shortcuts*'m'),'Callback','gui_configapproach([],true);');
uimenu(offline,'Label','Review/edit approach...','Accelerator',char(shortcuts*'r'),'Callback','gui_reviewapproach([],true);');
uimenu(offline,'Label','Save approach...','Accelerator',char(shortcuts*'s'),'Callback','gui_saveapproach');
uimenu(offline,'Label','Train new model...','Accelerator',char(shortcuts*'t'),'Callback','gui_calibratemodel','Separator','on');
uimenu(offline,'Label','Apply model to data...','Accelerator',char(shortcuts*'a'),'Callback','gui_applymodel');
uimenu(offline,'Label','Visualize model...','Accelerator',char(shortcuts*'v'),'Callback','gui_visualizemodel');
uimenu(offline,'Label','Run batch analysis...','Accelerator',char(shortcuts*'b'),'Callback','gui_batchanalysis','Separator','on');
uimenu(offline,'Label','Review results...','Accelerator',char(shortcuts*'i'),'Callback','gui_selectresults','Separator','on');

% Online Analysis menu
online = uimenu(parent,'Label','Online Analysis');
pipe = uimenu(online,'Label','Process data within...');
read = uimenu(online,'Label','Read input from...');
write = uimenu(online,'Label','Write output to...');
cm_read = uicontextmenu('Tag','bcilab_cm_read');
cm_write = uicontextmenu('Tag','bcilab_cm_write');

% for each plugin sub-directory...
dirs = dir(env_translatepath('functions:/online_plugins'));
for d={dirs(3:end).name}
    % find all files, their names, identifiers, and function handles
    files = dir(env_translatepath(['functions:/online_plugins/' d{1} '/run_*.m']));
    names = {files.name};
    idents = cellfun(@(n)n(1:end-2),names,'UniformOutput',false);
    % for each entry...
    for f=1:length(idents)
        try
            if ~exist(idents{f},'file') && ~isdeployed
                addpath(env_translatepath(['functions:/online_plugins/' d{1}])); end
            % get properties...
            props = arg_report('properties',str2func(idents{f}));
            % get category
            if strncmp(idents{f},'run_read',8);
                cats = [read,cm_read];
            elseif strncmp(idents{f},'run_write',9);
                cats = [write,cm_write];
            elseif strncmp(idents{f},'run_pipe',8);
                cats = pipe;
            end
            if isfield(props,'name')
                % add menu entry
                for cat=cats
                    uimenu(cat,'Label',[props.name '...'],'Callback',['arg_guidialog(@' idents{f} ');'],'Enable','on'); end
            else
                warning('env_showmenu:missing_guiname','The online plugin %s does not declare a GUI name; ignoring...',idents{f});
            end
        catch
            disp(['Could not integrate the online plugin ' idents{f} '.']);
        end
    end
end
uimenu(online,'Label','Clear all online processing','Callback','onl_clear','Separator','on');

% Settings menu
settings = uimenu(parent, 'Label','Settings');
uimenu(settings,'Label','Directory settings...','Callback','gui_configpaths');
uimenu(settings,'Label','Cache settings...','Callback','gui_configcache');
uimenu(settings,'Label','Cluster settings...','Callback','gui_configcluster');
uimenu(settings,'Label','Clear memory cache','Callback','env_clear_memcaches','Separator','on');

% Help menu
helping = uimenu(parent,'Label','Help');
uimenu(helping,'Label','BCI Paradigms...','Callback','env_doc code/paradigms');
uimenu(helping,'Label','Filters...','Callback','env_doc code/filters');
uimenu(helping,'Label','Machine Learning...','Callback','env_doc code/machine_learning');
scripting = uimenu(helping,'Label','Scripting');
uimenu(scripting,'Label','File input/output...','Callback','env_doc code/io');
uimenu(scripting,'Label','Dataset editing...','Callback','env_doc code/dataset_editing');
uimenu(scripting,'Label','Offline scripting...','Callback','env_doc code/offline_analysis');
uimenu(scripting,'Label','Online scripting...','Callback','env_doc code/online_analysis');
uimenu(scripting,'Label','BCILAB environment...','Callback','env_doc code/environment');
uimenu(scripting,'Label','Cluster handling...','Callback','env_doc code/parallel');
uimenu(scripting,'Label','Keywords...','Separator','on','Callback','env_doc code/keywords');
uimenu(scripting,'Label','Helpers...','Callback','env_doc code/helpers');
uimenu(scripting,'Label','Internals...','Callback','env_doc code/utils');
authoring = uimenu(helping,'Label','Plugin authoring');
uimenu(authoring,'Label','Argument declaration...','Callback','env_doc code/arguments');
uimenu(authoring,'Label','Expression functions...','Callback','env_doc code/expressions');
uimenu(authoring,'Label','Online processing...','Callback','env_doc code/online_analysis');
uimenu(helping,'Label','About...','Separator','on','Callback',@about);
uimenu(helping,'Label','Save bug report...','Separator','on','Callback','io_saveworkspace([],true)');
uimenu(helping,'Label','File bug report...','Callback','arg_guidialog(@env_bugreport);');

% toolbar (if not linked into the EEGLAB menu)
if ~(within_eeglab && ~forcenew)
    global tracking;
    cluster_requested = isfield(tracking,'cluster_requested') && ~isempty(tracking.cluster_requested);
    cluster_requested = hlp_rewrite(cluster_requested,false,'off',true,'on');
    
    ht = uitoolbar(parent,'HandleVisibility','callback');
    uipushtool(ht,'TooltipString','Load recording(s)',...
        'CData',load_icon('bcilab:/resources/icons/file_open.png'),...
        'HandleVisibility','callback','ClickedCallback','gui_loadset');
    uipushtool(ht,'TooltipString','New approach',...
        'CData',load_icon('bcilab:/resources/icons/approach_new.png'),...
        'HandleVisibility','callback','ClickedCallback','gui_newapproach','Separator','on');
    uipushtool(ht,'TooltipString','Load Approach',...
        'CData',load_icon('bcilab:/resources/icons/approach_load.png'),...
        'HandleVisibility','callback','ClickedCallback',@load_approach);
    uipushtool(ht,'TooltipString','Save approach',...
        'CData',load_icon('bcilab:/resources/icons/approach_save.png'),...
        'HandleVisibility','callback','ClickedCallback','gui_saveapproach');
    uipushtool(ht,'TooltipString','Modify approach',...
        'CData',load_icon('bcilab:/resources/icons/approach_edit.png'),...
        'HandleVisibility','callback','ClickedCallback','gui_configapproach([],true);');
    uipushtool(ht,'TooltipString','Review/edit approach',...
        'CData',load_icon('bcilab:/resources/icons/approach_review.png'),...
        'HandleVisibility','callback','ClickedCallback','gui_reviewapproach([],true);');
    uipushtool(ht,'TooltipString','Train new model',...
        'CData',load_icon('bcilab:/resources/icons/model_new.png'),...
        'HandleVisibility','callback','ClickedCallback','gui_calibratemodel','Separator','on');
    uipushtool(ht,'TooltipString','Load Model',...
        'CData',load_icon('bcilab:/resources/icons/model_load.png'),...
        'HandleVisibility','callback','ClickedCallback',@load_model);
    uipushtool(ht,'TooltipString','Save Model',...
        'CData',load_icon('bcilab:/resources/icons/model_save.png'),...
        'HandleVisibility','callback','ClickedCallback','gui_savemodel');
    uipushtool(ht,'TooltipString','Apply model to data',...
        'CData',load_icon('bcilab:/resources/icons/model_apply.png'),...
        'HandleVisibility','callback','ClickedCallback','gui_applymodel');
    uipushtool(ht,'TooltipString','Visualize model',...
        'CData',load_icon('bcilab:/resources/icons/model_visualize.png'),...
        'HandleVisibility','callback','ClickedCallback','gui_visualizemodel');
    uipushtool(ht,'TooltipString','Run batch analysis',...
        'CData',load_icon('bcilab:/resources/icons/batch_analysis.png'),...
        'HandleVisibility','callback','ClickedCallback','gui_batchanalysis');
    uipushtool(ht,'TooltipString','Read input from (online)',...
        'CData',load_icon('bcilab:/resources/icons/online_in.png'),...
        'HandleVisibility','callback','Separator','on','ClickedCallback',@click_read);
    uipushtool(ht,'TooltipString','Write output to (online)',...
        'CData',load_icon('bcilab:/resources/icons/online_out.png'),...
        'HandleVisibility','callback','ClickedCallback',@click_write);
    uipushtool(ht,'TooltipString','Clear online processing',...
        'CData',load_icon('bcilab:/resources/icons/online_clear.png'),...
        'HandleVisibility','callback','ClickedCallback','onl_clear');
    uitoggletool(ht,'TooltipString','Request cluster availability',...
        'CData',load_icon('bcilab:/resources/icons/acquire_cluster.png'),'HandleVisibility','callback','Separator','on','State',cluster_requested,'OnCallback','env_acquire_cluster','OffCallback','env_release_cluster');
    uipushtool(ht,'TooltipString','About BCILAB',...
        'CData',load_icon('bcilab:/resources/icons/help.png'),'HandleVisibility','callback','Separator','on','ClickedCallback',@about);
end

if within_eeglab && forcenew
    mainmenu = findobj('Tag','EEGLAB');
    % make the EEGLAB menu current again
    if ~isempty(mainmenu)
        figure(mainmenu); end
end


function about(varargin)
infotext = strvcat(...
    'BCILAB is an open-source toolbox for Brain-Computer Interfacing research.', ...
    'It is being developed by Christian Kothe at the Swartz Center for Computational Neuroscience,',...
    'Institute for Neural Computation (University of California San Diego).', ...
    ' ',...
    'Development of this software was supported by the Army Research Laboratories under', ...
    'Cooperative Agreement Number W911NF-10-2-0022, as well as by a gift from the Swartz Foundation.', ...
    ' ',...
    'The design was inspired by the preceding PhyPA toolbox, written by C. Kothe and T. Zander', ...
    'at the Berlin Institute of Technology, Chair Human-Machine Systems.', ...
    ' ',...
    'BCILAB connects to the following toolboxes/libraries:', ...
    '* AWS SDK (Amazon)', ...
    '* Amica (SCCN/UCSD)', ...
    '* Chronux (Mitra Lab, Cold Spring Harbor)', ...
    '* CVX (Stanford)', ...
    '* DAL (U. Tokyo)', ...
    '* DataSuite (SCCN/UCSD)', ...
    '* EEGLAB (SCCN/UCSD)', ...
    '* BCI2000import (www.bci2000.org)', ...
    '* Logreg (Jan Drugowitsch)', ...
    '* FastICA (Helsinki UT)', ...
    '* glm-ie (Max Planck Institute for Biological Cybernetics, Tuebingen)', ...
    '* glmnet (Stanford)', ...
    '* GMMBayes (Helsinki UT)', ...
    '* HKL (Francis Bach, INRIA)', ...
    '* KernelICA (Francis Bach, Berkeley)', ...
    '* LIBLINEAR (National Taiwan University)', ...
    '* matlabcontrol (Joshua Kaplan)', ...
    '* mlUnit (Thomas Dohmke)', ...
    '* NESTA (Caltech)', ...
    '* OSC (Andy Schmeder) and LibLO (Steve Harris)', ...
    '* PROPACK (Stanford)', ...
    '* PropertyGrid (Levente Hunyadi)', ...
    '* SparseBayes (Vector Anomaly)', ...
    '* SVMlight (Thorsten Joachims)', ...
    '* SVMperf (Thorsten Joachims)', ...
    '* Talairach (UTHSCSA)', ...
    '* Time-frequency toolbox (CRNS / Rice University)', ...
    '* t-SNE (TU Delft)', ...    
    '* UnLocBox (Nathanael Perraudin et al.)', ...
    '* VDPGM (Kenichi Kurihara)'); %#ok<REMFF1>

warndlg2(infotext,'About');


function click_read(varargin)
% pop up a menu when clicking the "read input from" toolbar button
tw = findobj('tag','bcilab_toolwnd');
cm = findobj('tag','bcilab_cm_read');
tpos = get(tw,'Position');
ppos = get(0,'PointerLocation');
set(cm,'Position',ppos - tpos(1:2), 'Visible', 'on');
set(cm,'Position',ppos - tpos(1:2), 'Visible', 'on');


function click_write(varargin)
% pop up a menu when clicking the "write output to" toolbar button
tw = findobj('tag','bcilab_toolwnd');
cm = findobj('tag','bcilab_cm_write');
tpos = get(tw,'Position');
ppos = get(0,'PointerLocation');
set(cm,'Position',ppos - tpos(1:2), 'Visible', 'on');
set(cm,'Position',ppos - tpos(1:2), 'Visible', 'on');


% run a script
function invoke_script(varargin)
[filename,filepath] = uigetfile({'*.m', 'MATLAB file'},'Select script to run',env_translatepath('bcilab:/'));
if ~isnumeric(filename)
    run_script([filepath filename],true); end


% load a toolbar icon
function cols = load_icon(filename)
[cols,palette,alpha] = imread(env_translatepath(filename));
if ~isempty(palette)
    error('This function does not handle palettized icons.'); end
cls = class(cols);
cols = double(cols);
cols = cols/double(intmax(cls));
cols([alpha,alpha,alpha]==0) = NaN;


% load a BCI model from disk
function load_model(varargin)
[filename,filepath] = uigetfile({'*.mdl', 'BCI Model'},'Select BCI model to load',env_translatepath('home:/.bcilab/models'));
if ~isnumeric(filename)
    contents = io_load([filepath filename],'-mat');
    for fld=fieldnames(contents)'
        tmp = contents.(fld{1});
        if isstruct(tmp) && isfield(tmp,'timestamp')
            tmp.timestamp = now;
            assignin('base',fld{1},tmp);
        end
    end
end

% load a BCI approach from disk
function load_approach(varargin)
[filename,filepath] = uigetfile({'*.apr', 'BCI Approach'},'Select BCI approach to load',env_translatepath('home:/.bcilab/approaches'));
if ~isnumeric(filename)
    contents = io_load([filepath filename],'-mat');
    for fld=fieldnames(contents)'
        tmp = contents.(fld{1});
        if isstruct(tmp) && isfield(tmp,'paradigm')
            tmp.timestamp = now;
            assignin('base',fld{1},tmp);
        end
    end
end
