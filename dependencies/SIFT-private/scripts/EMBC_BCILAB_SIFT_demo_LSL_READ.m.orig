%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REAL-TIME BCILAB/SIFT INTERFACE DEMO SCRIPT
%
% This script demonstrates a beta interface between BCILAB 
% (Christian Kothe: www.sccn.ucsd.edu/wiki/BCILAB) and SIFT (Tim Mullen:  
% www.sccn.ucsd.edu/wiki/SIFT) to perform real-time data extraction,
% cleaning, connectivity modeling and visualization.
% 
% This script was prepared for the following whitepaper/demonstration:
% [1] Tim Mullen, Christian Kothe, Yu Chi, Tzyy-Ping Jung,
%     and Scott Makeig. "Real-Time Estimation and 3D Visualization of Sparse 
%     Multivariate Effective Connectivity Using Wearable EEG." July 15th,
%     IEEE EMB/CAS/SMC Workshop on Brain-Machine-Body Interfaces. 34th 
%     Annual International IEEE EMBS Conference of the IEEE Engineering in 
%     Medicine and Biology Society (IEEE EMBS).
%
% Please refer to [1] when using any functionality from this demo.
% 
% Copyright Tim Mullen, Christian Kothe, July 2012, SCCN/INC/UCSD
%
% Author: Tim Mullen,       July 2012, SCCN/INC/UCSD
%         Christian Kothe,  July 2012, SCCN/INC/UCSD
%
% This function is part of the Source Information Flow Toolbox (SIFT)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Tip: to create an evaluable list of good channels (for flt_selchans) try this:
% hlp_tostring(setdiff({signal.chanlocs.labels},badchannels))

%% SET UP CONFIGURATION OPTIONS
MAX_AGE         = 1;
MAX_BUFFERED    = 5;
MIN_MEM_LIMIT   = 500; % min memory threshold (MB)
TRAIN_ONLY      = false;
ROTATE90        = false;

% Source reconstruction options (leave disabled)
COLOR_SOURCE_ROI = true;          % this will use special meshes for coloring ROIs
HEAD_MODEL_NAME  = 'resources:/headmodels/standard-Colin27-385ch.mat'; %'data:/mobilab/Cognionics_64_HeadModelObj_3751.mat'; %'data:/mobilab/Cognionics_64_HeadModelObj_11997.mat';  %'resources:/headmodels/standard-Colin27-385ch.mat';             % path to head model object for source reconstruction (relative to 'datapath'). Leave empty if you aren't doing source reconstruction
% Establish file paths
% NOTE: all paths are relative to 'datapath' which is a
% platform-independent path which itself can be relative to bcilab root
% folder (i.e. data:/ is the userdata folder in the bcilab root dir)
datapath         = 'data:/';       % this is relative to the BCILAB root dir
TrainingDataFile = 'Cognionics_64_Flanker_85_265.set'; %'calibration.xdf'; %'Cognionics_64_training.set'; %'Cognionics_64_Flanker.set'; %'Cognionics_64_Flanker_85_265.set'; 'Cognionics_64_SIMULATION_one_source.set'; %'Cognionics_64_Flanker_85_265.set';  %'Cognionics_64_Flanker.set';  %'Cognionics_64_Flanker_0_10.set'; %'Cognionics_64_Flanker.set'; %'Cognionics_64_training.set'; %'Cognionics_64_Flanker_85_265.set'; %'Cognionics_64_training.set'; %'calibration_mindo.xdf'; % %'calibration.xdf'; %'Cognionics_Pyramind_demo.set'; %'clean_reversed.xdf'; %'noisy.xdf'; %'Cognionics_Pyramind_demo.set';             % this is the relative path to the calibration dataset
TestingDataFile  = 'Cognionics_64_Flanker.set';  %'calibration_old1.xdf'; %'Cognionics_64_testing.set'; %'Cognionics_64_SIMULATION_manysources_nocsdsaved.set'; %'Cognionics_64_Flanker.set'; %'Cognionics_64_SIMULATION.set'; %'Cognionics_64_Flanker.set'; %'Cognionics_64_Flanker_85_265.set'; %'Cognionics_64_Flanker.set'; %'Cognionics_64_testing.set'; %'testing.xdf'; %'Cognionics_Pyramind_demo.set'; %'clean_reversed.xdf'; %'noisy.xdf'; %'Cognionics_Pyramind_demo.set';             % this is an optional path to a dataset to playback (if RUN_LSL = false)
GUI_CONFIG_NAME  = 'BMCFG_RECORD_STABILITY_TEST_VBLORETA.mat'; %'Cognionics_64_Pipeline_Demo_METACP_CFG.mat'; %'EMBC_PAPER_METACP_OPTS_NOSOURCES.mat'; %'DEMO_SOURCELOC_METACP_CFG_CombineROIs_nodelay_manyROIs_autochansel.mat'; %'SIMULATION_TEST_LORETA.mat'; %'DEMO_SOURCELOC_METACP_CFG_AllVertices.mat'; %'DEMO_SOURCELOC_METACP_CFG_CombineROIs.mat'; %'DARPA_DEMO_METACP_CFG_FEWCHANS.mat';             % relative path to a default pipeline configuration
GUI_BRAINMOVIE_CONFIG_NAME = 'DARPA_DEMO_BM_CFG.mat'; %'DEMO_SOURCELOC_BM_CFG.mat'; %'DARPA_DEMO_BM_CFG.mat';   % relative path to BrainMovie configuration

% Set up the name of the stream we will write to in the workspace
streamname              = 'EEGDATA';
LSL_SelectionProperty   = 'type';  % can also be 'name'
LSL_SelectionValue      = 'EEG';

% this is the name of the stream that carries SIFT models
INSTREAM_NAME = 'SIFT-Models';

%% initialize the input stream
stream_inlet = onl_lslreceive_init(INSTREAM_NAME,Inf,MAX_BUFFERED);

if ~isempty(MIN_MEM_LIMIT)
    % set up a timer for checking available memory
    tmr_checkmem = timer('Period', 2,'ExecutionMode','fixedRate', ...
        'TimerFcn',{@onl_lsl_bufmemcheck,stream_inlet,MIN_MEM_LIMIT},'StartDelay',0,'Tag','checkmem_timer');
end

%% load head model object
if ~isempty(HEAD_MODEL_NAME)
    hmObj = headModel.loadFromFile(env_translatepath(HEAD_MODEL_NAME));
    % load meshes for visualization
    surfData    = load(hmObj.surfaces);
    surfData    = surfData.surfData;
    scalpMesh   = surfData(1);
    csfMesh     = surfData(2);
    cortexMesh  = surfData(3);
    
    if ROTATE90
        scalpMesh.vertices = scalpMesh.vertices(:,[2 1 3]);
        csfMesh.vertices   = csfMesh.vertices(:,[2 1 3]);
        cortexMesh.vertices= cortexMesh.vertices(:,[2 1 3]);
        
        % center the head
        scalpMesh.vertices = bsxfun(@minus,scalpMesh.vertices,mean(scalpMesh.vertices));
        csfMesh.vertices = bsxfun(@minus,csfMesh.vertices,mean(csfMesh.vertices));
        cortexMesh.vertices = bsxfun(@minus,cortexMesh.vertices,mean(cortexMesh.vertices));
    end
end
%% initialize some variables
[benchmarking.preproc     ...
 benchmarking.modeling    ...
 benchmarking.brainmovie  ...
 benchmarking.viscsd      ...
 benchmarking.lsl] = deal(NaN);

newPipeline     = true;
newBrainmovie   = true;
specOpts        = [];
gobj            = [];


% %% initialize options
% if ~isempty(GUI_CONFIG_NAME) && exist(env_translatepath([datapath GUI_CONFIG_NAME]),'file')
%     % try to load a configuration file
%     io_load([datapath GUI_CONFIG_NAME]);
% else
%     % manually set some defaults
%     opts.miscOptCfg.dispBenchmark  = true;
%     opts.miscOptCfg.doSIFT         = [];
%     opts.holdPipeline              = false;
%     opts.exitPipeline              = false;
% 
%     opts.siftPipCfg.preproc.verb = 0;
% 
%     % set some defaults
%     opts.fltPipCfg = ...
%         struct('DataCleaning', ...
%                 struct('RetainPhases',true, ...
%                        'HaveBursts', true, ...
%                        'HaveChannelDropouts', false, ...
%                        'DataSetting', 'drycap'), ...
%                 'Rereferencing', {{}}); 
% 
%     opts.siftPipCfg = ...
%         struct('preproc', ...
%                     struct('normalize', ...
%                         struct('verb', 0, ...
%                                'method', {{'time'}})));
%                            
% end

% wait for an opts structure to arrive
fprintf('Waiting for a connection...\n');
haveDataChunk  = false;
haveOptsStruct = false;
while 1
    [newchunk,age] = onl_lslreceive(stream_inlet,-1);
    if isfield(newchunk,'data')
        % this is an EEG chunk
        eeg_chunk = newchunk;
        haveDataChunk = true;
    elseif isfield(newchunk,'siftPipCfg')
        % this is an opts chunk
        opts = newchunk;
        % pipeline changed
        newPipeline    = true;
        haveOptsStruct = true;
    end
    
    if haveDataChunk && haveOptsStruct
        break;
    end
    pause(0.1);
end
start(tmr_checkmem);
fprintf('Connected!\n');

% additional options regarding brainmovie node/edge adaptation
opts.adaptation.adaptationHL    = 10;  % HL of exp. win. MA in frames
opts.adaptation.updateInterval  = 1;
opts.adaptation.bufferTime      = 1;

% set the arg_direct flag to true
opts = arg_setdirect(opts,true);

%% initialize the Meta Control Panel
% if ~isempty(HEAD_MODEL_NAME)
%     eeg_chunk.srcpot = 1; % allow sources
% end
figHandles.MetaControlPanel = gui_metaControlPanel(eeg_chunk,opts);
set(figHandles.MetaControlPanel,'Name','MetaControlPanel: Visualizer');
waitfor(figHandles.MetaControlPanel,'UserData','initialized');

%% Initialize figures
figHandles.specDisplay      = [];
figHandles.BMDisplay        = [];
figHandles.BMControlPanel   = [];
figHandles.csdDisplay       = [];

%% try to load a BrainMovie configuration
try 
    io_load([datapath GUI_BRAINMOVIE_CONFIG_NAME]);
    BMCFG.BMopts.bmopts_suppl = {'title',{'Multivariate ','Granger Causality  '}};
    BMCFG.BMopts.caption = true;
catch
    BMCFG = struct([]);
end

BMRenderMode = 'init_and_render';

% Only run a pass through calibration data
if TRAIN_ONLY
    cleaned_data = exp_eval(flt_pipeline('signal',eeg_chunk,opts.fltPipCfg));
    return;
end

%% Main loop
% -------------------------------------------------------------------------
while ~opts.exitPipeline
    
    try
        
         if opts.holdPipeline
            % pause the pipeline
            pause(1);
            continue;
        end
        
        if newPipeline
            % upate the control panel
            fprintf('Pipeline changed\n');
            figHandles.MetaControlPanel = gui_metaControlPanel(eeg_chunk,opts);
            set(figHandles.MetaControlPanel,'Name','MetaControlPanel: Visualizer');
            newPipeline  = false;
        end
        
        % grab a chunk of data from the LSL stream
        % -----------------------------------------------------------------
        tmr_age = tic;
        [newchunk,age,ndiscarded] = onl_lslreceive(stream_inlet,-1,MAX_AGE);
        benchmarking.lsl = toc(tmr_age);
        if ndiscarded>0
            fprintf('Discarded %d samples\n',ndiscarded); end
        if isfield(newchunk,'data')
            % this is an EEG chunk
            eeg_chunk = newchunk;
            if eeg_chunk.pnts == 1
                eeg_chunk.xmin = age;
                eeg_chunk.xmax = age;
            end
        elseif isfield(newchunk,'siftPipCfg')
            % this is an opts chunk
            
            % override some configs using the values in this pipeline
            opts.siftPipCfg = newchunk.siftPipCfg;
            opts.fltPipCfg  = newchunk.fltPipCfg;
            opts.holdPipeline = newchunk.holdPipeline;
            opts.exitPipeline = opts.exitPipeline || newchunk.exitPipeline;
            opts.miscOptCfg.winLenSec = newchunk.miscOptCfg.winLenSec;
            
            % pipeline changed
            newPipeline = true;
            continue;
        end
        
        % visualize the current source density
        % -----------------------------------------------------------------
        if opts.fltPipCfg.psourceLocalize.arg_selection && opts.miscOptCfg.dispCSD.arg_selection
            viscsdbench = tic;
            if size(eeg_chunk.data,2) ~= size(eeg_chunk.srcpot_all,2)
                keyboard;
            end
            gobj = vis_csd(opts.miscOptCfg.dispCSD,'hmObj',eeg_chunk.hmObj,'signal',eeg_chunk,'gobj',gobj,'cortexMesh',eeg_chunk.dipfit.reducedMesh);
            benchmarking.viscsd = toc(viscsdbench);
            figHandles.csdDisplay = gobj.hFigure;
        else
            if ishandle(gobj)
                gobj = []; end
            benchmarking.viscsd = NaN;
        end
        
        % model the chunk via sift
        % ---------------------------------------------------------------------
        if opts.miscOptCfg.doSIFT.arg_selection
            if eeg_chunk.pnts == 0
                fprintf('No data!\n');
                pause(0.01);
                continue;
            end
                
            % visualize the brainmovie
            % -------------------------------------------------------------
            if opts.miscOptCfg.doSIFT.dispBrainMovie.arg_selection

                if newBrainmovie ...
                    | isempty(ishandle(figHandles.BMDisplay)) ...
                    | ~ishandle(figHandles.BMDisplay)

                    % initialize brainmovie
                    if ~isempty(BMCFG)
                        BMCFG.showNodeLabels.nodelabels = eeg_chunk.CAT.curComponentNames;
                        figHandles.BMControlPanel = gui_causalBrainMovie3D_online(eeg_chunk,eeg_chunk.CAT.Conn,struct('arg_direct',0),BMCFG,'timeRange',[]);
                    else
                        figHandles.BMControlPanel = gui_causalBrainMovie3D_online(eeg_chunk,eeg_chunk.CAT.Conn,struct('arg_direct',0),'showNodeLabels',{'nodelabels',eeg_chunk.CAT.curComponentNames});
                    end
                    waitfor(figHandles.BMControlPanel,'UserData','newBrainmovie');

                    % set the closing behavior
                    set(figHandles.BMControlPanel,'CloseRequestFcn','evalin(''base'',''figHandles.BMControlPanel = [];''); delete(gcbf);');

                    BM_CFG_CHANGED       = true;
                    newBrainmovie        = false;
                    figHandles.BMDisplay = [];
                    BMStateVars          = [];
                end

                if BM_CFG_CHANGED
                    BMRenderMode    = 'init_and_render';
                    BMStateVars     = [];
                    
                    % determine whether to adapt limits
                    adapt_nodeSizeDataRange  = isempty(BMCFG.BMopts.graphColorAndScaling.nodeSizeDataRange)  && opts.miscOptCfg.doSIFT.dispBrainMovie.adaptBrainMovieLimits;
                    adapt_edgeSizeDataRange  = isempty(BMCFG.BMopts.graphColorAndScaling.edgeSizeDataRange)  && opts.miscOptCfg.doSIFT.dispBrainMovie.adaptBrainMovieLimits;
                    adapt_nodeColorDataRange = isempty(BMCFG.BMopts.graphColorAndScaling.nodeColorDataRange) && opts.miscOptCfg.doSIFT.dispBrainMovie.adaptBrainMovieLimits;
                    adapt_edgeColorDataRange = isempty(BMCFG.BMopts.graphColorAndScaling.edgeColorDataRange) && opts.miscOptCfg.doSIFT.dispBrainMovie.adaptBrainMovieLimits;
                    
                    % clear state
                    nodeSizeState   = [];
                    nodeColorState  = [];
                    edgeSizeState   = [];
                    edgeColorState  = [];
                    
                    
                    % reset flag
                    BM_CFG_CHANGED = false;
                end
    
                % Handle special rendering of source meshes and colors
                % ---------------------------------------------------------
                if COLOR_SOURCE_ROI 
                    BG_COLOR = [0.1 0.1 0.1];
                    
                    if BMCFG.BMopts.Layers.custom.arg_selection
                        BMCFG.BMopts.Layers.custom.volumefile = eeg_chunk.dipfit.surfmesh;
                        BMCFG.BMopts.Layers.custom.meshcolor  = hlp_getROIVertexColorTable( ...
                                size(eeg_chunk.dipfit.surfmesh.vertices,1), ...
                                eeg_chunk.roiVertices,BG_COLOR,@(x)distinguishable_colors(x,BG_COLOR));
                    end
                    
                    % if the cortex mesh is a custom mesh with 'constant'
                    % coloring...
                    if strcmpi(BMCFG.BMopts.Layers.cortex.cortexres,'custommesh') ...
                       && strcmp(BMCFG.BMopts.Layers.cortex.cortexcolor.arg_selection,'Constant')
                            % ... replace the colors with ROI coloring
                            BMCFG.BMopts.Layers.cortex.cortexcolor.colormapping = ...
                                hlp_getROIVertexColorTable( ...
                                    size(eeg_chunk.dipfit.surfmesh.vertices,1), ...
                                    eeg_chunk.roiVertices,BG_COLOR,@(x)distinguishable_colors(x,BG_COLOR));
%                         {hmObj.atlas.label, hlp_getROIColorTable(hmObj.atlas.label,eeg_chunk.roiLabels,[0.5 0.5 0.5],[1 0 0])};
                    end
                end
                
                % Render the Brain Movie
                % ---------------------------------------------------------
                bmbench = tic;
                [BMcfg_tmp BMhandles BMout] = ...
                    vis_causalBrainMovie3D(eeg_chunk,eeg_chunk.CAT.Conn,BMCFG, ...
                        'timeRange',[], ...
                        'BMopts',hlp_mergeVarargin(BMCFG.BMopts, ...
                                    'figurehandle',figHandles.BMDisplay, ...
                                    'RenderMode',BMRenderMode, ...
                                    'speedy',true, ...
                                    'InternalStateVariables',BMStateVars));
                benchmarking.brainmovie = toc(bmbench);

                if isempty(BMStateVars)
                    set(figHandles.BMDisplay, 'MenuBar','figure', ...
                                'CloseRequestFcn','evalin(''base'',''figHandles.BMDisplay = [];''); delete(gcbf);', ...
                                'ToolBar','none', 'Name','BrainMovie3D');
                end
                    
                % save figure handle and state
                figHandles.BMDisplay = BMhandles.figurehandle;
                BMStateVars          = BMcfg_tmp.BMopts.vars;
                BMRenderMode         = 'render';

                % adapt the node and edge color/size limits
                % ---------------------------------------------------------
                if opts.miscOptCfg.doSIFT.dispBrainMovie.arg_selection ...
                   && opts.miscOptCfg.doSIFT.dispBrainMovie.adaptBrainMovieLimits
               
                    % adapt BrainMovie node/edge limits
                    if adapt_nodeSizeDataRange
                        [BMCFG.BMopts.graphColorAndScaling.nodeSizeDataRange, nodeSizeState] ...
                            = hlp_adaptMinMaxLimits('values',   BMout.NodeSize,     ...
                                                    'state',    nodeSizeState,      ...
                                                    'adaptOpts',opts.adaptation);
                    end
                    if adapt_nodeColorDataRange
                        [BMCFG.BMopts.graphColorAndScaling.nodeColorDataRange, nodeColorState] ...
                            = hlp_adaptMinMaxLimits('values',   BMout.NodeColor,     ...
                                                    'state',    nodeColorState,      ...
                                                    'adaptOpts',opts.adaptation);
                    end
                    if adapt_edgeSizeDataRange
                        [BMCFG.BMopts.graphColorAndScaling.edgeSizeDataRange, edgeSizeState] ...
                            = hlp_adaptMinMaxLimits('values',   BMout.EdgeSize,     ...
                                                    'state',    edgeSizeState,      ...
                                                    'adaptOpts',opts.adaptation);
                    end
                    if adapt_edgeColorDataRange
                        [BMCFG.BMopts.graphColorAndScaling.edgeColorDataRange, edgeColorState] ...
                            = hlp_adaptMinMaxLimits('values',   BMout.EdgeColor,     ...
                                                    'state',    edgeColorState,      ...
                                                    'adaptOpts',opts.adaptation);
                    end
                    
                end % adaptBrainMovieLimits
                
            else % dispBrainMovie == false
                % close the brainmovie figure if control panel is closed
                try
                    if isempty(figHandles.BMControlPanel) ...
                            || ~ishandle(figHandles.BMControlPanel) ...
                            && ishandle(figHandles.BMDisplay)
                        close(figHandles.BMDisplay);
                    end
                catch
                end
                benchmarking.brainmovie = NaN;
            end
                
            % display the VAR spectrum
            % -------------------------------------------------------------
            if opts.miscOptCfg.doSIFT.dispSpectrum
                try
                    if ishandle(figHandles.specDisplay)
                        % update the SpecViewer
                        vis_autospectrum(specOpts,'FigureHandle',figHandles.specDisplay,'stream',eeg_chunk);
                    else
                        % generate a GUI and initialize SpecViewer
                        [figHandles.specDisplay specOpts] = arg_guidialog(@vis_autospectrum,'Parameters',{'stream',eeg_chunk});

                        set(figHandles.specDisplay,'CloseRequestFcn','evalin(''base'',''opts.miscOptCfg.dispSpectrum=false; figHandles.specDisplay = [];''); delete(gcbf);');
                        specOpts.arg_direct = 0;
                    end
                catch e
                    disp(e.message);
                    opts.miscOptCfg.doSIFT.dispSpectrum = false;
                end
            else
                if ishandle(figHandles.specDisplay)
                    close(figHandles.specDisplay);
                end
            end

        else
            benchmarking.modeling = NaN;
        end
        
        if opts.miscOptCfg.dispBenchmark
            % TODO: implement this as a figure
            vis_benchmark(benchmarking);
        end

    catch err
         hlp_handleerror(err);
    end
   
    pause(0.0001);
    
end

%% all done, close all open figures
fnames = fieldnames(figHandles);
for f=1:length(fnames)
    if ishandle(figHandles.(fnames{f}))
        close(figHandles.(fnames{f}));
    end
end

%% close any streams
stream_inlet.close_stream;
stream_inlet.delete;
if ~isempty(MIN_MEM_LIMIT)
    stop(tmr_checkmem);
    delete(tmr_checkmem);
end