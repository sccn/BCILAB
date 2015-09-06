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

%% 
    
GUI_CONFIG_NAME = ''; 'DARPA_DEMO_METACP_CFG_LATEST.mat'; %'MetaCPopts_Pyramind.mat'; %MetaCPopts_Pyramind.mat'; %'MetaCPopts.mat';
GUI_BRAINMOVIE_CONFIG_NAME = 'DARPA_DEMO_BM_CFG.mat'; %'BMCFG.mat'

%% establish where we will read/write prefs, etc from
datapath = 'data:/';  % this is relative to the BCILAB root dir
TrainingDataFile = 'christian-calib4.xdf'; %'Cognionics_Pyramind_demo.set'; %'clean_reversed.xdf'; %'noisy.xdf'; %'Cognionics_Pyramind_demo.set';
TestingDataFile =  'christian-calib.xdf'; %'Cognionics_Pyramind_demo.set'; %'clean_reversed.xdf'; %'noisy.xdf'; %'Cognionics_Pyramind_demo.set';

%% Set up the name of the stream we will write to in the workspace
streamname              = 'EEGDATA';
LSL_SelectionProperty   = 'type';  % 'name'
LSL_SelectionValue      = 'EEG';   % 'Cogionics'
 
%% load calibration recording
raw = exp_eval(io_loadset([datapath TrainingDataFile],'markerchannel',{'remove_eventchns',false}));
flt_pipeline('update');

%% start LSL streaming
% run_readlsl('MatlabStream',streamname,'SelectionProperty',LSL_SelectionProperty,'SelectionValue',LSL_SelectionValue);

%% ... OR read from datafile
run_readdataset('MatlabStream',streamname,'Dataset',io_loadset([datapath TestingDataFile],'markerchannel',{'remove_eventchns',false}));

%% initialize some variables
[benchmarking.preproc     ...
 benchmarking.modeling    ...
 benchmarking.brainmovie] = deal(NaN);

newPipeline     = true;
newBrainmovie   = true;
specOpts = [];

%% initialize options
    
if ~isempty([datapath GUI_CONFIG_NAME]) && ~any(exist(env_translatepath([datapath GUI_CONFIG_NAME]),'file') == [0 7])
    % try to load a configuration file
    io_load([datapath GUI_CONFIG_NAME]);
else
    % manually set some defaults
    opts.miscOptCfg.dispBenchmark  = true;
    opts.miscOptCfg.doSIFT    = []; %struct(...
                                    %'dispBrainMovie',{'adaptBrainMovieLimits', true}, ...
                                    %'dispSpectrum',true, ...
                                    %'winLenSec',2);
    opts.holdPipeline   = false;
    opts.exitPipeline   = false;

    opts.siftPipCfg.preproc.verb = 0;

    % set some defaults
    opts.fltPipCfg = ...
        struct('DataCleaning', ...
                struct('RetainPhases',true, ...
                       'HaveBursts', true, ...
                       'HaveChannelDropouts', false, ...
                       'DataSetting', 'drycap'), ...
                'Rereferencing', {{}}); 

               
    opts.siftPipCfg = ...
        struct('preproc', ...
                    struct('normalize', ...
                        struct('verb', 0, ...
                               'method', {{'time'}})));
                           
end

% additional options regarding brainmovie node/edge adaptation
opts.adaptation.adaptationHL = 10;  % HL of exp. win. MA in frames
opts.adaptation.updateInterval = 1;
opts.adaptation.bufferTime = 1;


%% initialize the Meta Control Panel
figHandles.MetaControlPanel = gui_metaControlPanel(raw,opts);
waitfor(figHandles.MetaControlPanel,'UserData','initialized');

% request if you want to save the results. It will be saved to
 % datapath folder above
% if strcmp(questdlg('Do you want to save this configuration?','Save config?','Yes','No','No'),'Yes')
%     if isempty(GUI_CONFIG_NAME)
%         GUI_CONFIG_NAME = uiputfile;
%     end
%     io_save([datapath GUI_CONFIG_NAME],'opts');
% end

%% Initialize figures
figHandles.specDisplay      = [];
figHandles.BMDisplay        = [];
figHandles.BMControlPanel   = [];

%% try to load a BrainMovie configuration
try 
    io_load([datapath GUI_BRAINMOVIE_CONFIG_NAME]);
    BMCFG.BMopts.bmopts_suppl = {'title',{'Multivariate ','Granger Causality  '}};
    BMCFG.BMopts.caption = true;
catch
    BMCFG = struct([]);
end

BMRenderMode = 'init_and_render';

%% Main loop

while ~opts.exitPipeline
    
    try
    
        if newPipeline
            % create a new pipeline
            fprintf('Pipeline changed\n');
%             keyboard;
            cleaned_data = exp_eval(flt_pipeline('signal',raw,opts.fltPipCfg));
            pipeline = onl_newpipeline(cleaned_data,streamname);
            newPipeline = false;
        end

        if opts.holdPipeline
            % pause the pipeline
            pause(1);
            continue;
        end

        % grab a chunk of data from the stream
        % ---------------------------------------------------------------------
        tic;
        [eeg_chunk,pipeline] = onl_filtered(pipeline,round(opts.miscOptCfg.winLenSec*cleaned_data.srate));
        benchmarking.preproc = toc;

        % process the chunk via sift
        % ---------------------------------------------------------------------
        if opts.miscOptCfg.doSIFT.arg_selection
            tic;
            try
                % force the modeling window length to agree with chunk length
                opts.siftPipCfg.modeling.winlen = opts.miscOptCfg.winLenSec;
                % select a subset of channels
                if ~isempty(opts.miscOptCfg.doSIFT.channelSubset)
                    [dummy eeg_chunk] = evalc('hlp_scope({''disable_expressions'',1},@flt_selchans,''signal'',eeg_chunk,''channels'',opts.miscOptCfg.doSIFT.channelSubset,''arg_direct'',true)');
                    %[dummy eeg_chunk] = evalc('exp_eval(flt_selchans(''signal'',eeg_chunk,''channels'',opts.miscOptCfg.doSIFT.channelSubset))');
                end
                % run the SIFT modeling pipeline
                [eeg_chunk cfg] = hlp_scope({'is_online',1},@onl_siftpipeline,'EEG',eeg_chunk,opts.siftPipCfg);
                % [eeg_chunk cfg] = onl_siftpipeline('EEG',eeg_chunk,opts.siftPipCfg);
                
                if eeg_chunk.pnts == 0
                    fprintf('No data!\n');
                    pause(0.01);
                    continue;
                end
                
            catch err
                hlp_handleerror(err);
                pause(0.1);
                continue;
            end
            benchmarking.modeling = toc;
        


            % visualize the brainmovie
            % ---------------------------------------------------------------------
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
                    set(figHandles.BMControlPanel, ...
                        'CloseRequestFcn', ...
                        'evalin(''base'',''figHandles.BMControlPanel = [];''); delete(gcbf);');

                    BM_CFG_CHANGED = true;
                    newBrainmovie = false;
                    figHandles.BMDisplay = [];
                    BMStateVars = [];
                end

                if BM_CFG_CHANGED
                    BMRenderMode = 'init_and_render';
                    BM_CFG_CHANGED = false;
                    BMStateVars = [];

                    % reset the limits
                    numberOfRunsSoFar = 0;
                    
                    adapt_nodeSizeDataRange  = isempty(BMCFG.BMopts.graphColorAndScaling.nodeSizeDataRange);
                    adapt_edgeSizeDataRange  = isempty(BMCFG.BMopts.graphColorAndScaling.edgeSizeDataRange);
                    adapt_nodeColorDataRange = isempty(BMCFG.BMopts.graphColorAndScaling.nodeColorDataRange);
                    adapt_edgeColorDataRange = isempty(BMCFG.BMopts.graphColorAndScaling.edgeColorDataRange);
                end
                
                
                
                % update the scaling limits for brainmovie
                if opts.miscOptCfg.doSIFT.dispBrainMovie.arg_selection ...
                   && opts.miscOptCfg.doSIFT.dispBrainMovie.adaptBrainMovieLimits ...
                   && numberOfRunsSoFar > 0
               
                    if adapt_nodeSizeDataRange
                        BMCFG.BMopts.graphColorAndScaling.nodeSizeDataRange  = [lastNodeSizeMin lastNodeSizeMax];
                    end
                    if adapt_edgeSizeDataRange
                        BMCFG.BMopts.graphColorAndScaling.edgeSizeDataRange  = [lastEdgeSizeMin lastEdgeSizeMax];
                    end
                    if adapt_nodeColorDataRange
                        BMCFG.BMopts.graphColorAndScaling.nodeColorDataRange = [lastNodeColorMin lastNodeColorMax];
                    end
                    if adapt_edgeColorDataRange
                        BMCFG.BMopts.graphColorAndScaling.edgeColorDataRange = [lastEdgeColorMin lastEdgeColorMax];
                    end
                else
                    BMCFG.BMopts.graphColorAndScaling.nodeSizeDataRange  = [];
                    BMCFG.BMopts.graphColorAndScaling.edgeSizeDataRange  = [];
                    BMCFG.BMopts.graphColorAndScaling.nodeColorDataRange = [];
                    BMCFG.BMopts.graphColorAndScaling.edgeColorDataRange = [];

                    % initialize
                    [lastNodeSizeMin lastNodeSizeMax ...
                     lastEdgeSizeMin lastEdgeSizeMax ...
                     lastNodeColorMin lastNodeColorMax ...
                     lastEdgeColorMin lastEdgeColorMax] = deal(0);
                end

                % update the BrainMovie3D
                tic
                [BMcfg_tmp BMhandles BMout] = ...
                    vis_causalBrainMovie3D(eeg_chunk,eeg_chunk.CAT.Conn,BMCFG, ...
                        'timeRange',[], ...
                        'BMopts',hlp_mergeVarargin(BMCFG.BMopts, ...
                                    'figurehandle',figHandles.BMDisplay, ...
                                    'RenderMode',BMRenderMode, ...
                                    'speedy',true, ...
                                    'InternalStateVariables',BMStateVars));
                benchmarking.brainmovie = toc;

                % save figure handle and state
                figHandles.BMDisplay = BMhandles.figurehandle;
                BMStateVars          = BMcfg_tmp.BMopts.vars;
                BMRenderMode         = 'render';

                % adapt the node and edge color/size limits
                if opts.miscOptCfg.doSIFT.dispBrainMovie.arg_selection ...
                   && opts.miscOptCfg.doSIFT.dispBrainMovie.adaptBrainMovieLimits
               
                    [lastNodeSizeMin lastNodeSizeMax] ...
                        = hlp_scaleLimits(BMout.NodeSize,   ...
                                        lastNodeSizeMin,    ...
                                        lastNodeSizeMax,    ...
                                        numberOfRunsSoFar,  ...
                                        opts.adaptation);

                    [lastEdgeSizeMin lastEdgeSizeMax] ...
                        = hlp_scaleLimits(BMout.EdgeSize,   ...
                                        lastEdgeSizeMin,    ...
                                        lastEdgeSizeMax,    ...
                                        numberOfRunsSoFar,  ...
                                        opts.adaptation);

                    [lastNodeColorMin lastNodeColorMax] ...
                        = hlp_scaleLimits(BMout.NodeColor,  ...
                                        lastNodeColorMin,   ...
                                        lastNodeColorMax,   ...
                                        numberOfRunsSoFar,  ...
                                        opts.adaptation);

                    [lastEdgeColorMin lastEdgeColorMax] ...
                        = hlp_scaleLimits(BMout.EdgeColor,  ...
                                        lastEdgeColorMin,   ...
                                        lastEdgeColorMax,   ...
                                        numberOfRunsSoFar,  ...
                                        opts.adaptation);
                end

                numberOfRunsSoFar = numberOfRunsSoFar + 1;

                if numberOfRunsSoFar == 1
                    set(figHandles.BMDisplay, 'MenuBar','figure', ...
                        'CloseRequestFcn','evalin(''base'',''figHandles.BMDisplay = [];''); delete(gcbf);', ...
                        'ToolBar','none', 'Name','BrainMovie3D');
                end

            else
                % close the brainmovie figure as well
                if isempty(figHandles.BMControlPanel) ...
                        || ~ishandle(figHandles.BMControlPanel) ...
                        && ishandle(figHandles.BMDisplay)
                    close(figHandles.BMDisplay);
                end
            end

            % display the VAR spectrum
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

        end
        
        if opts.miscOptCfg.dispBenchmark
            % TODO: implement this as a figure
            vis_benchmark(benchmarking);
        end

    catch err
         hlp_handleerror(err);
    end
    
    pause(0.005);
    
end

%% all done, close all open figures
fnames = fieldnames(figHandles);
for f=1:length(fnames)
    if ishandle(figHandles.(fnames{f}))
        close(figHandles.(fnames{f}));
    end
end

