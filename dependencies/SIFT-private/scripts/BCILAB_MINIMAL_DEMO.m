%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REAL-TIME BCILAB-LSL DEMONSTRATION SCRIPT
%
% This script demonstrates how to set up a default BCILAB-LSL pipeline for 
% real-time data extraction and processing from an EEG device or file
%
% This script requires the following tools to be installed:
% BCILAB: https://www.sccn.ucsd.edu/wiki/BCILAB
% LSL:    https://code.google.com/p/labstreaminglayer/
% 
% Optional (if using SourceLocalization)
% MobiLab: http://sccn.ucsd.edu/wiki/MoBILAB
%
% Author: Tim Mullen,  October 14, 2013, SCCN/INC/UCSD
% Email:  tmullen@ucsd.edu
%
% .........................................................................
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
% .........................................................................
%
% If you use this script as a basis for your published work, kindly 
% consider acknowledging us in your manuscript. A possible example:
% "Thanks to Tim Mullen and Christian Kothe (SCCN, UCSD) for providing 
%  scripting support"
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Force EEGLAB to use double precision
pop_editoptions( 'option_single', false);

%% Define config options (USER MUST DEFINE THIS)
% -------------------------------------------------------------------------
% If runlsl = true, then stream 'online' from device
% If runlsl = false then stream data from opts.PlaybackDataFile (below)
opts.runlsl = false;           

% Run SIFT on the pipeline (following BCILAB filter pipeline)
opts.runSIFT = true;

% Set up the name of the LSL stream
opts.lsl.StreamName          = 'EEGDATA';
opts.lsl.SelectionProperty   = 'type';  % can also be 'name'
opts.lsl.SelectionValue      = 'EEG';

% Determine a fixed 'sliding' window length (seconds) for pulling data from LSL
% If winlen=0, pull the most recent (variable length) chunk of unprocessed data 
opts.winlen = 1;

% suppress console output in processing loop?
% if false, this can reusult in a LOT of text printed to console
opts.silence = true;

% Establish file paths
% ..........................................................................
% User Tip:
%
% All paths are relative to opts.datapath. This is platform-independent and 
% itself can be relative to bcilab root folder (i.e. data:/ is the userdata 
% folder in the bcilab root dir)
% ..........................................................................
opts.datapath         = 'data:/';
% (required) calibration data set2
opts.TrainingDataFile = 'Simulation_EEGLAB_XVII.set'; %'MyCalibrationFile.set'; 
% (optional) file to "playback" if opts.runlsl = false
opts.PlaybackDataFile = 'Simulation_EEGLAB_XVII.set'; %'MyPlaybackFile.set';      
% (optional) pipeline config file
opts.BCILAB_PipelineConfigFile = ''; 'MyFilteringPipeline.mat';  
% (option) SIFT pipeline config file
opts.SIFT_PipelineConfigFile = ''; 'MySiftFilteringPipeline.mat';  
% (optional) name of file containing MobiLab head model
opts.HeadModelDataFile= ''; %'NameOfHeadModelFile.mat';   % Look in <SIFT_RootFolder>\resources\headmodels\standard-Colin27-385ch.mat

%% Load head model object
% -------------------------------------------------------------------------
if ~isempty(opts.HeadModelDataFile)
    hmObj = headModel.loadFromFile(env_translatepath([opts.datapath opts.HeadModelDataFile]));
end

%% load calibration recording
% -------------------------------------------------------------------------
calibData = exp_eval_optimized(io_loadset([opts.datapath opts.TrainingDataFile], ...
                         'markerchannel',{'remove_eventchns',false}));
flt_pipeline('update');

%% Start LSL streaming
% -------------------------------------------------------------------------
if opts.runlsl == true
    run_readlsl('MatlabStream',opts.lsl.StreamName,                 ...
                'SelectionProperty',opts.lsl.SelectionProperty,     ...
                'SelectionValue',opts.lsl.SelectionValue);
else
    % ... OR read from datafile
    run_readdataset('MatlabStream',opts.lsl.StreamName,                         ...
                    'Dataset',io_loadset([opts.datapath opts.PlaybackDataFile], ...
                    'markerchannel',{'remove_eventchns',false}));
end

%% Set up the pipeline
% -------------------------------------------------------------------------
% load existing config file(s) (if it exists)
try    fltPipCfg = exp_eval(io_load([opts.datapath opts.BCILAB_PipelineConfigFile])); 
catch, disp('-- no existing pipeline --'); fltPipCfg = {}; end
try    siftPipCfg = exp_eval(io_load([opts.datapath opts.SIFT_PipelineConfigFile])); 
catch, disp('-- no existing pipeline --'); siftPipCfg = {}; end

% enforce usage of user-defined head model
if ~isempty(opts.HeadModelDataFile)
    fltPipCfg.psourceLocalize.hmObj = [opts.datapath opts.HeadModelDataFile]; end

% display GUI to allow user to configure the BCILAB pipeline...
fltPipCfg = arg_guipanel('Function',@flt_pipeline, ...
                         'Parameters',[{'signal',calibData} fltPipCfg], ...
                         'PanelOnly',false);

% display GUI to allow user to configure the SIFT pipeline...
if opts.runSIFT
    if isfield(calibData,'CAT')
        calibData = rmfield(calibData,'CAT');
    end
    siftPipCfg = arg_guipanel('Function',@onl_siftpipeline, ...
                             'Parameters',[{'EEG',calibData} siftPipCfg], ...
                             'PanelOnly',false);
    % force the modeling window length to agree with chunk length
    siftPipCfg.modeling.winlen = opts.winlen;
end
% ..........................................................................
% User Tip:
% 
% You might want to save the pipeline config (fltPipCfg) at this point for
% later reload (from 'opts.BCILAB_PipelineConfigFile')
% save(env_translatepath([opts.datapath opts.SIFT_PipelineConfigFile]),'-struct','siftPipCfg');
% save(env_translatepath([opts.datapath opts.BCILAB_PipelineConfigFile]),'-struct','fltPipCfg');
% ..........................................................................

% ...apply the pipeline to calibration data    
disp('-- Calibrating pipeline on training data (please wait) --');
cleaned_data = exp_eval(flt_pipeline('signal',calibData,fltPipCfg));

% ..........................................................................
% User Tip:
%
% Instead of using the GUI, you can also easily define pipelines in the 
% script itself. Here are some examples:
%
% % remove flatline channels:
% noflats      = flt_clean_flatlines(flt_seltypes(calibData,'EEG'));
% % apply a high-pass filter:
% highpassed   = flt_fir(noflats,[0.5 1],'highpass');
% % discard noisy channels:
% goodchannels = flt_clean_channels('signal',highpassed,'min_corr',0.35,'ignored_quantile',0.3);
%
% % you always need to make sure to apply the pipeline
% cleaned_data = exp_eval(goodchannels);  
% ..........................................................................


% ..........................................................................
% User Tip:
%
% After this step, cleaned_data will contain the results of having applied
% the entire pipeline to your calibration data. Thus, the steps up to this 
% point essentially define an offline data processing pipeline which you
% can use for your publications. It is also helpful at this stage to
% examine the contents of cleaned_data for deficiencies, which can help you 
% fine-tune your pipeline for later online use.
% ..........................................................................

% initialize the pipeline for streaming data
pipeline     = onl_newpipeline(cleaned_data,opts.lsl.StreamName);

% render a panel for viewing streams
gui_vis_filtered;

%% Main loop
% -------------------------------------------------------------------------
disp('-- Running pipeline on playback data --');
chunk_len = round(opts.winlen*cleaned_data.srate);
ConnData = [];
while true
    % grab a chunk of data from the stream and send through the filter pipeline
    [eeg_chunk,pipeline] = onl_filtered(pipeline, chunk_len, opts.silence);
    
    % process data chunk with SIFT
    if opts.runSIFT
        [eeg_chunk cfg] = hlp_scope({'is_online',1},@onl_siftpipeline,'EEG',eeg_chunk,siftPipCfg);
    end
    
    % .....................................................................
    % User Tip:
    % 
    % The structure 'eeg_chunk' now contains the results of applying the 
    % pipeline to the chunk of raw data acquired from the LSL stream. 
    % This is an augmented EEGLAB data structure.
    % 
    % - You can add code here to apply your own additional processing steps
    % - You'll also want to add code to store the results you'd like keep
    %
    % Hint:
    %  eeg_chunk.data contains processed/cleaned data
    %
    %  (if using source localization) 
    %  eeg_chunk.srcpot contains current density for ROIs
    %
    %  (if using source localization and keepFullCsd is enabled)
    %  eeg_chunk.srcpot_all contains current density for all mesh vertices
    %
    %  (if using ICA)
    %  eeg_chunk.icaact contains ica activations
    % 
    % (if using SIFT)
    % eeg_chunk.CAT  contains a SIFT object/datastructure
    %  CAT
    %   .Conn       Connectivity object
    %       .winCenterTimes
    %       .freqs
    %       .dims   dimensions for connectivity matrix 
    %               e.g. (chans_to, chans_from, freqs, time)
    %       .<measure_name>    Connectivity matrix (conforms to dims above)
    %   .MODEL      Model structure
    % 
    % .....................................................................
    
    
    % (optionally) give Matlab a little break...
    pause(0.0001);
end

