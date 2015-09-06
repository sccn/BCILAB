function [EEG g] = pre_prepData(varargin)
%
% Preprocess EEG dataset(s) for connectivity analysis. See [1] for
% mathematical details on preprocessing steps.
%
% Optional                      Information
% ------------------------------------------------------------------------------------------------------------------------------------
% VerbosityLevel:               Verbosity level. 0 = no output, 1 = text, 2 = graphical
%                               Possible values: 0,1,2
%                               Default value  : 2
%                               Input Data Type: real number (double)
%
% SignalType:                   Type of signal to analyze
%                               Possible values: 'Components','Channels'
%                               Default value  : 'Components'
%                               Input Data Type: string
%
%     | ConvertChanlocs2Dipfit: Create dipfit structure from chanlocs
%                               If channel locations are available and coregistered to an MRI (e.g. the MNI brain), this will
%                               construct a "source" dipfit structure which in turn will enable 3D visualization of network
%                               structure.
%                               NOTE: this will overwrite any existing dipfit structure in the current EEG set.
%                               Input Range  : Unrestricted
%                               Default value: 1
%                               Input Data Type: boolean
%     |
%     |     | chanlocs:         Channel locations structure
%                               Input Data Type: string
%     |
%     |     | ScalingFactor:    Scaling factor
%                               Each [X,Y,Z] channel coordinate will be multiplied by this factor. This is useful for converting
%                               between unit types (e.g. meters -> mm)
%                               Input Range  : Unrestricted
%                               Default value: 1000
%                               Input Data Type: real number (double)
%
% DifferenceData:               Differencing options
%                               A difference filter of order k applied to time series X is defined as [for i=1:k, for t=1:N, X(t) =
%                               X(t)-X(t-1)]
%                               Input Range  : Unrestricted
%                               Default value: 0
%                               Input Data Type: boolean
%
%     | VerbosityLevel:         Verbosity level. 0 = no output, 1 = text, 2 = graphical
%                               Possible values: 0,1,2
%                               Default value  : 1
%                               Input Data Type: real number (double)
%
%     | DifferencingOrder:      Differencing order
%                               Number of times to difference data
%                               Input Range  : [0  10]
%                               Default value: 1
%                               Input Data Type: real number (double)
%
% Detrend:                      Detrend or center each epoch
%                               Input Range  : Unrestricted
%                               Default value: 0
%                               Input Data Type: boolean
%
%     | SamplingRate:           Sampling Rate
%                               Input Data Type: string
%
%     | VerbosityLevel:         Verbosity level. 0 = no output, 1 = text, 2 = graphical
%                               Possible values: 0,1,2
%                               Default value  : 1
%                               Input Data Type: real number (double)
%
%     | DetrendingMethod:       Detrending options
%                               Linear: removes the least-squares fit of a straight line from each trial. Constant: removes the
%                               mean from each trial (centering)
%                               Possible values: 'linear','constant'
%                               Default value  : 'linear'
%                               Input Data Type: boolean
%
%     | Piecewise:              Use piecewise detrending
%                               Divide the data into (overlapping) segments and detrend each segment separately. Segment endpoints
%                               are merged and stitched together using a cubic spline function to minimize discontinuities at
%                               segment intersection points. This is useful for as an alternative to high-pass filtering for
%                               removing infraslow oscillations (e.g. SC drift) from the
%                               Input Range  : Unrestricted
%                               Default value: 0
%                               Input Data Type: boolean
%     |
%     |     | SegmentLength:    Length of each detrending segment (sec)
%                               Input Range  : [2.22045e-16           10]
%                               Default value: 0.33
%                               Input Data Type: real number (double)
%     |
%     |     | StepSize:         Step size between segment centers (sec)
%                               It is recommended to use at least 0.5*SegmentLength (50% overlap).
%                               Input Range  : [2.22045e-16           10]
%                               Default value: 0.0825
%                               Input Data Type: real number (double)
%
%     | Plot:                   Plot results for inspection
%                               Input Range  : Unrestricted
%                               Default value: 0
%                               Input Data Type: boolean
%
% NormalizeData:                Data normalization
%                               Normalize trials across time, ensemble, or both
%                               Input Range  : Unrestricted
%                               Default value: 1
%                               Input Data Type: boolean
%
%     | VerbosityLevel:         Verbosity level. 0 = no output, 1 = text, 2 = graphical
%                               Possible values: 0,1,2
%                               Default value  : 0
%                               Input Data Type: real number (double)
%
%     | Method:                 Normalize windows across time, ensemble, or both
%                               Possible values: 'time','ensemble'
%                               Default value  : 'time','ensemble'
%                               Input Data Type: boolean
%
%     | VerbosityLevel:         Verbosity level. 0 = no output, 1 = text, 2 = graphical
%                               Possible values: 0,1,2
%                               Default value  : 0
%                               Input Data Type: real number (double)
%
%     | Method:                 Normalize windows across time, ensemble, or both
%                               Possible values: 'time','ensemble'
%                               Default value  : 'time','ensemble'
%                               Input Data Type: boolean
%
% AmplitudeEnvelope:            Compute amplitude envelope
%                               This will add additional 'pseudochannels' which represent the amplitude envelopes of original
%                               channels over a specified frequency band. Envelopes are computed via a zero-phase filter (eegfilt)
%                               followed by a hilbert transform.
%                               Input Range  : Unrestricted
%                               Default value: 0
%                               Input Data Type: boolean
%
%     | AmplitudePassBand:      The [lo hi] pass-band to use for amplitude (Hz)
%                               Input Range  : Unrestricted
%                               Default value: 80  150
%                               Input Data Type: real number (double)
%
%     | NormalizeData:          Data normalization
%                               Normalize analytic amplitude trials across time, ensemble, or both
%                               Input Range  : Unrestricted
%                               Default value: 1
%                               Input Data Type: boolean
%     |
%     |     | VerbosityLevel:   Verbosity level. 0 = no output, 1 = text, 2 = graphical
%                               Possible values: 0,1,2
%                               Default value  : 0
%                               Input Data Type: real number (double)
%     |
%     |     | Method:           Normalize windows across time, ensemble, or both
%                               Possible values: 'time','ensemble'
%                               Default value  : 'ensemble'
%                               Input Data Type: boolean
%     |
%     |     | VerbosityLevel:   Verbosity level. 0 = no output, 1 = text, 2 = graphical
%                               Possible values: 0,1,2
%                               Default value  : 0
%                               Input Data Type: real number (double)
%     |
%     |     | Method:           Normalize windows across time, ensemble, or both
%                               Possible values: 'time','ensemble'
%                               Default value  : 'ensemble'
%                               Input Data Type: boolean
%
%     | PlotData:               Plot amplitude envelope
%                               Input Range  : Unrestricted
%                               Default value: 0
%                               Input Data Type: boolean
%
%     | Verbosity:              Verbose output
%                               Input Range  : Unrestricted
%                               Default value: 1
%                               Input Data Type: boolean
%
% BadDataSegments:              Intervals of bad data
%                               N x 2 matrix of [lo hi] intervals (seconds) of data within ea. trial to set to nan. Currently only
%                               compatible with vierra_morf MVAR method
%                               Input Range  : Unrestricted
%                               Default value: n/a
%                               Input Data Type: real number (double)
%
% BalanceTrials:                Equalize the number of trials between two conditions
%                               Input Range  : Unrestricted
%                               Default value: 0
%                               Input Data Type: boolean
%
% Output:
%
%   EEG:        Prepocessed EEG structure(s). EEG.CAT contains
%                   SIFT-specific data
%   args:           Argument specification structure(s).
%
%
% See Also: pop_pre_prepData()
%
% References:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual. Section 6.5.1
%   Available at: http://www.sccn.ucsd.edu/wiki/Sift
%
% Author: Tim Mullen 2010, SCCN/INC, UCSD.
% Email:  tim@sccn.ucsd.edu

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

g = struct([]);

% extract some defaults
EEG = arg_extract(varargin,{'EEG','ALLEEG'},1);
MyComponentNames = [];
MyChannelNames   = [];
% 
% if ~isempty(EEG)
%     
%     % determine the allowable signal types
%     defSigType = {'Channels'};
%     if ~any(cellfun(@isempty,{EEG.icaweights}))
%         % if EEG contains icaweights, allow components
%         defSigType = [defSigType {'Components'}];
%     end
%     if any(arrayfun(@(x) isfield(x,'srcpot'), EEG))
%         % if EEG contains the field 'srcpot', allow sources
%         defSigType = [defSigType {'Sources'}];
%     end
%     
%     clear EEG;
% else
%     defSigType = {'Channels'};
% end

defSigType = {'Channels','Sources','Components'};

verb = arg_extract(varargin,{'verb','VerbosityLevel'},[],0);

g = arg_define(varargin, ...
    arg_norep({'EEG','ALLEEG'},mandatory,[],'An EEGLAB dataset'), ...
    arg({'verb','VerbosityLevel'},int32(2),{int32(0) int32(1) int32(2)},'Verbosity level. 0 = no output, 1 = text, 2 = graphical'), ...
    arg_subswitch({'sigtype','SignalType','chantype'},hlp_getSigTypeArglist(defSigType,true),hlp_getSigTypeArglist(defSigType),'Type of signal to analyze. If ''Components'', data in EEG.icaact will be processed. If ''Channels'' EEG.data will be processed. If ''Sources'' EEG.srcpot will be processed.','suppress',{'componentsToKeep'}), ...
    arg({'varnames','VariableNames'},[],[],'Optional names for channels/components. This does not affect variable selection. If left empty, component indices or channel names (chanlocs) are used. The number of names provided must equal the number of components/channels selected','type','cellstr','shape','row'), ...
    arg_subtoggle({'diff','DifferenceData'},[],@pre_diffData,'Differencing options. A difference filter of order k applied to time series X is defined as [for i=1:k, for t=1:N, X(t) = X(t)-X(t-1)]','cat','Filtering','suppress',{'verb'}), ...
    arg_subtoggle({'detrend','Detrend'},[],@pre_detrend,'Detrend or center each epoch','cat','Filtering','suppress',{'verb','srate'}), ...
    arg_subtoggle({'normalize','NormalizeData'},{'method', {'time','ensemble'}, 'verb',verb},@pre_normData,'Data normalization. Normalize trials across time, ensemble, or both','cat','Normalization','suppress',{'verb'}), ...
    arg_subtoggle({'aamp','AmplitudeEnvelope'},[],@est_aamp,'Compute amplitude envelope. This will add additional ''pseudochannels'' which represent the amplitude envelopes of original channels over a specified frequency band. Envelopes are computed via a zero-phase filter (eegfilt) followed by a hilbert transform.','cat','Filtering','suppress',{'verb'}), ...
    arg({'resetConfigs','ResetConfigs'},false,[],'Clear existing config structures. These structures are used to set defaults for GUIs. This option is appropriate if you are re-preprocessing (or have modified) a dataset in a way that will fundamentally alter the structure of the dataset'), ...
    arg_nogui({'badsegments','BadDataSegments'},[],[],'Intervals of bad data. N x 2 matrix of [lo hi] intervals (seconds) of data within ea. trial to set to nan. Currently only compatible with vierra_morf MVAR method','cat','Data Selection'), ...
    arg_nogui({'newtrials','TrialsToKeep','TrialSubsetToUse'},[],[],'Indices of trials to keep','cat','Data Selection'), ...
    arg_nogui({'equalizetrials','BalanceTrials'},false,[],'Equalize the number of trials between two conditions','cat','Data Selection') ...
    );

% commit EEG variable to workspace
[data g] = hlp_splitstruct(g,{'EEG'});
arg_toworkspace(data);
clear data;

if g.verb
    fprintf('\n---------------------------------\n');
    fprintf('Pre-processing condition %s\n', ...
        fastif(isempty(EEG.condition),'[Untitled]',EEG.condition));
end

if g.verb==2
    % reset the color list
    hlp_getNextUniqueColor(hsv(10),'reset');
end

% make sure the data loaded
% EEG = eeg_checkset(EEG,'loaddata');

%% reconstruct ica activations if needed
if strcmpi(g.sigtype.arg_selection,'components') && isempty(EEG.icaact)
    EEG = hlp_icaact(EEG);
end

%% reconstruct times vector, if needed
% this will create EEG.pnts linearly spaced between EEG.xmin and EEG.xmax
% timepoints are converted to ms
if length(EEG.times) ~= EEG.pnts
    if g.verb
        fprintf('The length of EEG.times appears to be inconsistent with the number of time points (EEG.pnts).\nI will reconstruct the times vector using linear increments from EEG.xmin to EEG.xmax, converted to ms.\n');
    end
    EEG.times = linspace(EEG.xmin*1000,EEG.xmax*1000,EEG.pnts); % ms
end

%% create a new CAT structure and copy relevant fields from EEG

% extract any existing configs (so we can restore them later)
if ~g.resetConfigs && isempty(hlp_checkeegset(EEG,'configs'))
    configs = EEG.CAT.configs;
else
    configs = [];
end
switch lower(g.sigtype.arg_selection)
    
    case 'components'
        % specify default component names
        curComps = 1:size(EEG.icaweights,1);
        if isempty(g.varnames)
            MyComponentNames = strtrim(cellstr(num2str(curComps'))');
        else
            MyComponentNames = g.varnames;
        end
        % create dataset
        EEG.CAT = hlp_sift_emptyset(               ...
            'srcdata', EEG.icaact,                 ...
            'nbchan',  size(EEG.icaact,1),         ...
            'curComps',curComps,                   ...
            'curComponentNames', MyComponentNames, ...
            'times',   EEG.times,                  ...
            'pnts',    size(EEG.icaact,2),         ...
            'trials',  size(EEG.icaact,3)          ...
            );
    case 'channels'
        % specify default channel names
        if isempty(g.varnames)
            if isfield(EEG,'chanlocs') && ~isempty(EEG.chanlocs)
                MyChannelNames = {EEG.chanlocs.labels};
            else
                MyChannelNames = strtrim(cellstr(num2str((1:EEG.nbchan)'))');
            end
        else
            MyChannelNames = g.varnames;
        end
        % create dataset
        EEG.CAT = hlp_sift_emptyset(               ...
            'srcdata', EEG.data,                   ...
            'nbchan',  size(EEG.data,1),           ...
            'curComps',1:size(EEG.data,1),         ...
            'curComponentNames', MyChannelNames,   ...
            'times',   EEG.times,                  ...
            'pnts',    EEG.pnts,                   ...
            'trials',  EEG.trials                  ...
            );
    case 'sources'
        curComps = 1:size(EEG.srcpot,1);
        if isempty(g.varnames)
            MyComponentNames = strtrim(cellstr(num2str(curComps'))');
        else
            MyComponentNames = g.varnames;
        end
        % create dataset
        EEG.CAT = hlp_sift_emptyset(               ...
            'srcdata', EEG.srcpot,                 ...
            'nbchan',  size(EEG.srcpot,1),         ...
            'curComps',curComps,                   ...
            'curComponentNames', MyComponentNames, ...
            'times',   EEG.times,                  ...
            'pnts',    size(EEG.srcpot,2),         ...
            'trials',  size(EEG.srcpot,3)          ...
            );
end
% restore any previous configs (if they exist)
if ~isempty(configs)
    fn= fieldnames(configs);
    for k=1:length(fn)
        EEG.CAT.configs.(fn{k}) = configs.(fn{k});
    end
end
% do a little error-checking
if ~isempty(g.varnames) && length(g.varnames) ~= EEG.CAT.nbchan
    error('pre_prepData:badVarNamesLength','The number of entries in the list of %s names (''VariableNames'') must equal the number of %s', ...
        lower(g.sigtype.arg_selection(end:end-1)),lower(g.sigtype.arg_selection));
end

EEG.CAT.signalType = g.sigtype.arg_selection;


%% preprocess the data

% select trials
if ~isempty(g.newtrials)
    if g.verb, fprintf('Selecting trials...\n'); end
    EEG.CAT.srcdata = EEG.srcdata(:,:,g.newtrials);
    EEG.CAT.trials  = length(g.newtrials);
end

% detrend or center data
if g.detrend.arg_selection
    EEG.CAT.srcdata = pre_detrend('data',EEG.CAT.srcdata,'srate',EEG.srate,g.detrend,'verb',g.verb,'arg_direct',true);
end

% differencing
if g.diff.arg_selection
    EEG.CAT.srcdata = pre_diffData('data',EEG.CAT.srcdata,g.diff,'verb',g.verb,'arg_direct',true);
end

% compute band-limited amplitude envelope
% these are added as additional channels
if g.aamp.arg_selection
    EEG = est_aamp('EEG',EEG,g.aamp,'verb',g.verb,'arg_direct',true);
end

% remove bad segments of data
if isfield(g,'badsegments') && ~isempty(g.badsegments)
    for seg=1:size(g.badsegments,1)
        if g.verb, fprintf('Setting interval [%1.2f %1.2f] to NaN\n',g.badsegments(1),g.badsegments(2)); end
        [dummy pnts(1)] = min(abs(EEG.CAT.times-g.badsegments(1)*1000));
        [dummy pnts(2)] = min(abs(EEG.CAT.times-g.badsegments(2)*1000));
        EEG.srcdata(:,pnts(1):pnts(2),:) = NaN;
    end
    if g.verb, fprintf('Done!\n'); end
end

% normalize data
if g.normalize.arg_selection
    EEG.CAT.srcdata = pre_normData('data',EEG.CAT.srcdata,g.normalize,'verb',g.verb,'arg_direct',true);
end

% convert chanlocs to dipfit if desired
if strcmpi(g.sigtype.arg_selection,'channels') ...
        && g.sigtype.chanlocs2dipfit.arg_selection
    EEG.dipfit = hlp_chanlocs2dipfit('chanlocs',EEG.chanlocs,g.sigtype.chanlocs2dipfit);
end

% force double precision
if ~strcmpi(class(EEG.CAT.srcdata),'double')
   EEG.CAT.srcdata = double(EEG.CAT.srcdata);
end



% return an argument list containing subargs for the allowable signal types
function arglist = hlp_getSigTypeArglist(defSigType,defaultNameOnly)

if nargin < 2
    defaultNameOnly = false;
end

arglist = {};

if any(strcmpi(defSigType,'components'))
    arglist{end+1} = {'Components' {}};  %arg({'componentsToKeep','ComponentsToKeep'},[],[],'Vector of component indices to keep')
end

if any(strcmpi(defSigType,'channels'))
    arglist{end+1} = {'Channels' {arg_subtoggle({'chanlocs2dipfit','ConvertChanlocs2Dipfit'},{},@hlp_chanlocs2dipfit,['Create dipfit structure from chanlocs. If channel locations are available and coregistered to an MRI (e.g. the MNI brain), this will construct a "source" dipfit structure which in turn will enable 3D visualization of network structure.' sprintf('\n') 'NOTE: this will overwrite any existing dipfit structure in the current EEG set.'])}};
end

if any(strcmpi(defSigType,'sources'))
    arglist{end+1} = {'Sources' {}};  %arg({'componentsToKeep','ComponentsToKeep'},[],[],'Vector of row indices of sources to keep')
end

if defaultNameOnly
    % return only the name of the default signal type
    arglist = arglist{1}{1};
end
