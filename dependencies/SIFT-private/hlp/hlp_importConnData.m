function EEG = hlp_importConnData(varargin)
% Convert/Import connectivity matrices and (optionally) source locations
% into SIFT data structures and store in an EEGLAB data structure
% This can be used to visualize, in SIFT, connectivity data computed
% using another toolbox, such as SPM, GCCA, or Fieldtrip.

% Author: Tim Mullen, 2012, SCCN/INC, UCSD.
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

g= arg_define(varargin,...
    arg({'connDataStruct','ConnectivityDataStructure'},[],mandatory,'Connectivity data structure. This is a structure where each field ConnDataStruct.(<MethodName>) = [vars_to x vars_from x freqs x times] contains a 2D, 3D, or 4D connectivity matrix for a given MethodName. All matrices must be of equal size and first two dimensions must be equal size (num_vars_to==num_vars_from). If this describes an asymmetric (directed) graph, edges in the graph point from columns (vars_from) to rows (vars_to). If there is no frequency information, the 3rd dimension is a singleton (matrix is still 4-D), if there is no time information then the matrix is 3-dimensional (matlab always squeezes out trailing singleton dimensions). If there is neither frequency nor time information, matrix is 2D.'), ...
    arg({'winCenterTimes','WinCenterTimes'},[],[],'Time-stamps (sec). This is a vector of time-stamps for corresponding time-varying connectivity (winCenterTimes=[], if no time-varying connectivity)','shape','row'), ...
    arg({'erWinCenterTimes','EventRelatedWinCenterTimes'},[],[],'Event-related time-stamps (sec). This is a vector of "event-related" time-stamps for event-related analysis. Time t=0 is an event to which trials have been locked, t>0 is post-event and t<0 is pre-event. If empty (default), value is same as winCenterTimes.','shape','row'), ...
    arg({'freqs','Frequencies'},0,[],'Frequencies (Hz). This is a vector of frequencies for corresponding frequency-varying connectivity (freqs=0, if no frequency information)','shape','row'), ...
    arg({'winlen','SlidingWindowLength'},0,[],'Sliding Window Length (sec). This is the length of the sliding window (in sec) if you used a sliding-window analysis to get time-varying connectivity. Default is 0 (no sliding window used)'), ...
    arg({'variableNames','VariableNames','curComponentNames'},{},[],'Variable IDs. This is a cell array of ''user-friendly'' names for each variable -- e.g. anatomical names for each source, or channel names for each channel','type','cellstr','shape','row'), ...
    arg('MNI_XYZ',[],[],'MNI Coords (X,Y,Z) for each variable. An EEGLAB dipfit structure can also be provided. Otherwise, this is an optional matrix of size [num_vars x 3] where each column contains, respectively, X,Y,Z MNI coordinates for the variables (e.g. channels or sources). Subject-specific coordinates can also be provided, but then a subject-specific head model will need to be provided when rendering a brainmovie. This goes into the .posxyz field of a new EEGLAB dipfit structure','shape','matrix'), ...
    arg('EEG',[],[],'EEGLAB data structure. Optional EEGLAB data struct. SIFT data objects will be stored in here. If ommitted, an empty EEGLAB dataset will be created.') ...
    );

% Input checks and conversions
% -------------------------------------------------------------------------
fnames = fieldnames(g.connDataStruct);
if isempty(g.connDataStruct)
    error('Input argument ''ConnectivityDataStructure'' cannot be empty');
end

% get connectivity matrix dimensions
% and check that all connectivity matrices are of equal size
[numVars numVars numFreqs numTimes] = size(g.connDataStruct.(fnames{1}));
for k=2:length(fnames)
    if        ndims(g.connDataStruct.(fnames{k})) ~= ndims(g.connDataStruct.(fnames{1})) ...
       || any( size(g.connDataStruct.(fnames{k})) ~=  size(g.connDataStruct.(fnames{1})) )
        error('All connectivity matrices in ConnectivityDataStructure must be of equal size');
    end
end

% check time-stamp vectors
if isempty(g.winCenterTimes)
    g.winCenterTimes = 0;
end
if isempty(g.erWinCenterTimes)
    g.erWinCenterTimes = g.winCenterTimes;
end
if length(g.winCenterTimes)~=numTimes
    error('The number of time-stamps (length of ''WinCenterTimes'') does not agree with the corresponding dimension (dim 4) of the connectivity matrices');
end
if length(g.erWinCenterTimes)~=length(g.winCenterTimes)
    error('The number of event-related time-stamps (length of ''EventRelatedWinCenterTimes'') must equal the length of WinCenterTimes.');
end

% check frequency vector
if length(g.freqs)~=numFreqs
    error('The number of frequencies (length of ''Frequencies'') does not agree with the corresponding dimension (dim 3) of the connectivity matrices');
end

% check variable names
variableInds = 1:numVars;
if isempty(g.variableNames)
    g.variableNames = cellstr(num2str(variableInds'))';
elseif length(g.variableNames) ~= length(variableInds)
    error('The length of VariableNames must equal the number of variables');
end

% check MNI coords matrix
if ~isempty(g.MNI_XYZ) && ~isstruct(g.MNI_XYZ)
    if size(g.MNI_XYZ,1)~=numVars || size(g.MNI_XYZ,2) ~=3
        error('''MNI_XYZ'' must be either a dipfit structure or a matrix of size [num_vars x 3]');
    end
end

% Make Conn object
% -------------------------------------------------------------------------
Conn = struct(...
    'winCenterTimes'        , g.winCenterTimes     , ... % vector of time-stamps for corresponding time-varying connectivity (g.winCenterTimes=0, if no time-varying connectivity)
    'erWinCenterTimes'      , g.erWinCenterTimes   , ... % vector of "event-related" time-stamps for event-related analysis. Time t=0 is an event to which trials have been locked, t>0 is post-event and t<0 is pre-event. Otherwise, can be the same as g.winCenterTimes
    'freqs'                 , g.freqs                ... % vector of frequencies (Hz) for corresponding frequency-varying connectivity (g.freqs=0, if no frequency information)
    );
Conn = catstruct(Conn,g.connDataStruct);

% Make (abbreviated) model object
% -------------------------------------------------------------------------
MODEL = struct(...
    'winlen', 0 ...     % this is the only field you need from MODEL. It contains the length of the sliding window (in sec) if you used a sliding-window analysis to get time-varying connectivity. Otherwise, just set to 0
    );

% Make (abbreviated) SIFT data structure.
% -------------------------------------------------------------------------
CAT = hlp_sift_emptyset(...
    'curComps'          ,   variableInds    , ... % numbers for variables (e.g. channels or sources) here
    'curComponentNames' ,   g.variableNames   , ... % 'user-friendly' name for each variable -- e.g. anatomical names
    'nbchan'            ,   numVars         , ... % number of variables (e.g. numVars above)
    'Conn'              ,   Conn            , ... % Connectivity structure (see above)
    'MODEL'             ,   MODEL             ... % Model structure (see above)
    );


% Make EEGLAB data structure
% -------------------------------------------------------------------------
if isempty(g.EEG)
    % if g.EEG dataset not provided, create an empty one
    g.EEG = eeg_emptyset;
end
g.EEG.CAT = CAT; % add the the SIFT ('CAT') structure

% Make dipfit structure with source/channel coordinates
% -------------------------------------------------------------------------

% These coordinates are typically (although not necessarily) in MNI space.
% If you are working with channels and have a locations file, you can also
% import the channel locations using EEGLAB routines (c.f. pop_chanedit or
% pop_readlocs). You can then convert the channel locations to a dipfit
% structure using SIFT routine hlp_chanlocs2dipfit(). This will allow you
% to render a 'channel-level' BrainMovie
%
% Below is an example of a super-basic dipfit structure which contains MNI
% source coordinates [X Y Z] copied from X, Y, Z coordinate vectors
%
% See http://sccn.ucsd.edu/wiki/A08:_DIPFIT#DIPFIT_structure_and_functions
% for details on the fields in dipfit structure.
%
% NOTE: SIFT by default uses the freesurfer segmentations of the standard
%       (Colins) MNI brain for surface visualization. If you use your own
%       head model, you will need to supply this surface file via the
%       'Layers' option of vis_causalBrainMovie3D().
%       Also, if using default MNI brain and the network looks rotated
%       90 deg, try changing dipfit.coordformat to 'mni'

if ~isempty(g.MNI_XYZ)
    
    if isstruct(g.MNI_XYZ)
        % just copy the provided dipfit struct
        g.EEG.dipfit = g.MNI_XYZ;
    else
        % build a new dipfit struct
        dipfit.hdmfile  = '';
        dipfit.mrifile  = fullfile(fileparts(utl_whichfile('eeglab')),'plugins','dipfit2.2','standard_BEM','standard_mri.mat');
        dipfit.chanfile = '';
        dipfit.chansel  = 1:numVars;
        dipfit.coordformat = 'MNI';
        dipfit.coord_transform = [0 0 0 0 0 0 1 1 1];
        
        for i=1:numVars
            dipfit.model(i).posxyz      = g.MNI_XYZ(i,:);
            dipfit.model(i).momxyz      = [0 0 0];
            dipfit.model(i).rv          = 0;
            dipfit.model(i).select      = 1;
            dipfit.model(i).diffmap     = [];
            dipfit.model(i).sourcepot   = [];
            dipfit.model(i).datapot     = [];
        end
        g.EEG.dipfit = dipfit;
    end
end


EEG = g.EEG;



