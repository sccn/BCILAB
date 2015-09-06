% define some placeholder variables (for instructive purposes only)
numVars     = 3;   % 3 variables (e.g. channels or sources)
numFreqs    = 1;   % no frequency information or only 1 freq
numTimes    = 10;  % 10 time-stamps

% Define 3-D or 4-D connectivity matrix. Edges in the corresponding graph 
% point (information flows) from columns to rows. If there is no frequency 
% information, the 3rd dimension is a singleton (matrix is still 4-D), if 
% there is no time information then the matrix is 3-dimensional (matlab 
% always squeezes out trailing singleton dimensions).
ConnMatrix1 = rand(numVars,numVars,numFreqs,numTimes);  
ConnMatrix2 = rand(numVars,numVars,numFreqs,numTimes);

% this is the connectivity structure
Conn = struct(...
    'ConnMeasure1NameHere'  , ConnMatrix1 , ... % [numvars x numvars x numfreqs x numtimes] 3D or 4D connectivity (or other) matrix. Replace field name 'ConnMeasureNameHere' with the name of the connectivity measure as you want it to appear in figure plots (usually a short acronym such as GC or PDC). 
    'ConnMeasure2NameHere'  , ConnMatrix2 , ... % ... additional connectivity measures -- or other measures like spectral density matrix -- can be added as additional fields. All Measure matrices much have the same dimensionality. Augment with zeros if necessary.
    'winCenterTimes'        , 1:numTimes  , ... % vector of time-stamps for corresponding time-varying connectivity (winCenterTimes=0, if no time-varying connectivity)
    'erWinCenterTimes'      , 1:numTimes  , ... % vector of "event-related" time-stamps for event-related analysis. Time t=0 is an event to which trials have been locked, t>0 is post-event and t<0 is pre-event. Otherwise, can be the same as winCenterTimes
    'freqs'                 , 1:numFreqs    ... % vector of frequencies (Hz) for corresponding frequency-varying connectivity (freqs=0, if no frequency information)
    );

% this is the (abbreviated) model structure
MODEL = struct(...
        'winlen', 0 ...     % this is the only field you need from MODEL. It contains the length of the sliding window (in sec) if you used a sliding-window analysis to get time-varying connectivity. Otherwise, just set to 0
        );

% this is the (abbreviated) SIFT data structure. you can use 
% hlp_sift_emptyset() to create a consistent SIFT data structure template
CAT = hlp_sift_emptyset(...
        'curComps'          ,   [1 2 3]                     ,           ... % numbers for variables (e.g. channels or sources) here
        'curComponentNames' ,   {'var1','var2','var3'}      ,           ... % 'user-friendly' name for each variable -- e.g. anatomical names
        'nbchan'            ,   3                           ,           ... % number of variables (e.g. numVars above)
        'Conn'              ,   Conn                        ,           ... % Connectivity structure (see above)
        'MODEL'             ,   MODEL                                   ... % Model structure (see above)
        );


% this is an empty EEGLAB structure
EEG = eeg_emptyset;
EEG.CAT = CAT; % ... containing the SIFT ('CAT') structure

% Now you will have to put your source coordinates into a 'dipfit' structure
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

% some dummy MNI coordinates for 3 variables
MNI_SRC_COORDS_X = [76   82 -46];
MNI_SRC_COORDS_Y = [31  -10   1];
MNI_SRC_COORDS_Z = [-18 -14  26];

dipfit.hdmfile  = '';
dipfit.mrifile  = '';
dipfit.chanfile = '';
dipfit.chansel  = 1:numVars;
dipfit.coordformat = 'Spherical';
dipfit.coord_transform = [0 0 0 0 0 0 1 1 1];

for i=1:numVars
    dipfit.model(i).posxyz      = [MNI_SRC_COORDS_X(i) MNI_SRC_COORDS_Y(i) MNI_SRC_COORDS_Z(i)];
    dipfit.model(i).momxyz      = [0 0 0];
    dipfit.model(i).rv          = 0;
    dipfit.model(i).select      = 1;
    dipfit.model(i).diffmap     = [];
    dipfit.model(i).sourcepot   = [];
    dipfit.model(i).datapot     = [];
end
EEG.dipfit = dipfit;
        



