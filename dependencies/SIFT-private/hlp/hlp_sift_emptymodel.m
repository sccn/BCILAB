function Model = hlp_sift_emptymodel(varargin)
% create a new (empty) Model datastructure with default fields initialized
% varargin is an optional set of <name,value> pairs indicating fields and 
% associated values to store into the set
%
% MODEL
%   .morder             Autoregressive model order
%   .AR                 Autoregressive coefficient matrices [A1 A2 ... Ap]
%                       for model order p. This is of dimension [num_vars x num_vars*p]
%   .PE                 Prediction errors (noise covariance matrices) [E0 E1 .. Ep]
%                       This is of dimension [num_vars x num_vars*(p +1)]
%   .RC                 Reflection (parcorr) coefficients [num_vars x num_vars*(p +1)]
%   .mu                 Process mean [num_vars x 1]
%   .th                 Statistics structure (if using ARFIT)
%   .winStartTimes      Timestamps of window start points (sec)
%   .winstep            Sliding window step size (sec)
%   .winlen             Sliding window length (sec)
%   .algorithm          Algorithm used (see hlp_getMVARalgorithms)
%   .modelclass         Class of algorithm ('mvar')
%   .timeelapsed        CPU time required for each time window (sec)
%   .normalize          Data normalizion method used
%   .modelapproach      Approach used ('Segmentation VAR', 'SSM',...)
%   .taperFcn           Taper function used
%
% See Also: hlp_sift_emptyconn(), hlp_sift_emptyset()

if nargin > 1 && mod(length(varargin),2)
    error('SIFT:hlp_sift_emptymodel', ...
        ['Additional arguments must be a list of <name,value> pairs. ' ...
         'At least one name is missing its corresponding value.']);
end

Model.AR = {};
Model.PE = {};
Model.RC = {};
Model.mu = {};
Model.th = {};
Model.winStartTimes = []; %EEG.CAT.times(g.winStartIdx)/1000;
Model.morder        = [];
Model.winstep       = [];
Model.winlen        = [];
Model.algorithm     = '';
Model.modelclass    = '';
Model.timeelapsed   = [];
Model.normalize     = [];
Model.modelapproach = '';
Model.taperFcn      = '';

% append additional arguments
if nargin > 1
    Model = hlp_varargin2struct(varargin,Model);    
end


