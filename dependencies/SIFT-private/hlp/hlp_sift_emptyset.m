function CAT = hlp_sift_emptyset(varargin)
% create a new (empty) SIFT datastructure with default fields initialized
% varargin is an optional set of <name,value> pairs indicating fields and 
% associated values to store into the set
% 
% The SIFT datastructure is organized as follows:
%
% CAT
%   .signalType         Denotes the type of signal (channels or sources)
%   .srcdata            Data matrix. [num_vars x time x trials]
%   .nbchan             Number of variables (rows of srcdata)
%   .pnts               Number of time points (columns of srcdata)
%   .trials             Number of trials (pages of srcdata)
%   .times              Vector of timestamps (sec) for samples of srcdata
%   .curComps           Numeric indices of channels/sources in orig dataset
%   .curComponentNames  Names of components or channels
%   .MODEL              MODEL structure containing VAR model (see hlp_sift_emptymodel)
%   .Conn               Conn structure containing connectivity (see hlp_sift_emptyconn)
%   .configs            Configuration structures for previously applied
%                       processing steps via pop_* functions (this is a sort of history, 
%                       which can translated via hlp_cfg2comstr)
%
% See Also: hlp_sift_emptyconn(), hlp_sift_emptymodel()

if nargin > 1 && mod(length(varargin),2)
    error('SIFT:hlp_sift_emptyset', ...
        ['Additional arguments must be a list of <name,value> pairs. ' ...
         'At least one name is missing its corresponding value.']);
end


CAT.signalType  = '';
CAT.srcdata     = [];
CAT.nbchan      = 0;
CAT.pnts        = 0;
CAT.trials      = 0;
CAT.times       = [];

CAT.curComps    = [];
CAT.curComponentNames = {};

CAT.MODEL       = [];
CAT.Conn        = [];

% get names for config fields
names = hlp_microcache('popfcn',@getPopFcnSuffixes);
% initialize empty config struct for each pop function
for k=1:length(names)
    CAT.configs.(names{k}) = struct([]);
end

if nargin > 1
    CAT = hlp_varargin2struct(varargin,CAT);    
end


function names = getPopFcnSuffixes()
% get the names of all pop functions
% strip off 'pop' prefix and return function name

siftroot = hlp_getSiftRoot();
fpath    = [siftroot filesep 'pop' filesep];
popfns   = wildcardsearch(fpath,'pop_*.m',true,false);
popfns   = regexprep(popfns,'.*pop_','');
names    = regexprep(popfns,'\.m[~]?','');
names    = unique_bc(names);

