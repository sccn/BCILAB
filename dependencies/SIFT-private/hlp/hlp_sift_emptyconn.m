function Conn = hlp_sift_emptyconn(varargin)
% create a new (empty) Conn datastructure with default fields initialized
% varargin is an optional set of <name,value> pairs indicating fields and 
% associated values to store into the set
% 
% Conn
%   .winCenterTimes         Center of time windows (sec) used for model estimation
%   .erWinCenterTimes       Event-related time window centers (event: t=0) 
%   .freqs                  Connectivity frequencies
%   .dims                   Connectivity matrix dimensions. This should be
%                           some ordering of: {'var_to','var_from','freq','time'}
%   <.MethodName>           Named connectivity matrix with dimensions
%                           specified in .dims (see hlp_getValidConnMethods)
%
% See Also: hlp_sift_emptymodel(), hlp_sift_emptyset()

if nargin > 1 && mod(length(varargin),2)
    error('SIFT:hlp_sift_emptyconn', ...
        ['Additional arguments must be a list of <name,value> pairs. ' ...
         'At least one name is missing its corresponding value.']);
end

Conn.winCenterTimes     = [];
Conn.erWinCenterTimes   = [];
Conn.freqs              = [];
Conn.dims               = {'var_to','var_from','freq','time'};

% append additional arguments
if nargin > 1
    Conn = hlp_varargin2struct(varargin,Conn);    
end