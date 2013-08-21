function meta = utl_parse_metadata(varargin)
% Define a meta-data specification for a dataset.
% Data = utl_parse_metadata(Options...)
%
% This is an internal function that is used to allow for a convenient specification of (potentially 
% very complex) stream meta-data from GUIs and the command line - especially if a file or variable
% is available which already specifies a significant chunk of these data.
%
% In:
%   Options... : name-value pairs denoting meta-data fields; field names include:
%
%                'datasource' : Optionally, a dataset on disk or a MATLAB workspace variable, to
%                               serve as the source of meta-data.
%
%                'srate' : sampling rate; must be given, if not specified in the data source.
%
%                'chanlocs' : channel locations struct array / name cell array (or channel count,
%                             if using DataRiver); must be given, if not specified in the data
%                             source.
%
%                ... and optionally any other field that may be encountered in an EEGLAB set.
%
% Out:
%   Data : A struct with meta-data fields
%
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2010-11-19

% read core options
opts = arg_define(varargin, ...
    arg({'datasource','DataSource'},'',[],'Source variable/file. File name, structure, cell array of name-value pairs, or workspace variable name / expression. If this is left unspecified, sampling rate and channels must be given directly.','type','char','shape','row'), ...
    arg({'srate','SamplingRate'},512,[],'Sampling rate of the input data. Depends on the amplifier & DataRiver settings. If specified ,takes precedence of the sampling rate in the data source.'), ...
    arg({'chanlocs','ChannelLabels'},[],[],'Channel locs/labels. Cell-string array of channel labels in the input stream, usually according to the 10-20 system for EEG. If specified, takes precedence of the channels in the data source. If specified as a number k, taken as the first k DataRiver channels.','type','expression','shape','row'));

% parse the meta-data structure, beginning with the datasource
meta = opts.datasource;
if ischar(meta) && ~isempty(meta)
    % if given as string, first try to evaluate it in the workspace
    try
        meta = evalin('base',meta);
    catch
        % if that fails, try to load it as a file name
        try
            meta = exp_eval(io_loadset(meta)); 
        catch
            error('The given meta-data string could not be interpreted (neither as a file nor a workspace variable).');
        end
    end
elseif iscell(meta)
    % cell arrays are interpreted as name-value pairs
    meta = hlp_varargin2struct(meta);
elseif isempty(meta)
    meta = struct();
elseif isfield(meta,{'head','parts'})
    meta = exp_eval(meta);
elseif ~isstruct(meta)
    error('The given meta-data cannot be interpreted.');
end

% next, add all other fields in in_metadata
meta = hlp_varargin2struct(rmfield(opts,'datasource'),meta);

if ~isfield(meta,'chanlocs') 
    error('chanlocs must be specified.'); end
if ~isfield(meta,'srate') 
    error('srate must be specified.'); end

% remove all meta-data fields that would change during accumulation (i.e. are not static)
meta = rmfield(meta,intersect(fieldnames(meta)',{'data','icaact','events','epoch','xmax','pnts','urevent'}));

% auto-generate default datariver channel labels...
if isnumeric(meta.chanlocs)
    meta.chanlocs = cell(1,meta.chanlocs);
    for c=1:length(meta.chanlocs)
        group = ceil(c/32)-1; % 0-based group index (A1..A32, B1..B32, ...)
        meta.chanlocs{c} = [char('A'+group) num2str(c-group*32)];
    end
end
