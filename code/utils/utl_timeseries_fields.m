function fields = utl_timeseries_fields(signal)
% Get the time-series fields of te given signal.
% function Fields = utl_timeseries_fields(Signal)
%
% This function returns the field names in the given signal that are 
% carrying time-series information (and which therefore should be filtered,
% buffered, etc.). The result is a combination of a set of hard-coded fields
% and the names listed in the field .tracking.timeseries_fields. Only fields
% that are present in the signal are returned.
%
% In:
%   Signal : a signal for which time-series field names shall be looked up
%
% Out:
%   Fields : Cell array of field names (row vector). It is assumed that the second dimension of
%            these fields is the time axis along which buffering and/or filtering should happen; any
%            number of other dimensions may be present for any field returned by this function.
%
% See also:
%   utl_register_field
%
%                                       Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                       2013-08-16

% define a database of cached queries and associated results
% (this approach is more than 10x as fast as using unique/intersect every time)
persistent keys;
persistent values;

try
    % generate the query (assuming that the signal usually has the right fields)
    try
        query = [fieldnames(signal)' {'|'} signal.tracking.timeseries_fields];
    catch
        query = fieldnames(signal)';
    end

    % turn into a string
    query = [query{:}];

    try
        % look up from the database (assuming that it is already contained)
        fields = values{strcmp(keys,query)};
    catch
        % not yet in DB: actually determine the timeseries fields
        if isfield(signal,'tracking') && isfield(signal.tracking,'timeseries_fields')
            fields = get_timeseries_fields(fieldnames(signal),signal.tracking.timeseries_fields);
        else
            fields = get_timeseries_fields(fieldnames(signal),[]);
        end
        % initialize the DB if necessary, otherwise append
        if ~iscell(keys)
            keys = {query};
            values = {fields};
        else
            keys{end+1} = query;
            values{end+1} = fields;
        end
    end
catch e
    % diagnose the error
    if ~isstruct(signal)
        error('The given argument must be a signal struct.'); end
    if ~isscalar(signal)
        error('The given argument must be a 1x1 struct.'); end
    if ~isfield(signal,'tracking')
        if ~any(isfield(signal,{'data','srate'}))
            if all(isfield(signal,{'head','parts'}))
                error('The given input is an unevaluated expression but this function expects a signal data structure. You may have to first evaluate the data using exp_eval.');
            else
                error('The given data structure does not appear to be a valid signal struct.');
            end
        else
            error('The given signal is lacking the required .tracking field.'); 
        end
    else
        rethrow(e);
    end
end
 
function fields = get_timeseries_fields(field_names,ts_fields)
% This function performs the actual computation

% generate initial list of candidates
candidates = {'data','icaact','srcpot','stamps'};

% append whatever is registered in the signal's .tracking.timeseries_fields
candidates = unique([candidates(:); ts_fields(:)]);

% restrict to those that are actually present in the signal
fields = intersect(candidates(:),field_names);
fields = fields(:)';

