function result_fields = utl_registered_fields(signal, fieldclasses)
% Get registered fields of the desired classes for given signal.
% function Fields = utl_registered_fields(Signal, FieldClasses)
%
% This function returns the field names in the given signal that are 
% of the desired classes (e.g., 'timeseries'). The result is a combination of a set of hard-coded 
% fields and the names listed in the field .tracking.<classname>_fields. Only fields
% that are present in the signal are returned.
%
% In:
%   Signal : a signal for which time-series field names shall be looked up
%
%   FieldClasses: string or cell-string array of field classes that shall be returned
%
% Out:
%   Fields : Cell array of field names (row vector).
%
% See also:
%   utl_register_field
%
%                                       Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                       2013-08-16

% define a database of cached queries and associated results
% (this approach is more than 10x as fast as using unique/intersect every time)
persistent db;

if ischar(fieldclasses)
    fieldclasses = {fieldclasses}; end

% for each desired field class...
result_fields = {};
for fc = fieldclasses
    fieldclass = fc{1};
    fieldclass_fields = [fieldclass '_fields'];
    try
        % generate the query (assuming that the signal usually has the right fields)
        try
            query = [fieldnames(signal)' {'|'} signal.tracking.(fieldclass_fields)];
        catch
            query = fieldnames(signal)';
        end

        % turn into a string
        query = [query{:}];

        try
            % look up from the database (assuming that it is already contained)
            fields = db.(fieldclass).values{strcmp(db.(fieldclass).keys,db.(fieldclass).query)};
        catch
            % not yet in DB: actually determine the timeseries fields
            if isfield(signal,'tracking') && isfield(signal.tracking,fieldclass_fields)
                fields = get_registered_fields(fieldnames(signal),signal.tracking.(fieldclass_fields),fieldclass);
            else
                fields = get_registered_fields(fieldnames(signal),[],fieldclass);
            end
            % init DB record if necessary
            if ~isfield(db, fieldclass)
                db.(fieldclass) = struct('keys',{{}},'values',{{}}); end
            % append to DB
            db.(fieldclass).keys{end+1} = query;
            db.(fieldclass).values{end+1} = fields;
        end
        result_fields = [result_fields fields]; %#ok<AGROW>
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
end

 
function fields = get_registered_fields(field_names,registered_fields,fieldclass)
% This function performs the actual computation

% generate initial list of candidates
switch fieldclass
    case 'timeseries'
        candidates = {'data','icaact','srcpot'};
    case 'timeaxis'
        candidates = {'times','stamps'};
    case 'samplingrate'
        candidates = {'srate'};
    otherwise
        candidates = {};
end

% append whatever is registered in the signal's .tracking.timeseries_fields
candidates = unique([candidates(:); registered_fields(:)]);

% restrict to those that are actually present in the signal
fields = intersect(candidates(:),field_names);
fields = fields(:)';
