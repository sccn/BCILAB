function [action,result] = utl_memoize_lookup(exp,settings,varargin)
% Check for memoizability and/or availability of the given expression.
% [Action,Result] = utl_memoize_lookup(KeyExpression,CandidateLocations,ExecutionContext)
% 
% In:
%  KeyExpression : an expression that uniquely identifies the object to be looked up;
%                  this expression is also required to be identical to the one stored in the
%                  looked-up object's .tracking.expression field
% 
%  CandidateLocations : cell array of location specifiers, given as name-value pairs. possible names
%                       are currently 'disk' and 'memory', and the value is expected to be an
%                       expression that evaluates to 0 or 1, to indicate whether the given location
%                       is enabled for storage/retrieval or not. Most straightforward is to use 0 or
%                       1 for always-on or always-off, but the expression may contain the symbol
%                       @expression that will be substituted by KeyExpression and then evaluated.
%                       The locations are probed in order, so the fastest locations (RAM, SSD, local
%                       drive) should come first in the list.
%
%  ExecutionContext : execution context (with field .stack as returned by dbstack)
%
% Out:
%   Action  : the action to be taken by the caller, one of:
%             'skip': do not memoize the result
%             'return' : result contains the object that shall be returned as first output argument
%             'memoize': result contains a cell array of store locations to be used for storing the
%                        result, after it has eventually been computed 
%                        - store locations begining with the filesep indicate disk locations, 
%                          relative to the appropriate memoization domain 
%                        - store locations beginning with a . indicate memory locations,
%                          relative to the in-memory cache data structure
% 
%   Result  : if action is 'skip', this is {}
%             if action is 'return', this is a cell array of output arguments to return
%             if action is 'memoize', this is a cell array of locations where to store the result
%                                     once it has been computed (using utl_memoize_commit)
%
% Notes:
%   If a location further down in the list has the object in question but preceding ones do not, the
%   object is implicitly commited (copied) to these locations.
%
% See also:
%   utl_memoize_commit
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-05-23
dp;

action = 'skip';
result = {};

% check if cache is present
global tracking;
if ~isfield(tracking,'cache')
    return; end

% optionally look up the settings for the current scope
if nargin < 2 || isequal(settings,[])
    settings = hlp_resolve('memoize',{},varargin{:}); end

% validate settings
if ~iscellstr(settings(1:2:end))
    error('The given cache settings must be a cell array of the form {''location1'',criterion, ''location2'',criterion, ...} but were: %s.',hlp_tostring(settings)); end

for i=1:2:length(settings)
    
    % check whether caching is enabled for this expression/location/scope combo and get a hash of it if so (otherwise false)
    location = settings{i};
    criterion = settings{i+1};
    exp_size_mb = getfield(whos('exp'),'bytes')/2^20;
    if exp_size_mb > 10
        disp_once('Note: your expression is %.1f MB; this can slow down processing.\n',exp_size_mb); end
    data_hash = hlp_microcache('cache_lookup',@hash_or_false,criterion,exp);
    if ~data_hash
        continue; end
        
    % perform actual lookup
    switch location
        case 'memory'            
            % --- look up from memory ---
            % compute the storage location according to a fieldname schema
            hash_field = ['x' data_hash];
            if isfield(tracking.cache,'data') && isfield(tracking.cache.data,hash_field)
                try
                    % get the object and perform sanity checks
                    obj = tracking.cache.data.(hash_field);
                    if ~iscell(obj) || isempty(obj) || ~isfield(obj{1}, 'tracking') || ~isfield(obj{1}.tracking,'expression')
                        error('NOTE: Memory cache record %s contains unsupported data; assumed {expression}, got %s',hlp_tostring(obj)); end
                    if ~utl_same(obj{1}.tracking.expression,exp) && ~isequal_weak(obj{1}.tracking.expression,exp)
                        error('INFO: Got a hash conflict on memory cache record %s:\n* requested: %s\n* retrieved: %s',hash_field,exp_fullform(exp),exp_fullform(obj{1}.tracking.expression)); end
                    
                    % update recently used time and commit the data to any earlier cache locations
                    tracking.cache.times.(hash_field) = cputime;
                    if ~isempty(result) && strcmp(action,'memoize')
                        for r=1:length(result)
                            utl_memoize_commit(obj,result{r}); end
                    end
                    
                    % return the record
                    result = obj;
                    action = 'return';
                    return;
                catch e
                    disp_once(e.message);
                end
            end
            
            % not in memory cache: remember to memoize later at this location
            action = 'memoize';
            result{end+1} = ['.' hash_field]; %#ok<AGROW>

        case 'disk'
            % --- look up from disk ---
            if ~isfield(tracking.cache,'disk_paths')
                continue; end
            hash_path = [filesep data_hash(1:2) filesep data_hash(3:end) '.sto'];
            % for each candidate disk path where the record might be found...
            for p = fieldnames(tracking.cache.disk_paths)'
                try
                    filename = [tracking.cache.disk_paths.(p{1}).dir hash_path];
                    if exist(filename,'file')                            
                        try
                            % try to load the file...
                            t0 = tic; fprintf('retrieving %s...',filename);                            
                            io_load(filename,'obj');
                            fprintf('%.1f seconds.',toc(t0));
                        catch e
                            % failed to load existing file: check if currently being written to
                            fileage = (now - getfield(dir(filename),'datenum'));
                            if fileage < 1/(24*60)
                                % file was written to within less than one minute, assume that it's still incomplete
                                continue; 
                            else
                                % file is corrupted (e.g., simultaneously written to by two MATLAB's); delete it
                                delete(filename);
                                if exist(filename,'file')
                                    error('WARNING: Found corrupted disk cache record (%s) that could not be deleted.',filename);
                                else
                                    error('NOTE: Deleted corrupted disk cache record (%s); load error was: %s',filename,e.message);
                                end
                            end
                        end

                        % validate the loaded data
                        if ~exist('obj','var')
                            error('NOTE: Disk cache record %s does not contain required variable named ''obj''.',filename); end                            
                        if ~iscell(obj) || isempty(obj) || ~isfield(obj{1}, 'tracking') || ~isfield(obj{1}.tracking,'expression')
                            error('NOTE: Disk cache record %s contains unsupported data; assumed {expression}, got %s',filename,hlp_tostring(obj)); end
                        if ~utl_same(obj{1}.tracking.expression,exp) && ~isequal_weak(obj{1}.tracking.expression,exp)
                            error('INFO: Got a hash conflict on disk cache record %s:\n* requested: %s\n* retrieved: %s',filename,exp_fullform(exp),exp_fullform(obj{1}.tracking.expression)); end

                        % update some statistics on the read speed of this location
                        stats = struct('size',{getfield(whos('obj'),'bytes')},'time',{toc(t0)});
                        try
                            tracking.cache.disk_paths.(p{1}).readstats(end+1) = stats;
                        catch
                            tracking.cache.disk_paths.(p{1}).readstats = stats;
                        end
                        
                        % commit the data to any earlier cache locations
                        if ~isempty(result) && strcmp(action,'memoize')                        
                            for r=1:length(result)
                                utl_memoize_commit(obj,result(r),stats.size); end
                        end

                        % return it
                        action = 'return';
                        result = obj;
                        return;
                    end
                catch e
                    disp_once('Error during cache lookup: %s',hlp_handleerror(e));
                end
            end

            % value not in disk cache: remember to memoize later at this location
            action = 'memoize';
            result{end+1} = hash_path; %#ok<AGROW>
        otherwise
            disp_once('WARNING: the given location type (%s) is not supported by the cache (valid names are ''disk'' and ''memory'').',location);
    end
end

% find out whether a given memoization expression for some expression is enabled, and obtain the expression's hash if so...
function hash = hash_or_false(memo_exp,exp)
if exp_eval(utl_releasehold(utl_replacerepeated(memo_exp,{exp_rule(@expression,exp)})),inf)
    hash = num2str(145342 + hlp_fingerprint(exp));
else
    hash = false;
end