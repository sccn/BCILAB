function [action,result] = utl_memoize_lookup(exp,memo_ctrl,ctx)
% Check for memoizability and/or availability of the given expression.
% [Action,Result] = utl_memoize_lookup(Key-Expression,Memo-Locations)
% 
% In:
%  Key-Expression : an expression that uniquely identifies the object to be looked up;
%                   this expression is also required to be identical to the one stored in the looked-up object's .tracking.expression field
% 
%  Memo-Locations : cell array of location specifiers, given as name-value pairs. possible names are currently 'disk' and 'memory', and the value
%                   is expected to be an expression that evaluates into 0 or 1, to indicate whether the given location is enabled for 
%                   storage/retrieval or not. The expression may use the symbol @expression to refer to the Key-Expression. The locations 
%                   are looked up in order, so the fastest locations should come first in the list.
%
% Out:
%   Action  : the action to be taken by the caller, one of:
%             'return' : result contains the object that shall be returned as first output argument
%             'memoize': result contains a cell array of store locations to be used for storing the result, when it is eventually computed
%                        - store locations begining with the fileep indicate disk locations, relative to the appropriate memoization domain
%                        - store locations beginning with a . indicate memory locations, relative to the in-memory cache data structure
%             'skip': do not memoize the result
% 
%   Result  : the payload to which the action applies (either the successfully looked up result or the memoization 
%             locations for a subsequent commit)
%
% Notes:
%   If a location further down in the list has the object in question but preceding ones do not, the object is implicitly commited 
%   (copied) to these locations.
%
% See also:
%   utl_memoize_commit
%
%                                                     Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                                     2010-05-23

global tracking;

if ~exist('memo_ctrl','var') || (isempty(memo_ctrl) && ~iscell(memo_ctrl))
    memo_ctrl = hlp_resolve('memoize',[],ctx); end

action = 'skip';
result = {};
% first check whether the current expression is matched by anything in memoize
if iscell(memo_ctrl)
    for i=1:2:length(memo_ctrl)
        % evaluate whether memoization is enabled right now, and get a hash of exp, if so, and [] otherwise
        exp_hash = hlp_microcache('cache_lookup',@hash_if_enabled,memo_ctrl{i+1},exp);
        if exp_hash
            % memoization is enabled for the following location
            location = memo_ctrl{i};
            if strcmp(location,'disk')
                % DISK-BASED MEMOIZATION: compute the storage location according to some path schema
                memo_path = [filesep exp_hash(1:2) filesep exp_hash(3:end) '.mat'];
                if ~isfield(tracking,'cache')
                    % cache disabled...
                    return; end
                % try all known disk paths for lookup
                for p = fieldnames(tracking.cache.disk_paths)'
                    try
                        % try to load the file
                        location = tracking.cache.disk_paths.(p{1});
                        file_location = [location.dir memo_path];
                        t0 = tic;
                        
                        load(file_location,'obj');
                        disp(['loaded ' file_location '.']);
                        
                        % delete the file if it is invalid for some reason (e.g. the machine crashed during saving)
                        if ~exist('obj','var') && strcmp('.mat',file_location(end-3:end))
                            delete(file_location); 
                            error('error');
                        end
                        % check whether the expression matches
                        if ~iscell(obj) || isempty(obj) || ~isfield(obj{1}, 'tracking') || ~isfield(obj{1}.tracking,'expression')
                            warning('BCILAB:exp_beginfun:memo_error','memoized object is not conformat; reverting...');
                            error('error');
                        end
                        if ~utl_same(obj{1}.tracking.expression,exp) && ~isequal_weak(obj{1}.tracking.expression,exp)
                            disp('exp_beginfun: hash conflict between ');
                            disp(['* retrieved: ' exp_fullform(obj{1}.tracking.expression)]);
                            disp(['* requested: ' exp_fullform(exp)]);
                            error('error');
                        end
                        % got the correct result: make the function return it
                        if ~isempty(result) && strcmp(action,'memoize')
                            % ... but first commit it to the other higher-priority memory locations that don't yet have it
                            for r=1:length(result)
                                utl_memoize_commit(obj,result{r}); end
                        end                        

                        result = obj;                       
                        action = 'return';
                        try
                            % ... and record some statistics on the read speed of this location...
                            objinfo = whos('obj');
                            location.readstats(end+1) = struct('size',{objinfo.bytes},'time',{toc(t0)});
                            % write back the updated location record
                            tracking.cache.disk_paths.(p{1}) = location;
                        catch,end
                        return;
                    catch,end
                end
                
                % value not yet stored: we want to memoize after the function's body has been evaluated, so remember the location data
                action = 'memoize';
                result{end+1} = memo_path;
            elseif strcmp(location,'memory')
                % MEMORY-BASED MEMOIZATION: compute the storage location according to a fieldname schema
                memo_field = ['x' exp_hash];
                % check if the requested object is present
                try
                    % found: return
                    obj = tracking.cache.data.(memo_field);
                    % check whether the expression matches
                    if ~iscell(obj) || isempty(obj) || ~isfield(obj{1}, 'tracking') || ~isfield(obj{1}.tracking,'expression')
                        warning('BCILAB:exp_beginfun:memo_error','memoized object is not conformat; reverting...');
                        error('error');
                    end
                    if ~utl_same(obj{1}.tracking.expression,exp) && ~isequal_weak(obj{1}.tracking.expression,exp)
                        disp('exp_beginfun: hash conflict during lookup.');
                        error('error');
                    end
                    tracking.cache.times.(memo_field) = cputime; % for the least-recently used policy
                    % got the correct result: make the function return it
                    if ~isempty(result) && strcmp(action,'memoize')
                        % ... but first commit it to the other higher-priority memory locations that don't yet have it
                        for r=1:length(result)
                            utl_memoize_commit(obj,result{r}); end
                    end
                    result = obj;
                    action = 'return';
                    return;
                catch,end
                % not found: remember to memoize later
                action = 'memoize';
                result{end+1} = ['.' memo_field];
            end
        end
    end
end

% find out whether a given memoization expression for some expression is enabled, and obtain the expression's hash if so...
function hash = hash_if_enabled(memo_exp,exp)
if exp_eval(utl_releasehold(utl_replacerepeated(memo_exp,{exp_rule(@expression,exp)})),inf)
    hash = num2str(145342 + hlp_fingerprint(exp));
else
    hash = [];
end