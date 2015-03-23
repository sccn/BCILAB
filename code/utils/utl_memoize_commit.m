function utl_memoize_commit(obj,locations,inputbytes) %#ok<INUSD>
% Commit an object to a memory location.
% utl_memoize_commit(Object,StoreLocations,InputBytes)
%
% In:
%   Object   : the object to be memoized; required to have a .tracking.expression field, which
%              uniquely identifies the object if a .tracking.computation_time field is present, the
%              disk cache location may vary accordingly (as specified in env_startup)
%             
%   StoreLocations : Cell array of locations where to store the object; this is the Result value
%                    returned by utl_memoize_lookup, when its returned Action was 'memoize'.
%
%   InputBytes : size of the input to the calculation that produced the object; used to determine
%                whether caching is beneficial or not
%
% See also:
%   utl_memoize_lookup
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-05-23
dp;

show_diagnostics = false;   % whether to display diagnostics on caching decisions
required_speedup = 1.5;     % this is the minimum required speedup by caching below which we don't cache (only for disk cache)

% check inputs
if nargin < 3 
    error('Three input arguments are required.'); end
if ~iscellstr(locations)
    error('The given Locations argument must be a cell array of strings.'); end
objbytes = getfield(whos('obj'),'bytes');

% check if cache present
global tracking;
if ~isfield(tracking,'cache')
    return; end

% for each storage location...
for loc = locations
    id = loc{1};
    try
        switch id(1)
            case '.'
                % --- store in memory cache ---
                capacity = tracking.cache.capacity;
                if objbytes < capacity
                    if ~isfield(tracking.cache,'sizes') || ~isstruct(tracking.cache.sizes)
                        tracking.cache.sizes = struct(); end
                    if ~isfield(tracking.cache,'times') || ~isstruct(tracking.cache.times)
                        tracking.cache.times = struct(); end
                    if ~isfield(tracking.cache,'data') || ~isstruct(tracking.cache.data)
                        tracking.cache.data = struct(); end
                    % free least-recently-used entries to make up space
                    while capacity < sum(cell2mat(struct2cell(tracking.cache.sizes))) + objbytes
                        % delete oldest field in a manner that's safe under concurrent access
                        cache = tracking.cache;
                        fnames = fieldnames(cache.times);
                        oldest_field = fnames{argmin(cell2mat(struct2cell(cache.times)))};
                        [cache.data,cache.times,cache.sizes] = deal(rmfield(cache.data,oldest_field),rmfield(cache.times,oldest_field),rmfield(cache.sizes,oldest_field));
                        tracking.cache = cache;
                    end
                    % insert into cache
                    id = id(2:end);
                    [tracking.cache.data.(id),tracking.cache.sizes.(id),tracking.cache.times.(id)] = deal(obj,objbytes,cputime);
                end
                
            case filesep
                % --- store in disk cache ---
                if ~isfield(tracking.cache,'disk_paths')
                    return; end
                
                % pick the disk path where we want to store based on how long the computation took;
                % we pick the location that has the largest minimum computation time that is still 
                % shorter than our object's required computation time (the assumption is that longer-
                % running computations are more valuable and cache locations that store only long-
                % running computations are more persistent, so we store in the most persistent place
                % applicable to our data.
                
                % get the min compute time thresholds & names for each cache loc
                location_thresholds = cellfun(@(c)c.time,struct2cell(tracking.cache.disk_paths));
                location_names = fieldnames(tracking.cache.disk_paths);
                if isfield(obj{1}.tracking,'computation_time')
                    % object has a compute time: mask out all thresholds that are not applicable...
                    location_thresholds(location_thresholds > obj{1}.tracking.computation_time) = 0;
                    % ... and pick the highest one that survived (if any)
                    if ~any(location_thresholds)
                        return; end
                    location_name = location_names{argmax(location_thresholds)};
                else
                    % object has no compute time: store at the lowest-threshold location
                    location_name = location_names{argmin(location_thresholds)};
                end                                
                
                % get the location where we want to store
                location_info = tracking.cache.disk_paths.(location_name);
                
                % test if it makes sense to cache this data set, depending on computation time and disk speeds
                write_speed = [location_info.writestats.time] / [location_info.writestats.size];
                read_speed = [location_info.readstats.time] / [location_info.readstats.size];
                write_time_estimate = objbytes * write_speed;
                read_time_estimate = objbytes * read_speed;
                input_read_time_estimate = 0; % inputbytes * read_speed; --> we don't have a good estimator of this since we don't know whether the input will be in memory cache or not
                time_without_caching = (input_read_time_estimate + obj{1}.tracking.computation_time) * tracking.cache.reuses;
                time_with_caching = write_time_estimate + read_time_estimate * tracking.cache.reuses;
                caching_worthwhile = time_with_caching < time_without_caching / required_speedup;
                if show_diagnostics
                    fprintf('Cache diagnostics: write_speed=%f; read_speed=%f; write_time_est=%.3f; read_time_est=%.3f; inp_read_time_est=%.3f; comp_time=%.3f\n',write_speed,read_speed,write_time_estimate,read_time_estimate,input_read_time_estimate,obj{1}.tracking.computation_time);
                    fprintf('                   time_w/o_cache=%.3f; time_w_cache=%.3f,docache=%i\n',time_without_caching,time_with_caching,caching_worthwhile);
                end
                % if we don't make it sufficiently faster, we skip the caching
                if ~caching_worthwhile
                    return; end
                
                % remove old files if we're exceeding the capacity
                if location_info.space_checker.getTotalSpace
                    free_space = location_info.space_checker.getFreeSpace;
                    % check if we exceed the quota and have to delete old data
                    if free_space - objbytes < location_info.free
                        % get an updated view of the files in the cache
                        if ~isfield(location_info,'snapshot')
                            location_info.snapshot = struct(); end
                        % define a function to retain only directories/files whoes name doesn't begin with a '.'
                        sanitize = @(dirlist) dirlist(cellfun(@(d)d(1)~='.',{dirlist.name}));
                        % find all cache branches (the two-digit subdirectories); this call is cached by the OS
                        curdirs = sanitize(dir(location_info.dir));
                        % update our records for those branches that we don't have yet or whose timestamps are obsolete
                        dirnames = cellfun(@(n)['dir' n],{curdirs.name},'UniformOutput',false);
                        for d = 1:length(curdirs)
                            fname = dirnames{d};
                            if ~isfield(location_info.snapshot,fname) || location_info.snapshot.(fname).datenum ~= curdirs(d).datenum
                                location_info.snapshot.(fname) = curdirs(d);
                                % get the list of all contained files
                                files = sanitize(dir([location_info.dir filesep curdirs(d).name]));
                                % and make them full path names...
                                for f=1:length(files)
                                    files(f).path = [location_info.dir filesep curdirs(d).name filesep files(f).name]; end
                                location_info.snapshot.(fname).files = files;
                                % if the directory has no files, delete it...
                                if isempty(files)
                                    try rmdir([location_info.dir filesep curdirs(d).name]); catch,end
                                end
                            end
                        end
                        location_info.snapshot = rmfield(location_info.snapshot,setdiff(fieldnames(location_info.snapshot),dirnames));
                        % get a flat list of all files
                        allfiles = [];
                        for fn = dirnames
                            files = location_info.snapshot.(fn{1}).files;
                            if ~isempty(files)
                                if isempty(allfiles)
                                    allfiles = files;
                                else
                                    allfiles(end+1:end+length(files)) = files;
                                end
                            end
                        end
                        % sort them according to date
                        [x,newest_to_oldest] = sort([allfiles.datenum],'descend'); %#ok<ASGLU>
                        % delete old files as long as ours doesn't yet fit into memory
                        while location_info.space_checker.getFreeSpace - objbytes < location_info.free
                            % if none, return (then, our current data set does not fit into the cache)
                            if isempty(newest_to_oldest)
                                return; end
                            % delete the oldest one
                            try delete(allfiles(newest_to_oldest(end)).path); catch,end
                            newest_to_oldest = newest_to_oldest(1:end-1);
                        end
                    end
                else
                    % this happens when getTotalSpace returns 0
                    disp_once(['Note: the disk space on ' location_info.dir ' cannot be queried; disabling capacity constraint for the cache.']);
                end
                                
                % save result to disk
                storepath = [location_info.dir id];
                t0 = tic; fprintf('committing %s...',storepath);
                io_save(storepath,'obj','-serialized','-makedirs','-attributes','''+w'',''a''');
                fprintf('%.1f seconds.',toc(t0));
                
                % and record some statistics on the write speed to this location
                stats = struct('size',{objbytes},'time',{toc(t0)});
                try
                    location_info.writestats(end+1) = stats;
                catch
                    location_info.writestats = stats;
                end
                
                % write back the updated location record
                tracking.cache.disk_paths.(location_name) = location_info;
                
            otherwise
                error('The given type of cache location is unrecognized: %s (expected .fieldname or %sfilename).',id,filesep);
        end
    catch e
        disp_once('WARNING: got an error while writing to the cache: %s',hlp_handleerror(e));
    end
end
