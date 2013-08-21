function utl_memoize_commit(obj,id,inputbytes)
% Commit an object to a memory location.
% utl_memoize_commit(Object,Location,InputBytes)
%
% In:
%   Object   : the object to be memoized; required to have a .tracking.expression field, which
%              uniquely identifies the object if a .tracking.computation_time field is present, the
%              disk cache location may vary accordingly (as specified in env_startup)
%             
%   Location : Memory location; expected to originate from a prior call to utl_memoize_lookup, which
%              returns it as a result (attached to the memoize action).
%
%   InputBytes : size of the input to the calculation that produced the object; used to determine
%               whether caching is beneficial or not
%
% See also:
%   utl_memoize_lookup
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-05-23

global tracking;
objinfo = whos('obj');
if id(1) == '.'
    % in-memory memoization
    id = id(2:end);
    % here, we are subject to a quota. delete least-recently used items until there is enough capacity to fit the object
    if ~isfield(tracking,'cache')
        tracking.cache = struct(); end
    if ~isfield(tracking.cache,'sizes')
        tracking.cache.sizes = struct(); end
    if ~isfield(tracking.cache,'times')
        tracking.cache.times = struct(); end
    while 1
        if ~isempty(fieldnames(tracking.cache.sizes))
            cache_size = sum(cell2mat(struct2cell(tracking.cache.sizes)));
        else
            cache_size = 0;
        end
        if cache_size + objinfo.bytes < tracking.cache.capacity
            % the object fits into the cache: store and break
            [tracking.cache.data.(id),tracking.cache.sizes.(id),tracking.cache.times.(id)] = deal(obj,objinfo.bytes,cputime);
            break;
        end
        % object does not fit: try to delete the least-recently used object from cache, if any
        if isempty(fieldnames(tracking.cache.times))
            break; end
        try
            [tmp,I] = sort(cell2mat(struct2cell(tracking.cache.times))); %#ok<ASGLU>
            fn = fieldnames(tracking.cache.times); fn = fn{I(1)};
            [tracking.cache.data,tracking.cache.times,tracking.cache.sizes] = deal(rmfield(tracking.cache.data,fn),rmfield(tracking.cache.times,fn),rmfield(tracking.cache.sizes,fn));
        catch
        end
    end
elseif id(1) == filesep
    % on-disk memoization
    try
        % get the time thresholds for the cache locations & sort them
        [times,I] = sort(cellfun(@(c)c.time,struct2cell(tracking.cache.disk_paths)),'descend');
        tags = fieldnames(tracking.cache.disk_paths);
        tags = tags(I);
        domain_idx  = [];
        if isfield(obj{1}.tracking,'computation_time')
            for i=1:length(times)
                if obj{1}.tracking.computation_time >= times(i)
                    % required computation time threshold for this domain exceeded (i.e. the domain applies):
                    % take the appropriate dir to store results
                    domain_idx = i;
                    break;
                end
            end
            if isempty(domain_idx)
                % no domain applied: skip this location
                return; end
        else
            % take the last (fastest-computation) domain if the computation time is not available for some reason
            domain_idx = length(tags);
        end
        % get the location where we want to store
        location = tracking.cache.disk_paths.(tags{domain_idx});
        
        % test if it makes sense to cache this data set, depending on computation time and disk speeds
        try
            write_speed = regress([location.writestats.time]',[location.writestats.size]');
            read_speed = regress([location.readstats.time]',[location.readstats.size]');
            write_time_estimate = objinfo.bytes * write_speed;
            read_time_estimate = objinfo.bytes * read_speed;
            input_read_time_estimate = inputbytes * read_speed;
            time_without_caching = (input_read_time_estimate + obj{1}.tracking.computation_time) * tracking.cache.reuses;
            time_with_caching = write_time_estimate + read_time_estimate * tracking.cache.reuses;
            % if we don't make it at least 2x faster, we skip the caching
            if time_without_caching < time_with_caching * 2
                return; end                
        catch,end
        
        % check if we exceed the quota and have to delete old data (only supported starting with MATLAB 2007b)
        try 
            free_space = location.space_checker.getFreeSpace;
            % we use this to determine whether disk space can be queried for this location (it's 0 if not)
            total_space = location.space_checker.getTotalSpace; 
            if ~total_space
                disp_once(['Note: the disk space on ' location.dir ' cannot be queried; disabling capacity constraint for the cache.']); end
            if total_space ~= 0 && (free_space - objinfo.bytes < location.free)
                % get an updated view of the files in the cache
                if ~isfield(location,'snapshot')
                    location.snapshot = struct(); end
                % define a function to retain only directories/files whoes name doesn't begin with a '.'
                properize = @(dirlist) dirlist(cellfun(@(d)d(1)~='.',{dirlist.name}));
                % find all cache branches (the two-digit subdirectories); this call is cached by the OS
                curdirs = properize(dir(location.dir));
                % update our records for those branches that we don't have yet or whose timestamps are obsolete                                
                dirnames = cellfun(@(n)['dir' n],{curdirs.name},'UniformOutput',false);
                for d = 1:length(curdirs)
                    fname = dirnames{d};
                    if ~isfield(location.snapshot,fname) || location.snapshot.(fname).datenum ~= curdirs(d).datenum
                        location.snapshot.(fname) = curdirs(d);
                        % get the list of all contained files
                        files = properize(dir([location.dir filesep curdirs(d).name]));
                        % and make them full path names...                        
                        for f=1:length(files)                            
                            files(f).path = [location.dir filesep curdirs(d).name filesep files(f).name]; end
                        location.snapshot.(fname).files = files;
                        % if the directory has no files, delete it...
                        if isempty(files)
                            try rmdir([location.dir filesep curdirs(d).name]); catch,end
                        end
                    end
                end
                location.snapshot = rmfield(location.snapshot,setdiff(fieldnames(location.snapshot),dirnames));                
                % get a flat list of all files
                allfiles = [];
                for fn = dirnames
                    files = location.snapshot.(fn{1}).files;
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
                while location.space_checker.getFreeSpace - objinfo.bytes < location.free
                    % if none, return (then, our current data set does not fit into the cache)
                    if isempty(newest_to_oldest)
                        return; end
                    % delete the oldest one
                    try delete(allfiles(newest_to_oldest(end)).path); catch,end
                    newest_to_oldest = newest_to_oldest(1:end-1);
                end
            end
        catch,end
        
        % create directories
        filepath = [location.dir id];
        io_mkdirs(filepath,{'+w','a'});
        t0 = tic;
        % save
        if objinfo.bytes >= (2000*1024*1024)
            save(filepath,'obj','-v7.3');
        else
            save(filepath,'obj');
        end
        % make writable
        warning off MATLAB:FILEATTRIB:SyntaxWarning
        fileattrib(filepath,'+w','a');
        % record some statistics on the write speed to this location...
        location.writestats(end+1) = struct('size',{objinfo.bytes},'time',{toc(t0)});
        % write back the updated location record
        tracking.cache.disk_paths.(tags{domain_idx}) = location;
    catch, end
else
    error('unknown memoization location specified.');
end
