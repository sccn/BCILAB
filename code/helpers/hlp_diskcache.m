function varargout = hlp_diskcache(options, f, varargin)
% Cache results of function invocations.
% Results... = hlp_diskcache(Settings, Function, Arguments...)
%
% This function maintains a disk cache of function results in a user-specified folder and 
% if a result had already been computed before, it will be immediately looked up instead of being
% computed again.
%
% This function should only be used for computations that take a long enough time to justify the 
% overhead since disk I/O can be relatively slow, since results can fill up the disk quickly, and 
% since there are important safety considerations (see below).
%
% In:
%   Settings : settings that determine where and what to cache (either a cell array of name-value 
%              pairs or a struct). The most important options are:
%              * folder: the folder relative to which the results are saved (default: '.')
%              * freemem: the amount of memory that should remain free on the disk (in GB if > 1, 
%                otherwise a fraction of total space) (default: 80)
%              * bypass: bypass caching system (default: false)
%              See Advanced Options below for further options, and see Pre-defining and Recalling 
%              Settings for a way to separate the settings declaration from the call-site that 
%              invokes hlp_diskcache.
%
%   Function  : function handle to compute a result from some arguments
%
%   Arguments... : arguments to pass to the function
%
% Out:
%   Results... : return values of the function for the given arguments
%
%
% Safety Notes:
%   1) Only *referentially transparent* functions can be cached safely; if a function can give different 
%      outputs for the same arguments, this can lead to subtle bugs; examples include functions that
%      refer to global state (global variables, files on disk). This can be fixed by turning all
%      dependencies of the function into arguments. Also, if a function's desired behavior includes
%      side effects (e.g., creating or updating files), it cannot be cached safely. This can be
%      worked around by caching only the core function that performs the actual computation (if
%      any). When applying functions to "smart" objects, make sure that the call does not
%      inadvertently fall under these categories (e.g., reads, creates or updates files or global
%      variables).
%
%   2) If the execution of a function changes, an obsolete result might be looked up from the cache.
%      Therefore the code of the passed-in function is checksummed against the code that calculated
%      the original result, however *none* of the functions called by that function will be checksummed
%      since it cannot readily be made fast enough. Therefore, if a dependent function changes and
%      that change affects results in the cache, you will want to delete or disable the cache. This is 
%      especially important during debugging sessions.
%
%
% Usability Notes:
%   1) If you are debugging a function whose results are being cached, *bypass* the cache temporarily
%      (you can do this either below in the first code line, or at the call site by passing in 'bypass').
%
%   2) If you do not want the cache to be invalidated when you change a function and you know what
%      you are doing, you can add a version line to your function that you increment whenever you
%      make a change that renders previous results obsolete. This requires very serious programmer 
%      discipline and cannot possibly be enforced when multiple people edit the code at random.
%
%   3) Make sure that your disk is fast enough for the caching to make sense; ideally you want an SSD.
%      Do not use hlp_diskcache for small jobs that are very fast to compute -- use hlp_microcache 
%      instead, which is in-memory (not persistent across MATLAB restarts) and extremely lightweight.
%
%
% Advanced Options:
%      The following further settings can be passed for Settings:
%      * maxsize: don't save result files that are larger than this in bytes (default: Inf)
%      * minsize: don't save result files that are smaller than this in bytes (default: 0)
%      * mintime: don't save results that took less than this time to compute, in seconds (default: 0)
%      * versiontag: syntax of the optional version tag in your function's code (default: '$(funcname)_version<\S+>')
%                    code versioning can be disabled by setting the versiontag to false (discouraged)
%      * subdir: the cache sub-directory relative to the folder; created if missing (default: 'diskcache')
%      * exactmatch_cutoff: if the input is larger than this many bytes, it will not be stored with the result
%                           and will not be compared byte-for-byte during lookup, in bytes (default: 1000000)
%                           note that the hash is usually strong enough to make a byte-for-byte check unnecessary
%      * spot_hashing: if true, faster hashing will be performed on a subset of the data for speed
%                      (default: false)
%      * cleanup: clean up old cache entries when running out of disk space (default: true)
%      * serialize: use serialization for the result, faster than raw save/load for large files (default: true)
%      * permissions: cell array of file permissions to use for created directories and files, as in fileattrib 
%                     (default: {'+w','a'})
%      * bypass_if_folder_missing: if the given folder does not exist, the cache will be bypassed; 
%                                  otherwise the folder will be created if missing (default: false)
%                                  note: this can be useful to prevent inadvertent littering of directories 
%                                  with cache files when running from a different installation
%      * overwrite_files: overwrite existing files (can create broken files when multiple processes write
%                         to the same files in parallel) (default: false)
%      * load_only: if true and no result is in the cache, then this function will return the string
%                   'hlp_diskcache:notfound' (default: true)
%
%
% Pre-defining and Recalling Settings:
%   If Settings is passed in as a string instead of a cell array of name-value pairs or a struct,
%   it is taken as the name of a settings profile (similar to a cache "domain" in hlp_microcache).
%   Note that the profile name should be a valid MATLAB field name.
%   
%   Pre-defining settings for a profile:
%       To assign settings for a named cache profile, call:
%       > hlp_diskcache('myprofile','name',value,'name',value, ...)
%       Where myprofile is the name of the profile for which settings shall be assigned
%       and the names/values are the settings to assign to it (what would normally be passed as a 
%       cell array). Note that settings are not persistent across MATLAB runs, so you'd need to 
%       put them into an initializer or startup function. By default a "clear all" does not clear the 
%       settings (this can be disabled by setting clearable_settings to true in the code below).
%
%   Recalling settings from a profile:
%       To recall settings from a profile in a call of hlp_diskcache, just use the name of the
%       profile instead of the Settings cell array. It is permitted to recall from a profile that
%       has not been defined before (in this case, all defaults will be assumed).
%   
%
% Examples:
%   % if this line is executed for the first time, it is as slow as magic(2000)
%   % the result will be written into a sub-directory of ~/myresults
%   m = hlp_diskcache({'folder','./myresults'},@magic,2000);
%
%   % if it is executed a second time, it is likely much faster than m=magic(2000)
%   % if the result is found on disk
%   m = hlp_diskcache({'folder','./myresults'},@magic,2000);
%
%   % it is also possible to assign settings separately for a named 'profile', and then later recall 
%   % them:
%   hlp_diskcache('myprofile','folder','./myresults','freemem',10);
%   m = hlp_diskcache('myprofile',@magic,2000);
%
%
% See also:
%  hlp_microcache
%
% Depends on:
%  hlp_serialize, hlp_deserialize, hlp_cryptohash
%
%
%                                 Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                 2013-04-15

% Copyright (C) Christian Kothe, SCCN, 2013, christian@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License as published by the Free Software Foundation; either version 2 of the
% License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
% even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License along with this program; if not,
% write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
% USA
dp;

bypass = false;             % whether to bypass caching (switch for debugging)
clearable_settings = false; % whether settings should be clearable by "clear all"

archive_version = 1.0;      % version of the archive format
persistent settings;
persistent have_translatepath;
if isempty(have_translatepath)
    have_translatepath = exist('env_translatepath','file'); end

% parse options
if isstruct(options) || (iscell(options) && mod(length(options),2) == 0 && iscellstr(options(1:2:end)))
    % options are directly given as cell array or struct
    options = assign_defaults(options);
elseif ischar(options) || iscell(options) && mod(length(options),2) == 1 && iscellstr(options([1 2:2:end]))
    % options is referring to the name of a profile (possibly followed by option overrides)
    if ischar(options)
        profilename = options;
        overrides = {};
    else
        profilename = options{1};
        overrides = options(2:end);
    end
    if ~isvarname(profilename)
        error('The given profile name is not a valid MATLAB variable name: %s',profilename); end
    % make sure that it exists
    if ~isfield(settings,profilename)
        settings.(profilename) = assign_defaults({}); end
    if isa(f,'function_handle')
        % recall options from it
        options = settings.(profilename);
        if ~isempty(overrides)
            % have more options
            for o=1:2:length(overrides)
                options.(overrides{o}) = overrides{o+1}; end
        end
    elseif ischar(f)
        % assign options to it
        varargin = [{f} varargin];
        for k=1:2:length(varargin)
            settings.(profilename).(varargin{k}) = varargin{k+1}; end
        if ~clearable_settings
            mlock; end
        return;
    else
        error('Unrecognized syntax: the second argument must be either a function handle or a string.');
    end
else
    error('Unrecognized syntax: the Settings argument must be either a struct, cell array, or string/');
end 

% make path platform-specific
if ~have_translatepath
    options.folder = strrep(strrep(options.folder,'\',filesep),'/',filesep);
else
    options.folder = env_translatepath(options.folder); 
end

% optionally bypass the caching
if bypass || options.bypass || (options.bypass_if_folder_missing && ~exist(options.folder,'dir'))
    [varargout{1:nargout}] = f(varargin{:});
    return;
end

uid_struct = trim_expression({nargout,f,func_version(f,options.versiontag),varargin});
if ~options.spot_hashing
    % get a unique identifier of the computation
    uid = hlp_serialize(uid_struct);
    % get a short hash value of that
    hash = hlp_cryptohash(uid);
else
    % get the UID via fingerprinting
    uid = hlp_fingerprint(uid_struct,false);
    % get a short hash value of that
    hash = hlp_cryptohash(uid);
end
    
% get the file name under which this would be cached
cachedir = [options.folder filesep options.subdir];
filename = [cachedir filesep hash(1:2) filesep hash(3:end) '.mat'];

if exist(filename,'file')    
    try
        % try to load from disk
        result = load(filename);
        % do some sanity checks
        if floor(result.settings.archive_version) > floor(archive_version)
            disp('Note: The file was saved with a newer major archive version. Performing a safe fallback...');
        elseif ~all(isfield(result,{'settings','hash','uid','varargout'}))
            disp('The cached file is apparently malformed (missing some required fields). Performing a safe fallback...'); 
        elseif ~isequal(result.hash,hash)
            disp('Note: the hash does not match that of the file on disk. Performing a safe fallback...');
        elseif ~isempty(result.uid) && ~isequal(result.uid,uid)
            disp('Note: two results yielded the same MD5 hash; this should be a very rare event. Performing a safe fallback...');
        else
            % deserialize data if necessary
            if result.settings.is_serialized
                result.varargout = hlp_deserialize(result.varargout); end
            % all went well: return result
            varargout = result.varargout;
            return;
        end
    catch e
        error_message('Could not look up result from disk',e);
        disp('Performing a safe fallback...');
    end
elseif options.load_only
    fprintf('hlp_diskcache: result not found.\n');
    varargout = {'hlp_diskcache:notfound'};
    return;
end

result.settings = options;
if numel(uid) <= result.settings.exactmatch_cutoff
    result.uid = uid;
else
    result.uid = [];
    clear uid;
end
result.settings.is_serialized = result.settings.serialize;
result.settings.archive_version = archive_version;
result.hash = hash;

% (re)calculate the result
start_time = tic;
[varargout{1:nargout}] = f(varargin{:});
computation_time = toc(start_time);

% prepare result for writeback
if result.settings.is_serialized
    try
        result.varargout = hlp_serialize(varargout);
    catch e
        error_message('Could not serialize the result (try to disable the option serialize)',e);
        disp('Saving result in unserialized form (can be slow).')
        result.varargout = varargout;
    end
else
    result.varargout = varargout;
end

% do some size & time checks
stats = whos('result');
resultsize = stats.bytes;
if resultsize > options.maxsize
    % result too big to store
    return;
elseif resultsize < options.minsize
    % result too small to store
    return;
elseif computation_time < options.mintime
    % computation too short to be worth it
    return;
else
    
    % ensure that we have enough space
    total_space = disk_total_space(filename);
    if options.freemem > 0 && total_space > 0    
        free_space = disk_free_space(filename);
        if options.freemem < 1
            ensured_space = total_space * options.freemem;
        else
            ensured_space = options.freemem*(2^9);
        end
        if total_space && free_space-resultsize < ensured_space
            if ~options.cleanup
                return; end
            % generate a list of all result files
            allfiles = struct();
            records = dir(cachedir);
            for d = 1:length(records)
                record = records(d);
                if record.isdir && length(record.name) == 2 && all(record.name~='.')
                    % get a list of all MATAB files in this subdir
                    files = dir([cachedir filesep record.name filesep '*.mat']);
                    for f=1:length(files)
                        files(f).path = [cachedir filesep record.name filesep files(f).name]; 
                        allfiles(end+1) = files(f); %#ok<AGROW>
                    end
                    % try to remove dirs that are empty
                    if isempty(files)
                        try rmdir([cachedir filesep record.name]); catch,end; end %#ok<CTCH>
                end
            end
            % delete old files as long as ours doesn't yet fit into memory
            [dummy,newest_to_oldest] = sort([allfiles.datenum],'descend'); %#ok<ASGLU>
            while ~isempty(newest_to_oldest) && disk_free_space(filename) - resultsize < ensured_space
                try delete(allfiles(newest_to_oldest(end)).path); catch,end %#ok<CTCH>
                newest_to_oldest = newest_to_oldest(1:end-1);
            end
        end
    else
        persistent message_shown; %#ok<TLEV>
        if isempty(message_shown)
            disp('Note: cannot determine free disk space on your platform, trying to cache results anyway. This message will only be shown once.');
            message_shown = true;
        end
    end
    
    % save the result
    try
        % ensure that the target directory exists
        make_directories(filename,options.permissions);
        % save result file
        fldnames = fieldnames(result);
        if exist(filename,'file') && ~options.overwrite_files
            return; end
        if resultsize >= (2000*1024*1024)
            save(filename,'-struct','result',fldnames{:},'-v7.3');
        else
            save(filename,'-struct','result',fldnames{:});
        end
        % finalize file permissions
        if ~isempty(options.permissions)
            warning off MATLAB:FILEATTRIB:SyntaxWarning
            try
                fileattrib(filename,options.permissions{:});
            catch e
                error_message('Note: could not set permissions for result file',e);
            end
        end
    catch e
        error_message('Could not cache result on disk',e);
    end
end


function options = assign_defaults(options)
% Assign default settings to an options struct / name-value pair list
if iscell(options)
    options = cell2struct(options(2:2:end),options(1:2:end),2); end
if ~isfield(options,'folder')
    options.folder = '.'; end
if ~isfield(options,'maxsize')
    options.maxsize = Inf; end
if ~isfield(options,'minsize')
    options.minsize = 0; end
if ~isfield(options,'mintime')
    options.mintime = 0; end
if ~isfield(options,'freemem')
    options.freemem = 80; end
if ~isfield(options,'versiontag')
    options.versiontag = '$(funcname)_version<\S+>'; end
if ~isfield(options,'serialize')
    options.serialize = true; end
if ~isfield(options,'subdir')
    options.subdir = 'diskcache'; end
if ~isfield(options,'exactmatch_cutoff')
    options.exactmatch_cutoff = 1000000; end
if ~isfield(options,'spot_hashing')
    options.spot_hashing = false; end
if ~isfield(options,'cleanup')
    options.cleanup = true; end
if ~isfield(options,'permissions')
    options.permissions = {'+w','a'}; end
if ~isfield(options,'bypass')
    options.bypass = false; end
if ~isfield(options,'bypass_if_folder_missing')
    options.bypass_if_folder_missing = false; end
if ~isfield(options,'overwrite_files')
    options.overwrite_files = false; end
if ~isfield(options,'load_only')
    options.load_only = false; end


function x = trim_expression(x)
% Recursively trim partially evaluated parts of a data structure containing expressions.
% In particular, x.tracking.expression is replaced by x.
if isfield(x,'tracking') && isfield(x.tracking,'expression')
    x = trim_expression(x.tracking.expression);
elseif iscell(x)
    x = cellfun(@trim_expression,x,'UniformOutput',false);
elseif isfield(x,{'head','parts'})
    x.parts = trim_expression(x.parts);
end


function v = func_version(func,versiontag)
% Get a version identifier of a MATLAB function; can be any of the following
%  * cell array of version strings of a MATLAB function, if present
%  * MD5 hash of the file if unversioned.
%  * string form of the input if there is no accessible file (e.g., anonymous function),
%    or if the versiontag is passed in as false
try
    if ischar(func)
        filename = which(func);
    else
        filename = getfield(functions(func),'file');
    end
catch
    filename = '';
end
func = char(func);
if isequal(versiontag,false) || strncmp(char(func),'@',1)
    v = char(func);
else
    if ~isempty(filename)
        % open the source file
        f = fopen(filename,'r');
        try
            % read the code
            code = fread(f,Inf,'uint8=>char')';
            % check if it contains the version descriptor tag
            v = regexp(code,strrep(versiontag,'$(funcname)',func),'match');
            % otherwise we just hash the entire code
            if isempty(versiontag) || isempty(v)
                v = hlp_cryptohash(code); end
            fclose(f);
        catch %#ok<CTCH>
            try
                fclose(f);
            catch %#ok<CTCH>
            end
            v = func;
        end
    else
        % otherwise use the string representation as version
        v = func;
    end
end

function res = disk_total_space(path)
% Get the amount of total space on the disk (can be 0 if the check fails).
f = java.io.File(path);
res = f.getTotalSpace();


function res = disk_free_space(path)
% Get the amount of free space on the disk (can be 0 if the check fails).
f = java.io.File(path);
res = f.getFreeSpace();


function error_message(msg,e)
% Display a formatted error message with traceback.
if exist('hlp_handleerror','file')
    disp([msg ': ']);
    hlp_handleerror(e);
else
    disp([msg ': ' e.message]);
end


function res = strsplit(str,delims)
% Split a string according to some delimiter(s).
pos = find(diff([0 ~sum(bsxfun(@eq,str(:)',delims(:)),1) 0]));
res = cell(~isempty(pos),length(pos)/2);
for k=1:length(res)
    res{k} = str(pos(k*2-1):pos(k*2)-1); end


function make_directories(filepath,attribs)
% Create directories recursively for a given file path.
paths = strsplit(filepath,filesep);
% find the base directory where the directory creation begins
if filepath(1) == filesep
    % Unix, absolute
    curpath = filesep; first = 1;
elseif ~isempty(strfind(paths{1},':')) && ispc
    % Windows, absolute
    curpath = [paths{1} filesep]; first = 2;
else
    % relative
    curpath = [pwd filesep]; first = 1;
end
% determine where to stop
if filepath(end) == filesep
    last = length(paths);
else
    last = length(paths)-1;
end
% walk...
for i=first:last
    if ~exist([curpath paths{i}],'dir')
        if ~mkdir(curpath,paths{i})
            error(['Unable to create directory ' filepath]);
        else
            % set attributes
            if ~isempty(attribs)
                warning off MATLAB:FILEATTRIB:SyntaxWarning
                try
                    fileattrib([curpath paths{i}],attribs{:});
                catch e
                    error_message(['Note: could not set permissions for created directory (' filepath ')'],e);
                end
            end
        end
    end
    curpath = [curpath paths{i} filesep]; %#ok<AGROW>
end
