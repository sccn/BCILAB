function env_buildslave(varargin)
% Run as a build slave: recompile the toolbox whenever it has changed.
%
% In:
%   ScanInterval : the interval at which the toolbox is scanned for changes, in seconds
%                  (default: 30)
%
%   WaitPeriod : the period that must have passed between the last change and a rebuild
%                (default: 120)
%
%   ScanDirectories : the directories and files to scan for changes
%                     (default: {'bcilab:/code', 'bcilab:/*.m'})
%
%   ProjectName : name/location of the project to monitor
%                 (default: 'bcilab:/build')
%
%   SpinInterval : interval at which the function is querying whether the build has
%                  finished (default: 5)
%
%   LogFile : name/location of the logfile, if any
%             (default: 'bcilab:/build/buildslave.log')
%
%   IssueBuild : function to be called to issue a new build (default: 'env_compile_bcilab')
%
%   PostBuild : function to be called after a successful build (default: [])
%               (called with the name of the binary)
%
%   WaitPreCompile : time to wait until we expect that the compilation has ramped up
%                    (default: 30)
%
%   WaitPostCompile : time to wait until we expect that the compilation has completed (after mcc's have finished)
%                     (default: 30)
%
%   WaitPostConflict : time to wait after we have had a failed build (e.g. due to concurrent editing)
%                      (default: 120)
%
% Notes:
%   This function will likely not run on Win64, as it requires a process monitoring command to work.
%
% See also:
%   env_compile_bcilab, deploytool
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2011-03-13

% read options
o = hlp_varargin2struct(varargin, ...
    {'scan_interval','ScanInterval'},15, ...
    {'wait_period','WaitPeriod'},120, ...
    {'scan_dirs','ScanDirectories'},{'bcilab:/code','bcilab:/*.m'}, ...
    {'project_name','ProjectName'},'bcilab:/build', ...
    {'spin_interval','SpinInterval'},5, ...
    {'log','LogFile'},'home:/.bcilab/logs/buildslave.log', ...
    {'issue_build','IssueBuild'},@env_compile_bcilab, ...
    {'post_build','PostBuild'},[], ...
    {'wait_precompile','WaitPreCompile'},30, ...
    {'wait_postcompile','WaitPostCompile'},30, ...
    {'wait_postconflict','WaitPostConflict'},120);

% sanitize inputs
if ischar(o.scan_dirs)
    o.scan_dirs = {o.scan_dirs}; end
for d=1:length(o.scan_dirs)
    o.scan_dirs{d} = env_translatepath(o.scan_dirs{d}); end
o.log = env_translatepath(o.log);
if ischar(o.issue_build)
    o.issue_build = str2func(o.issue_build); end
if ischar(o.post_build)
    o.post_build = str2func(o.post_build); end

% name and location of the project file
o.project_name = env_translatepath(o.project_name);
[proj_dir,proj_tag] = fileparts(o.project_name); proj_dir(end+1) = filesep;


% name and location of the binary
binary_dir = [proj_dir proj_tag filesep 'distrib' filesep];
binary_name = proj_tag;
if ispc
    binary_ext = '.exe';
else
    binary_ext = '';
end
binary_path = [binary_dir binary_name binary_ext];


% names of marker files
mrk_notsynched = [o.project_name filesep 'not_synched'];
mrk_synching = [o.project_name filesep  'synching'];
mrk_synched = [o.project_name filesep  'synched'];


% check for status of the process monitoring tool...
if ~is_process_running('matlab')
    error('Process monitoring does not work on your platform.'); end

% turn off a few warnings
if ispc
    warning off MATLAB:FILEATTRIB:SyntaxWarning; end
warning off MATLAB:DELETE:FileNotFound;

% run...
quicklog(o.log,'===== Now running as build slave =====');
cleaner = onCleanup(@cleanup);
while 1
    % get the date of the binary
    bindate = timestamp(binary_path);
    
    % get the most recent date of scanned directories
    codedate = max(cellfun(@timestamp,o.scan_dirs));
    
    % check if we are still sync'ed
    if codedate > bindate
        quicklog(o.log,'out of sync');
        rmfile(mrk_synched);
        makefile(mrk_notsynched);
        rmfile(mrk_synching);
    end
    
    % check if wait period expired and there is no other compilation running...
    if (codedate - bindate)*24*60*60 > o.wait_period && ~is_process_running('mcc')
        % issue a recompilation...
        quicklog(o.log,'recompiling...');
        makefile(mrk_synching);
        rmfile(mrk_notsynched);
        rmfile(mrk_synched);
        
        % rename the binary...
        try
            ds = datestr(clock); ds(ds == ':' | ds == ' ') = '-';
            if exist(binary_path,'file')
                movefile(binary_path,[binary_dir binary_name '-old-' ds binary_ext]); end
        catch
            quicklog(o.log,'cannot move binary %s',binary_path);
        end
        
        % issue a rebuild
        try
            clear functions;
            o.issue_build();
            % wait until mcc has ramped up
            pause(o.wait_precompile);
        catch e
            quicklog(o.log,evalc('env_handleerror(e)'));
            rethrow(e);
        end
        
        % wait until mcc's are done ...
        while is_process_running('mcc')
            pause(o.spin_interval); end
        
        % wait until the packaging has finished...
        pause(o.wait_postcompile);
        
        quicklog(o.log,'finished.');
        
        % check whether the code has been edited during the build
        newdate = max(cellfun(@timestamp,o.scan_dirs));
        if (newdate > codedate) && codedate
            quicklog(o.log,'code has been edited concurrently.');
            try
                % if so, delete the binary
                delete(binary_path);
            catch
                quicklog(o.log,'Cannot delete failed binary. This is a serious condition; exiting...');
                return;
            end
        end
        
        if exist(binary_path,'file')
            % binary present
            quicklog(o.log,'last compilation successful.');
            if ~isempty(o.post_build)
                % run the postbuild step
                try
                    quicklog(o.log,'running post-build step...');
                    o.post_build(binary_path);
                    quicklog(o.log,'post-build step completed successfully...');
                catch e
                    quicklog(o.log,evalc('env_handleerror(e)'));
                    rethrow(e);
                end
            end
            makefile(mrk_synched);
            rmfile(mrk_notsynched);
            rmfile(mrk_synching);
        else
            % compilation must have been unsucessful...
            quicklog(o.log,'last compilation failed...');
            makefile(mrk_notsynched);
            rmfile(mrk_synched);
            rmfile(mrk_synching);
            % wait for a moment to let the dust settle
            pause(o.wait_postconflict);
        end
    else
        if codedate > bindate
            quicklog(o.log,'waiting for edits to finish...'); end
        pause(o.scan_interval);
    end
end


% calc the most recent time stamp of a given file system location
function ts = timestamp(loc)
if any(loc=='*') || ~exist(loc,'dir')
    % assume that it's a file mask / reference
    infos = dir(loc);
else
    % assume that it's a directory reference
    infos = cellfun(@dir,hlp_split(genpath(loc),pathsep),'UniformOutput',false);
    infos = vertcat(infos{:});
end
if ~isempty(infos)
    % take the maximum time stamp
    if ~isfield(infos,'datenum')
        [infos.datenum] = celldeal(cellfun(@datenum,{infos.date},'UniformOutput',false)); end
    ts = max([infos.datenum]);
else
    % no file: set to the beginning of time...
    ts = 0;
end


% create a status marker file...
function makefile(name)
try
    fid = fopen(name,'w+');
    fclose(fid);
    fileattrib(name,'+w','a');
catch
    if fid ~= -1
        try fclose(fid); catch,end
    end
end


% remove a status marker file...
function rmfile(name)
try
    delete(name);
catch
end


% check if a process is running using OS APIs
function tf = is_process_running(name)
if ispc
    cmd = 'wmic process get description';
else
    cmd = 'ps -A';
end
[errcode,content] = system(cmd);
if errcode
    error('Process monitoring non-operational on your platform. This means you''ll need to compile manually using env_compile_bcilab.'); end
tf = ~isempty(strfind(lower(content),lower(name)));


% clean up any stray MCC tasks
function cleanup(varargin)
if ispc
    system('taskkill /F /IM mcc.exe');
else
    system('killall mcc');
end
