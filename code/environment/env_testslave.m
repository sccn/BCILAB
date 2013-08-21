function env_testslave(varargin)
% Run as a test slave: run tests for the toolbox whenever it has changed.
%
% In:
%   ScanInterval : the interval at which the toolbox is scanned for changes, in seconds
%                  (default: 30)
%
%   WaitPeriod : the period that must have passed between the last change and a re-test
%                (default: 120)
%
%   ScanDirectories : the directories and files to scan for changes until re-tests are triggered
%                     (default: {'bcilab:/code', 'bcilab:/*.m'})
%
%   TestDirectory : the directory where the test cases are located
%                   (default: 'bcilab:/resources/testcases')
%
%   IssueTest : function to be called to issue a test run (default: @env_test_bcilab)
%               (assumed to return data related to the test summary)
%
%   PostTest : function to be called after a test run (default: [])
%              (called with the results of IssueTest)
%
%   LogFile : name/location of the logfile, if any
%             (default: 'bcilab:/build/buildslave.log')
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2011-03-17

% read options
o = hlp_varargin2struct(varargin, ...
    {'scan_interval','ScanInterval'},15, ...
    {'wait_period','WaitPeriod'},120, ...
    {'scan_dirs','ScanDirectories'},{'bcilab:/code','bcilab:/*.m'}, ...
    {'issue_test','IssueTest'},@env_test_bcilab, ...
    {'post_test','PostTest'},[], ...
    {'log','LogFile'},sprintf('bcilab:/resources/testcases/log_%s.log',hlp_hostname), ...
    {'report','ReportFile'},'bcilab:/resources/testcases/report_%s-%s.log');

% sanitize inputs
if ischar(o.scan_dirs)
    o.scan_dirs = {o.scan_dirs}; end
for d=1:length(o.scan_dirs)
    o.scan_dirs{d} = env_translatepath(o.scan_dirs{d}); end
o.log = env_translatepath(o.log);
o.report = env_translatepath(o.report);
if ischar(o.post_build)
    o.post_build = str2func(o.post_build); end

% names of marker files
mrk_untested = [proj_dir 'untested'];
mrk_testing = [proj_dir 'testing'];
mrk_tested = [proj_dir 'tested'];

% turn off a few warnings
if ispc
    warning off MATLAB:FILEATTRIB:SyntaxWarning; end
warning off MATLAB:DELETE:FileNotFound;

% run...
quicklog(o.log,'===== Now running as test slave =====');
while 1
    % get the date of the last test (if any)
    testdate = timestamp(mrk_tested);
    
    % get the most recent date of scanned directories
    codedate = max(cellfun(@timestamp,o.scan_dirs));
    
    % check if we are still up-to-date
    if codedate > testdate
        quicklog(o.log,'out of sync');
        rmfile(mrk_tested);
        makefile(mrk_untested);
        rmfile(mrk_testing);
    end
    
    % check if wait period expired ...
    if (codedate - testdate)*24*60*60 > o.wait_period
        % issue a re-test...
        quicklog(o.log,'re-testing...');
        makefile(mrk_testing);
        rmfile(mrk_untested);
        rmfile(mrk_tested);
                
        % issue a re-test
        try
            results = o.issue_test();
        catch e
            quicklog(o.log,evalc('env_handleerror(e)'));
            rethrow(e);
        end
        
        quicklog(o.log,'finished.');
        
        % check whether the code has been edited during the build
        newdate = max(cellfun(@timestamp,o.scan_dirs));
        if (newdate > codedate) && codedate
            quicklog(o.log,'code has been edited concurrently.');
            % if so, revert to untested
            makefile(mrk_untested);
            rmfile(mrk_tested);
            rmfile(mrk_testing);            
        else
            quicklog(o.log,'testing completed.');
            if ~isempty(o.post_test)
                % run the post-testing step
                try
                    quicklog(o.log,'running post-testing step...');
                    o.post_test(results);
                    quicklog(o.log,'post-testing step completed successfully...');
                catch e
                    quicklog(o.log,evalc('env_handleerror(e)'));
                    rethrow(e);
                end
            end            
            % test completed successfully
            makefile(mrk_tested);
            rmfile(mrk_untested);
            rmfile(mrk_testing);
            % write/update report, if applicable
            ds = datestr(clock); ds(ds == ':' | ds == ' ') = '-';
            if ischar(results)
                quicklog(sprintf(o.report,hlp_hostname,ds),results); end
        end
    else
        if codedate > testdate
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


