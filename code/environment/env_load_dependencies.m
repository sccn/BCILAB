function env_load_dependencies(dependency_dir,autocompile)
% Execute dependency loaders in the given directory tree.
% env_load_dependencies(Directory)
%
% This function recognizes 3 types of dependency loader scripts in the dependencies (sub-)folders.
%  1. Scripts called env_exec.m. Each such script is executed by env_load_dependencies.
%
%  2. Scripts called env_add.m. Each such function is executed, and if the 'ans' variable is true 
%     after executing the script (it is set to true before running the script), and no exception 
%     had occured, the containing directory is added using addpath.
%
%  3. Scripts called env_compile.m. The presence of such a function indicates to the loader function
%     that any mex or Java source files in the given directory should be automatically compiled.
%     Additional compiler inputs (such as defines, libraries and include directories), if any, may 
%     be declared as Name=value; statements in the file itself; These will be passed as named arguments
%     to hlp_trycompile (examples can be found throughout the existing dependencies).
%     A special parameter is Skip, which can be set to false on platforms on which the compilation
%     shall be skipped (perhaps because it is known to fail).
%
% In:
%   Directory : directory tree with dependencies to be loaded; by default, searches for a directory
%               named 'dependencies', relative to the current path or the path of this file.
%
%   Compile : whether to act on env_compile.m files (default: true)
%
% Example:
%   % load the dependencies in the C:\mydependencies folder.
%   env_load_dependencies('C:\mydependencies')
%
% Notes:
%   If you intend to compile your toolbox, you cannot use empty .m files. Therefore, your marker files
%   must contain at least a space character to satisfy the compiler.
%
% See also:
%   env_startup
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2011-03-06

if ~exist('dependency_dir','var')
    % other possible directories where the dependencies could be located:
    possible_dirs = {pwd,fileparts(mfilename('fullpath'))};
    while ~exist('dependency_dir','var')
        % for each of the possible directories...
        for d=1:length(possible_dirs)
            % check whether it contains a 'dependencies' directory
            curdir = possible_dirs{d};
            if exist([curdir filesep 'dependencies'],'dir')
                dependency_dir = [curdir filesep 'dependencies'];
                break;
            end
            % otherwise peel off the last directory
            fsp = find(curdir==filesep,1,'last');
            if ~isempty(fsp)
                possible_dirs{d} = curdir(1:fsp-1); end
        end
        if cellfun('isempty',possible_dirs)
            error('Did not find a directory called ''dependencies''. Please specify it.'); end
    end
end

if ~exist('autocompile','var')
    autocompile = true; end

if autocompile
    % check if compilation is suppressed on this platform
    try
        fname = [hlp_homedir filesep '.bcilab' filesep 'bcilab-dont-compile'];
        fid = fopen(fname,'r');
        platforms = fgetl(fid);
        fclose(fid);
        if strfind(platforms,[hlp_hostname ';'])
            autocompile = false; end
    catch
    end
end

load_dependencies(genpath(dependency_dir),autocompile);

% now, make sure that the build-* directories are included after their ancestor directories...
if ~isdeployed
    % get current path order
    added_paths = strsplit(path,pathsep)';
    added_paths = added_paths(cellfun(@(s)strncmp(dependency_dir,s,length(dependency_dir)),added_paths));
    
    super_dirs = {}; % this is a cell array of super directories with build dirs that may have to be moved
    sub_dirs = [];   % this is an array of indices (into added_paths) of sub-directories that need 
                     % to be moved in front of the respective super dirs (indexed the same)
    % we start from the end so that we don't have to change indices
    for p=length(added_paths):-1:1
        try
            cur = added_paths{p};
            k = strfind(cur,[filesep 'build-']);
            if ~isempty(k)
                % if this is a build directory, add to the potential-move list
                super_dirs{end+1} = cur(1:k(end)-1);
                sub_dirs(end+1) = p;
            end
            % check if this is one of the super dirs
            match = strcmp(cur,super_dirs);
            if any(match)
                % move all the respective sub-dirs to right before this position
                indices = sub_dirs(match);
                added_paths = vertcat(added_paths(1:p-1),added_paths(indices),added_paths(setdiff(p:end,indices)));
                super_dirs = super_dirs(~match);
                sub_dirs = sub_dirs(~match);
            end
        catch
            warning('bcilab:dependency_glitch','Glitch during mex/build path reordering.');
            super_dirs = {};
            sub_dirs = [];
        end
    end
    
    % rebuild path
    rmpath(added_paths{:});
    addpath(added_paths{:});
end



% from the given path list, add all paths that contain an active env_add.m
% file (active: either empty or runs without exceptions and ans is not being 
% set to false), and execute all env_exec.m files
function load_dependencies(pathlist,autocompile)
compile_problems = {};
paths = strsplit(pathlist,pathsep);
for px = paths(end:-1:1)
    p = px{1};
    if exist([p filesep 'env_add.m'],'file')
        % env_add file: check if "active"
        [dummy,isactive] = evalc('runscript([p filesep ''env_add.m''])');  %#ok<ASGLU>
        % immediately add all java .jar packages
        dirinfos = dir([p filesep '*.jar']);
        for m={dirinfos.name}
            if isdeployed
                warning off MATLAB:javaclasspath:jarAlreadySpecified; end
            evalc('javaaddpath([p filesep m{1}]);'); 
        end
    else
        isactive = false;
    end
    % env_exec file: run it
    if exist([p filesep 'env_exec.m'],'file')
        runscript([p filesep 'env_exec.m']); end
    % env_compile file: invoke compilation, if necessary
    if autocompile && exist([p filesep 'env_compile.m'],'file')
        settings = get_settings([p filesep 'env_compile.m']);
        % the Skip field determines whether the current platform should be skipped
        if isfield(settings,'Skip')
            docompile = ~settings.Skip;
        else
            docompile = true;
        end
        if ~isequal(settings,false) && docompile
            if ~hlp_trycompile('Directory',p, 'Style','eager', settings);
                compile_problems{end+1} = p; end
        end
    end
    if isactive
        % finally add to MATLAB & Java path
        if ~isdeployed
            addpath(p); end
        % also add java files here (note: this may clear global variables)
        if ~isempty(dir([p filesep '*.class']))
            if isdeployed
                warning off MATLAB:javaclasspath:jarAlreadySpecified; end
            javaaddpath(p); 
        end
    end
end

% handle compilation errors gracefully
if ~isempty(compile_problems) && ~isdeployed
    disp('There were compiling problems for the dependencies in:');
    for k=1:length(compile_problems)
        disp(['  ' compile_problems{k}]); end
    fprintf('\n');
    disp('This means that the affected features of BCILAB will be disabled.');
    if ispc
        disp('These issues are usually resolved by installing a (free) Microsoft Visual Studio Express version that is between five years to a half year older than your MATLAB release.');
        disp('If you install a new compiler, you might have to type "mex -setup" in the MATLAB command line, before MATLAB recognizes it.');
    elseif ismac
        disp('These issues are usually resolved by installing the Xcode application that comes with your OS (with all install options to include gcc compiler support).');
    elseif isunix
        disp('These issues are usually resolved by a installing a version of the gcc compiler package that is between three years to a half year older than your MATLAB release.');
    end
    yn = lower(input('Should BCILAB try to recompile in the future (or skip recompilation if not)? [y/n]','s'));
    if any(strcmp(yn,{'y','n','yes','no'}))
        if yn(1) == 'n'
            if ~exist([hlp_homedir filesep '.bcilab'],'dir')
                if ~mkdir(hlp_homedir,'.bcilab'); 
                    disp('Cannot create directory .bcilab in your home folder.'); end
            end            
            fname = [hlp_homedir filesep '.bcilab' filesep 'bcilab-dont-compile'];
            disp(['BCILAB will add/update a file called ' fname ' which contains the names of computers on which you chose to disable compilation.']);
            fid = fopen(fname,'a');
            if fid ~= -1
                try
                    fwrite(fid,[hlp_hostname ';']);
                    fclose(fid);
                catch
                    try fclose(fid); catch,end
                    disp('There were permission problems while accessing this file. Did not update it.');
                end
            else
                disp('The file name appears to be invalid on your platform; did not update it.');
            end
        end
    else
        disp('Please type ''y'' or ''n'' next time.');
    end
    fprintf('\n');
end


% run the given script (in its own scope...)
% and return the result of the last line
function ans = runscript(filename__)
try
    ans = 1; %#ok<NOANS>
    fprintf('Loading %s...\n',fileparts(filename__));
    run_script(filename__); % Note: This function is equivalent to run(), but works in deployed mode, 
                            %       too. If the code fails here, you probably don't have this 
                            %       function in your path. You may retrieve it (it's in the 'misc' 
                            %       directory of BCILAB), or replace it here by run().
catch e
    disp(['Error running dependency loader ' filename__ '; reason: ' e.message]);
    disp('The dependency will likely not be fully operational.');
    ans = 0; %#ok<NOANS>
end


% run the given settings script (in its own scope...)
% and return the collected defined variables (or false)
function settings = get_settings(filename__)
try
    evalc('run_script(filename__)');
    infos = whos();
    settings = struct();
    for n = {infos.name}
        settings.(n{1}) = eval(n{1}); end
catch
	disp(['The settings script ' filename__ ' contains an error.']);
    settings = false;
end


% Split a string according to some delimiter(s). Not as fast as hlp_split (and doesn't fuse 
% delimiters), but doesn't need bsxfun().
function strings = strsplit(string, splitter)
ix = strfind(string, splitter);
strings = cell(1,numel(ix)+1);
ix = [0 ix numel(string)+1];
for k = 2 : numel(ix)
    strings{k-1} = string(ix(k-1)+1:ix(k)-1); end
strings = strings(~cellfun('isempty',strings));
