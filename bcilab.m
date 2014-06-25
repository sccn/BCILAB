function bcilab(varargin)
% BCILAB startup function
%
% The startup function determines startup options according to a config file and then starts the
% BCILAB toolbox based on config options. If a file named bcilab_config.m exists in the user's home
% directory (under .bcilab) it will take precedence, otherwise the file of the same name in the
% bcilab folder will be used. The file name (or path) can also be specified directly as the first
% argument to this function.
%
% In addition, startup options can also be passed directly to this function, in the form of name-
% value pairs. These options selectively override those that are specified in the config file.
%
% Generally, the way to load bcilab is to cd into the BCILAB path, and execute this function.
% The BCILAB path should not be added to the default MATLAB path, especially *never* recursively.
%
% In:
%   ConfigScript : optional file name of a configuration script (assigning values to config 
%                  variables, see bcilab_config.m for an example)
%
%   Options... : optional name-value pairs to override options; for possible options, 
%                see env_startup
%
% Notes:
%   The config files are regular MATLAB script and may include conditional statements (depending on
%   the platform used), etc. However, to allow modification of these files via the GUI, it is
%   recommended to keep them simple.
%
% Examples:
%   % start the toolbox using default settings (found in bcilab_config.m)
%   cd /my/path/to/bcilab-0.9-beta2b; bcilab
%
%   % start the toolbox, passing a particular config script
%   cd /my/path/to/bcilab-0.9-beta2b; bcilab('myconfig.m')
%
%   % start the toolbox, passing a particular config script, and also override a few parameters directly
%   cd /my/path/to/bcilab-0.9-beta2b; bcilab('myconfig.m','data','/data','worker',true)
%
% See also:
%   env_startup
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-10-30

disp('starting BCILAB...');
if mod(nargin,2) == 0
    % an even number of arguments (or none) was passed -- therefore, no config script was given, 
    % so assume the default one    
    configscript = 'bcilab_config';
else
    % an odd number of arguments was passed 
    % --> the first argument is the name of the config script
    if ~ischar(varargin{1})
        error('If an odd number of parameters is provided then the first parameter must be the name of a config script.'); end
    configscript = varargin{1};
    % ... and the rest are the options.
    varargin = varargin(2:end);
end
if ~iscellstr(varargin(1:2:end))
    error('The arguments to bcilab.m must be name-value pairs.'); end

% check if the user's home directory (.bcilab sub-directory) contains a bcilab_config.m, and prefer 
% that if present
home = char(java.lang.System.getProperty('user.home'));
if exist([home filesep '.bcilab' filesep configscript '.m'],'file')
    configscript = [home filesep '.bcilab' filesep configscript '.m']; end

if isdeployed && ~isempty(varargin)
    try
        varargin = cellfun(@eval,varargin,'UniformOutput',false);
    catch
        disp('When passing command-line arguments, make sure that strings are enclosed in '''' characters.');
        disp('Otherwise, all of them will be assumed to be strings.');
    end
end


% sanitize and run the config script
if ~any(configscript == filesep)
    % relative path given
    if isdeployed
        if ~any(configscript == '.')
            configscript = [configscript '.m']; end
        % we assume that the configscript lies in a directory at or above the binary
        tmpdir = [pwd filesep];
        delims = strfind(tmpdir,filesep);
        for k=length(delims):-1:1
            scriptpath = [tmpdir(1:delims(k)) configscript];
            if exist(scriptpath,'file')
                configscript = scriptpath;
                break; 
            end
        end
        if ~exist(scriptpath,'file')
            error('Config script %s not found in any path relative to the binary.',configscript); end
    else
        configscript = which(configscript); 
    end
end
fprintf('running config script %s...\n',configscript);
addpath([fileparts(mfilename('fullpath')) filesep 'code' filesep 'misc']);
run_script(configscript);

% collect the variables as defined in the config script into a struct
infos = whos();
for n = {infos.name}
    settings.(n{1}) = eval(n{1}); end

% load the toolbox (using code/environment/env_startup)
% (using the settings from the config script, possibly overridden by name-value pairs in varargin)
disp('running startup function...');
if ~isdeployed
    cd([fileparts(mfilename('fullpath')) filesep 'code' filesep 'environment']); end

try
    env_startup(settings,varargin{:},'configscript',configscript);
catch
    % store a crash report for debugging in deployed cases
    report.settings = settings;
    report.varargin = varargin;
    report.lasterror = lasterror; %#ok<*LERR>
    save bcilab_crashreport report
    rethrow(lasterror);
end
