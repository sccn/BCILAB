function io_saveworkspace(filename,bugreport)
% Save workspace to a file.
% io_saveworkspace(filename)
%
% In:
%   Filename : workspace file name
%
%   Bugreport : whether to save a full bug report (default: false)
%
% See also:
%   io_loadworkspace
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-11-19

global tracking;

% set default
if ~exist('filename','var') || isempty(filename)
    filename = 'bcilab:/resources/workspaces/'; end
if ~ischar(filename)
    error('The given file name must be a string.'); end
if size(filename,1) ~= 1
    error('The given file name must be a non-empty row vector of characters.'); end
if ~exist('bugreport','var') || isempty(bugreport)
    bugreport = false; end

% translate to platform-specific path
filename = env_translatepath(filename);

if exist(filename,'dir')
    % save to given directory
    [filename,filepath] = uiputfile({'*.mat','Workspace files (*.mat)'},'Save workspace',filename);
    filename = [filepath filename];
end

if isempty(filename) || ~ischar(filename)
    % user pressed cancel: nothing to do
else
    % read out specs
    sysspec.computer = computer;
    sysspec.matversion = evalc('ver');
    try sysspec.gccversion = hlp_gcc_version; catch,end
    try sysspec.bciversion = env_version; catch,end
    assignin('base','sysspec',sysspec);
    consolelog = {};
    if bugreport
        if isfield(tracking,'logfile')
            try
                diary off;
                f = fopen(tracking.logfile,'rt');
                while 1
                    line = fgetl(f);
                    if ~ischar(line), break, end
                    consolelog{end+1} = line;
                end
                fclose(f);
                assignin('base','consolelog',consolelog);
            catch
                try diary(tracking.logfile); catch,end
                try fclose(f); catch,end
            end
        end
    end
    fprintf('Saving workspace...');
    evalin('base',sprintf('io_save(''%s'');',filename));
    fprintf('done.\n');
end
