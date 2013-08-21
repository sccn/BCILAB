function io_loadworkspace(filename)
% Load workspace from a file.
% io_loadworkspace(filename)
%
% Brings up an appropriate file chooser dialog, if necessary.
%
% In:
%   Filename : workspace file name
%
% See also:
%   io_saveworkspace
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-11-19

% set default
if ~exist('filename','var') || isempty(filename)
    filename = 'bcilab:/resources/workspaces'; end

% translate to platform-specific path
filename = env_translatepath(filename);

if exist(filename,'dir')
    % load from given directory
    [filename,filepath] = uigetfile({'*.mat','Workspace files (*.mat)'},'Load workspace',filename);
    filename = [filepath filename];
elseif ~exist(filename,'file')
    % open file chooser from current path
    [filename,filepath] = uigetfile({'*.mat','Workspace files (*.mat)'},'Load workspace');
    filename = [filepath filename];
end

if isempty(filename) || ~ischar(filename)
    % user pressed cancel: nothing to do
elseif exist(filename,'file')
    fprintf('Loading workspace...');
    evalin('base',sprintf('clear; io_load(''%s'');',filename));
    fprintf('done.\n');
else
    error(['Could not find ' filename]);
end
