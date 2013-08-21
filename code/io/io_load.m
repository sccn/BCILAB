function res = io_load(varargin)
% Like MATLAB's builtin load(), except that the filename may be platform-independent.
% Contents = io_load(Filename,Arguments...)
%
% In:
%   Filename    :   Platform-independent file name.
%
%   Arguments...:   optional arguments to load, e.g. the list of variables, or load options
% 
% Out:
%   Contents    :   structure of contents of the file, if desired (if not specified, 
%                   variables are loaded into the workspace)
%
% Example:
%   % load a file
%   io_load('data:/myfile.mat');
%
%   % load only variable 'lastmodel' from the given file
%   io_load('bcilab:/resources/workspaces/myworkspace.mat', 'lastmodel')
%
% See also:
%   load, io_save, env_translatepath
%                                                     
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2010-06-21

if ~isempty(varargin{1}) && varargin{1}(1) == '-' && ~any(varargin{1}=='.')
    contents = load(varargin{1},env_translatepath(varargin{2}),varargin{3:end});
else
    contents = load(env_translatepath(varargin{1}),varargin{2:end});
end

if nargout > 0
    res = contents;
else
    for fn = fieldnames(contents)'
        assignin('caller',fn{1},contents.(fn{1})); end    
end