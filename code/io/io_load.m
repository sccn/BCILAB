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

fname = varargin{1};
if length(fname)>4 && strcmpi(fname(end-3:end),'.sto')
    % load from .sto file
    if strncmp(fname,'cluster:',8)
        % load over the network
        bytes = par_fileclient(env_translatepath(fname(9:end)),'Verbosity',2);
    else
        % load from disk
        f = fopen(env_translatepath(fname),'r');
        try
            bytes = fread(f,inf,'*uint8');
            fclose(f);
        catch e
            if ~exist(fname,'file')
                error('The file named "%s" does not exist.',fname);
            else
                rethrow(e);
            end
        end
    end
    contents = hlp_deserialize(bytes);
else
    if ~isempty(varargin{1}) && varargin{1}(1) == '-' && ~any(varargin{1}=='.')
        contents = load(varargin{1},env_translatepath(varargin{2}),varargin{3:end});
    else
        contents = load(env_translatepath(varargin{1}),varargin{2:end});
    end
end
    
if isfield(contents,'is_serialized__')
    contents = rmfield(contents,'is_serialized__');
    for fn = fieldnames(contents)'
        contents.(fn{1}) = hlp_deserialize(contents.(fn{1})); end
end

if nargout > 0
    res = contents;
else
    for fn = fieldnames(contents)'
        assignin('caller',fn{1},contents.(fn{1})); end    
end