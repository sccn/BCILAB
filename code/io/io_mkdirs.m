function io_mkdirs(filepath,attribs)
% Create directories recursively and optionally set permissions.
% io_mkdirs(Filepath,Attributes)
%
% In:
%   Filepath   : path to the file (or directory) in question.
%   Attributes : optional cell array of attributes to use for the created directories, length 1 or 2
%                attribs contains the second or second and third argument to be passed to the fileattrib() function,
%                to set attributes for the current user or specific user groups, e.g., {'+w','a'}.
%
% Notes:
%   x/y/z is treated as a file path (and only x/y/ is created)
%   x/y/z/ is treated as a directory path
%
%   Absolute platform-specific paths (e.g. /xxx/ or C:\xxx\) should be avoided for script 
%   portability reasons (see also env_translatepath); the preferred way is to use paths relative to the store path, the 
%   data path, or the current working directory: 'data:/xxx/', 'store:/xxx/', or 'xxx/'
%
% Examples:
%   % make sure that a particular directory exists
%   io_mkdirs('bcilab:/resources/defaults/')
%
%   % as before, but make it writeable for all users upon creation, if on UNIX
%   io_mkdirs('bcilab:/resources/defaults/',{'+w','a'})
%
% See also:
%   mkdir, fileattrib
%
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2010-03-29

if ~exist('attribs','var')
    attribs = []; end
if length(attribs) > 2
    attribs = attribs{1:2}; end
if ~isempty(attribs) && ~iscellstr(attribs)
    error('file attributes must be a cell array of strings.'); end

filepath = env_translatepath(filepath);
paths = hlp_split(filepath,filesep);
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

% determine where to last
if filepath(end) == filesep 
    last = length(paths); 
else
    last = length(paths)-1; 
end

for i=first:last    
    if ~exist([curpath paths{i}],'dir')
        if ~mkdir(curpath,paths{i})
            error(['unable to create directory ' filepath]);
        else
            % set attributes
            if ~isempty(attribs)
                warning off MATLAB:FILEATTRIB:SyntaxWarning
                fileattrib([curpath paths{i}],attribs{:}); 
            end
        end
    end
    curpath = [curpath paths{i} filesep];
end
