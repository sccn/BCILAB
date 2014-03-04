function fname = io_saveset(signal,fname,overwrite) %#ok<INUSL>
% Save a data set to disk.
% Filename = io_loadset(Signal,Filename,Overwrite)
%
%
% In:
%   Signal : EEGLAB dataset structure to save
%   
%   Filename : name of the file; platform-independent path preferred.
%
%   Overwrite : whether to overwrite the data if it already exists (default: false).
%
% Out:
%   Filename : the file name used to save the data
%
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2014-02-05

if nargin < 3
    overwrite = false; end
if ~ischar(fname)
    error('The given file name must be a string.'); end
if size(fname,1) ~= 1
    error('The given file name must be a non-empty row vector of characters.'); end

% append file ending if necessary
if length(fname)<4 || (~strcmp(fname(end-3:end),'.set') && ~strcmp(fname(end-3:end),'.mat'))
    fname = [fname '.set']; end
    
% translate to platform-dependent path
[fp,fn,fx] = fileparts(env_translatepath(fname));

% create directories, if necessary (with appropriate permissions)
if ~exist(fp,'dir')
    io_mkdirs([fp filesep],{'+w','a'}); end

% check for overwrite
if overwrite && exist([fp filesep fn fx],'file')
    return; end

% save
fprintf('Saving set to %s...',[fp filesep fn fx]);
if strcmp(fname(end-3:end),'.set')
    evalc('pop_saveset(signal,''filepath'',fp,''filename'',[fn fx]);');
else
    save(fname,'signal','-mat');
end
fprintf('done.\n');