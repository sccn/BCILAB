function fname = io_saveset(signal,fname,overwrite)
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

EEG = signal;

if nargin < 3
    overwrite = false; end
if ~ischar(fname)
    error('The given file name must be a string.'); end
if size(fname,1) ~= 1
    error('The given file name must be a non-empty row vector of characters.'); end

% append file ending if necessary
if length(fname)<4 || (~strcmp(fname(end-3:end),'.set') && ~strcmp(fname(end-3:end),'.mat'))
    fname = [fname '.set']; end

fprintf('Saving set to %s...',fname);
io_save(fname,'-mat','-makedirs',quickif(~overwrite,'-nooverwrite',''),quickif(strcmp(fname(end-3:end),'.set'),'EEG','signal'));
fprintf('done.\n');
