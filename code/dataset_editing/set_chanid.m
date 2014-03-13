function [remaining,notfound] = set_chanid(data,lookup)
% Translate channel names into indices in a data set's channels, or in a chanlocs structure.
% [Indices,Not-Found] = set_chanid(Dataset/Chanlocs,Lookup)
%
% In:
%   Dataset/Chanlocs :   data set, or EEGLAB chanlocs structure, or cell array of channel labels
%
%   Lookup           :   cell-string array of channel names to look up, or numeric indices
%
% Out:
%   Indices          :   The indices of the specified channel names, in the order of appearance in
%                        the supplied data set/chanlocs.
%
%   Not-Found        :   indices, into lookup, of channels that were not found in the data set
%
% Note:
%   Differs from eeg_chaninds in that it does not abort when some of the channels are not present.
%   Also, the channel indices are ordered as they appear in Lookup.
%
%   Parameters cannot be passed by name to this function.
%
% Examples:
%   % get the indices of the channels 'C3' and 'Cz'
%   indices = set_chanid(eeg,{'C3','Cz'})
%
%   % get the indices of the channels 'C3' and 'Cz' from the .chanlocs field of a data set
%   indices = set_chanid(eeg.chanlocs,{'C3','Cz'})
%
%   % get the indices of the channels 'C3','Fuyooz','Cz', which returns the indices for the present 
%   % channels and returns the indices of all channels that were not found in the second output
%   [indices,notfound] = set_chanid(eeg,{'C3','Fuyooz','Cz'})
%
%   % get the indices of all present channels
%   indices = set_chanid(eeg,[])
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-05-08

% get the channel labels from the data
if ~isfield(data,'chanlocs') && isfield(data,'tracking')
    data = exp_eval(data); end
if isfield(data,'chanlocs') && isfield(data.chanlocs,'labels')
    labels = {data.chanlocs.labels};
elseif isfield(data,'labels')
    labels = {data.labels};
elseif iscell(data)
    labels = data;
else
    error('Chanlocs data structure or EEGLAB data set struct expected.');
end

if ~iscellstr(labels)
    error('The given channel labels must all be strings.'); end
labels = lower(labels);

% check for duplicates in the labels array
tmp = sort(labels);
duplicates = strcmp(tmp(1:end-1),tmp(2:end));
if any(duplicates)
    error('Your channel labels must all be unique but the following duplicates were found: %s.',hlp_tostring(tmp(duplicates))); end

if iscell(lookup)
    if ~iscellstr(lookup)
        error('If channel labels are provided as cell array they must all be strings.'); end
    lookup = lower(lookup);
    if isequal(lookup,labels)
        % fast path
        remaining = 1:length(lookup);
        notfound = [];
    else
        % got channel names: look them up from the chanlocs, but ordered according to lookup
        [x,a,b] = intersect(lookup,labels); %#ok<ASGLU>
        [x,I] = sort(a); remaining = b(I); %#ok<ASGLU>
        notfound = setdiff(1:length(lookup),a);
    end
elseif isnumeric(lookup)
	% got channel indices: filter out the invalid channels
    if isempty(lookup)
        lookup = 1:size(data.data,1); end;
	remaining = lookup(1:min(length(lookup),size(data.data,1)));
    notfound = lookup((1+size(data.data,1)):end);
end

% check for duplicates in the lookup array
tmp = sort(lookup);
duplicates = strcmp(tmp(1:end-1),tmp(2:end));
if any(duplicates)
    error('Your channel labels/indices must all be unique but the following duplicates were found: %s.',hlp_tostring(tmp(duplicates))); end
