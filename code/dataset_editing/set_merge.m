function [data,idxmap] = set_merge(varargin)
% Merge epoched EEGLAB data sets across trials or time in a fault-tolerant way.
% [Merged,IndexMap] = set_merge(Set-#1, Set-#2, ...)
%
% In:
%   Set-#k   : data set #k
%
% Out:
%   Merged   : The merged epoched data set
%
%   IndexMap : a mapping from trial index (in the merged set), to data set index 
%              (in the list of sets supplied)
%
% Notes:
%   This function merges data sets either across trials (if epoched) or across time (if continuous)
%   in a manner that supports inconsistent meta-data. For data that is known to be consistent, use
%   set_joinepos or set_concat, which are much faster.
%
% Examples:
%   % merge data sets eegA, eegB and eegC across trials
%   eeg = set_merge(eegA,eegB,eegC)
%
% See also:
%   set_concat, set_joinepos
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-03-31

% set_merge_version<1.0> -- for the cache

if ~exp_beginfun('editing') return; end

declare_properties('name','MergeSets','independent_channels',true,'independent_trials',true);

idxmap = {};
data = {};
if ~isempty(varargin)
    % identify non-empty sets & create the index map
    for i=1:length(varargin)
        if isfield(varargin{i},'data')
            data{end+1} = varargin{i}; %#ok<AGROW>
            idxmap{end+1} = i*ones(1,size(varargin{i}.data,3)); %#ok<AGROW>
        end
    end
    
    % do the merging
    data = merge(data{:});
    idxmap = [idxmap{:}];
end

exp_endfun;

% merge recursively to avoid growing big arrays incrementally
function X = merge(varargin)
if nargin > 1
    A = merge(varargin{1:floor(nargin/2)});
    B = merge(varargin{(floor(nargin/2)+1):end});
    try        
        X = pop_mergeset(A,B);
    catch err
        disp(['The two data sets ' hlp_tostring(A) ' and ' hlp_tostring(B) ' cannot be merged; reason:']);
        rethrow(err);
    end
else
    X = varargin{1};
end
