function result = set_concat(varargin)
% Concatenate continuous signals across time.
% Result = set_joinepos(Set1, Set2, ...)
%
% In:
%   SetK   : The k'th data set to concatenate.
%
% Out:
%   Result : A new data set that is the concatenation of all input sets. The following changes are made:
%            * .data and all other time-series fields are concatenated across time (2nd dimension)
%            * .event is joined and .latency fields are updated appropriately
%            * .xmax/.pnts are updated
%
% Notes:
%   This function returns a new data set with meta-data set to that of the first input set, and the
%   time series field sjoined across all sets. No checks for meta-data consistency are done. There
%   is a heavy-duty function for merging inconsistent sets called set_merge, which can merge cats
%   and dogs. This function does not attempt to keep miscellaneous EEGLAB meta-data consistent,
%   including: setname,filename,filepath,subject,group,condition,session,comments,urevent,reject,stats,history,etc
%
% Examples:
%   % concatenate data sets eegA, eegB and eegC across time
%   eeg = set_concat(eegA,eegB,eegC)
%
% See also:
%   set_joinepos, set_merge
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-03-31

% set_joinepos_version<1.0> -- for the cache

if ~exp_beginfun('editing') return; end

declare_properties('name','Concatenate','independent_channels',true,'independent_trials',false);

if ~isempty(varargin)
    if any(cellfun(@(x)isfield(x,'epoch') && ~isempty(x.epoch),varargin))
        error('Only continuous data can be concatenated with set_concat -- use set_joinepos for epoched data.'); end
    result = varargin{1};
    if length(varargin) > 1
        % concatenate time series fields
        for field = utl_timeseries_fields(result)
            data = cellfun(@(x)x.(field{1}),varargin,'UniformOutput',false);
            result.(field{1}) = cat(2,data{:}); 
            if isempty(result.(field{1}))
                result.(field{1}) = []; end
        end
        % count events, epochs and samples in each set
        event_count = cellfun(@(x)length(x.event),varargin);
        sample_count = cellfun(@(x)x.pnts,varargin);
        % concatenate .event and .epoch fields
        event = cellfun(@(x)x.event,varargin,'UniformOutput',false); result.event = [event{:}];
        % shift event latencies based on cumulative sample counts
        if ~isempty(result.event)
            [result.event.latency] = arraydeal([result.event.latency]+replicate(cumsum(sample_count)-sample_count,event_count)); end
        % update misc fields
        [result.nbchan,result.pnts,result.trials,extra_dims] = size(result.data); %#ok<NASGU>
        result.xmax = result.xmin + (result.pnts-1)/result.srate;
    end
else
    result = eeg_emptyset;
end

exp_endfun;

function result = replicate(values,counts)
% Replicate each element Values(k) by Count(k) times.
result = zeros(1,sum(counts));
k = 0;
for p=find(counts)
    result(k+(1:counts(p))) = values(p);
    k = k+counts(p);
end
