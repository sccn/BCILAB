function result = set_joinepos(varargin)
% Join epoched signals across epochs.
% Result = set_joinepos(Set1, Set2, ...)
%
% In:
%   SetK   : The k'th data set to join.
%
% Out:
%   Result : A new data set with trials from all sets joined. The following changes are made:
%            * .data and all other time-series fields are joined across trials (3rd dimension)
%            * .epoch is joined and its .event field is updated appropriately
%            * .event is joined and its .epoch and .latency fields are updated appropriately
%            * .xmax/.trials are updated
%
% Notes:
%   This function returns a new data set with meta-data set to that of the first input set, and the
%   trials joined across all sets. No checks for meta-data consistency are done. There is a
%   heavy-duty function for merging inconsistent sets called set_merge, which can merge cats and
%   dogs. This function does not attempt to keep miscellaneous EEGLAB meta-data consistent, including:
%   setname,filename,filepath,subject,group,condition,session,comments,urevent,reject,stats,history,etc
%
%   To concatenate continuous sets across time, use set_concat.
%
% Examples:
%   % merge data sets eegA, eegB and eegC across epochs
%   eeg = set_joinepos(eegA,eegB,eegC)
%
% See also:
%   set_concat, set_merge
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-03-31
dp;

% set_joinepos_version<1.0> -- for the cache

if ~exp_beginfun('editing') return; end

declare_properties('name','JoinTrials','independent_channels',true,'independent_trials',true);

% input validation
for k=1:length(varargin)
    utl_check_fields(varargin{k},{'data','event','epoch','pnts','xmin','xmax','srate'},'input','signal');
    if ~isempty(varargin{k}.event)
        if ~isfield(varargin{k}.event,'latency')
            error('A dataset passed to set_joineos is lacking the .event.latency field.'); end
        if ~isfield(varargin{k}.event,'epoch')
            error('A dataset passed to set_joineos is lacking the .event.epoch field.'); end
        latency_numels = cellfun('prodofsize',{varargin{k}.event.latency});
        if any(latency_numels == 0)
            error('One or more of the events in the given data set have an empty .latency field, which is not permitted.'); end
        if any(latency_numels ~= 1)
            error('One or more of the events in the given data set have a .latency value that is not a scalar, which is not permitted.'); end
    end
    if ~isempty(varargin{k}.epoch)
        if ~isfield(varargin{k}.epoch,'event')
            error('A dataset passed to set_joineos is lacking the .epoch.event field.'); end
        epoch_numels = cellfun('prodofsize',{varargin{k}.event.epoch});
        if any(epoch_numels == 0)
            error('One or more of the events in the given data set have an empty .epoch field, which is not permitted.'); end
        if any(epoch_numels ~= 1)
            error('One or more of the events in the given data set have an .epoch value that is not a scalar, which is not permitted.'); end
    end
end

if ~isempty(varargin)
    result = varargin{1};
    if length(varargin) > 1
        % concatenate time series fields
        for field = utl_timeseries_fields(result)
            data = cellfun(@(x)x.(field{1}),varargin,'UniformOutput',false);
            try
                result.(field{1}) = cat(3,data{:}); 
            catch e
                % concatenation failed: produce a reasonable error message
                if ~isempty(data)
                    sizes1 = cellfun('size',data,1); sizes1 = sizes(sizes1 ~= 0);
                    sizes2 = cellfun('size',data,2); sizes2 = sizes(sizes2 ~= 0);
                    if length(sizes1) > 1 && ~all(sizes1==sizes1(1))
                        error('The time-series field .%s must have the same number of channels in each data set.',field{1}); end
                    if length(sizes2) > 1 && ~all(sizes2==sizes2(1))
                        error('The time-series field .%s must have the same number of tine points in each data set.',field{1}); end
                end
                size_info = hlp_tostring(cellfun(@size,data,'UniformOutput',false));
                error('Concatenation of time-series fields failed with error: %s (the data sizes were %s -- make sure that they are mutually compatible).',e.message,size_info);
            end
            if isempty(result.(field{1}))
                result.(field{1}) = []; end
        end
        % count events, epochs and samples in each set
        event_count = cellfun(@(x)length(x.event),varargin);
        epoch_count = cellfun(@(x)length(x.epoch),varargin);
        sample_count = cellfun(@(x)x.pnts,varargin).*epoch_count;
        % concatenate .event and .epoch fields
        event = cellfun(@(x)x.event,varargin,'UniformOutput',false); result.event = [event{:}];
        epoch = cellfun(@(x)x.epoch,varargin,'UniformOutput',false); result.epoch = [epoch{:}];
        % shift event latencies based on cumulative sample counts
        if ~isempty(result.event)
            [result.event.latency] = arraydeal([result.event.latency]+replicate(cumsum(sample_count)-sample_count,event_count));
            % shift event/epoch cross-references based on cumulative counts
            [result.event.epoch] = arraydeal([result.event.epoch]+replicate(cumsum(epoch_count)-epoch_count,event_count));
            [result.epoch.event] = chopdeal([result.epoch.event]+replicate(cumsum(event_count)-event_count,event_count),cellfun('length',{result.epoch.event}));
        end
        % update misc fields
        [result.nbchan,result.pnts,result.trials,extra_dims] = size(result.data); %#ok<NASGU>
        result.xmax = result.xmin + (result.pnts-1)/result.srate;
    end
else
    result = struct('setname','','filename','','filepath','','subject','','group','','condition','','session',[],'comments','','nbchan',0,...
        'trials',0,'pnts',0,'srate',1,'xmin',0,'xmax',0,'times',[],'data',[],'icaact',[],'icawinv',[],'icasphere',[],'icaweights',[], ...
        'icachansind',[],'chanlocs',[],'urchanlocs',[],'chaninfo',[],'ref',[],'event',[],'urevent',[],'eventdescription',{{}}, ...
        'epoch',[],'epochdescription',{{}},'reject',[],'stats',[],'specdata',[],'specicaact',[],'splinefile','','icasplinefile','', ...
        'dipfit',[],'history','','saved','no','etc',[]);
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
