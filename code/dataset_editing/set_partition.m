function res = set_partition(varargin)
% Obtain a subset of a data set, using an index set.
% Signal = set_partition(Signal, IndexSet, EpochBounds)
%
% A partitioning function (as required by the generic functions utl_crossval, utl_searchmodel,
% utl_nested_crossval) that takes a data set (either epoched or continuous), and index set (either 
% of the raw samples, the target events, or of the epochs), and returns a reduced (subset) data set.
% 
% In:
%   Signal   : continuous or epoched input data set
%
%   IndexSet : either [] or a vector of indices;
%               * if [], the function returns the cardinality (highest allowed index) of the index set,
%               * if a vector of indices, the function returns the data set reduced to those indices
%                 if the indices refer to target events, but the data is continuous, then the 
%                 partition function tries to respect the implied epoch bounds when partitioning
%
%   EpochBounds : only required when partitioning a continuous data set based on target events; 
%                 this is done with the help of implied epoch bounds (as in set_makepos)
%
% Out:
%   Result : Reduced data set or index set cardinality.
%
% Examples:
%   % for an epoched data set or a continuous set with target markers, get the number of trials
%   numtrials = set_partition(eeg,[])
%
%   % for an epoched data set or a continuous set with target markers, retain only trials 20:40
%   eeg = set_partition(eeg,20:40)
%
%   % for a continuous data set without target markers, get the number of samples
%   numsamples = set_partition(eeg,[])
%
%   % for a continuous data set without target markers, retain only samples 1:10000 and 50000:100000
%   eeg = set_partition(eeg,[1:10000 50000:100000])
%
% See also:
%   set_makepos, set_selepos, set_selinterval, utl_crossval
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-01

% set_partition_version<1.0> -- for the cache

if ~exp_beginfun('editing', 'memoize',0, 'fingerprint_create','passthrough') return; end

declare_properties('independent_channels',true,'independent_trials',false);

arg_define(varargin, ...
    arg_norep({'signal','Signal'}), ...
    arg({'idxset','IndexSet'},[],[], 'Index set. This is either [] (indicating that the function should return the index set cardinality), or a vector of indices (indicating that the function shall partition the data set accordingly).'), ...
    arg({'epoch_bounds','EpochBounds'},[],[], 'Epoch bounds. Only required when partitioning on a continuous data set based on target events; this is performed using the epoch bounds (as in set_makepos).'));

is_continuous = isempty(signal.epoch) && (isempty(signal.data) || size(signal.data,3) == 1); %#ok<NODEF>

if isempty(idxset) %#ok<NODEF>
    % --- calc index set size ---
    
    if ~is_continuous
        % epoched data set: return the number of epochs
        res = size(signal.data,3);
    else        
        % continuous data set
        if ~isempty(signal.event) && isfield(signal.event,'target')
            targetmask = ~cellfun('isempty',{signal.event.target});
            if any(targetmask)
                % return the number of target events
                res = nnz(targetmask);
                return
            end
        end
        % otherwise return the number of samples
        res = size(signal.data,2);
    end
else
    % --- partition the data ---
    
    if ~is_continuous
        % set is epoched: strictly partition in terms of epochs
        res = exp_eval(set_selepos(signal,idxset));
    else
        % set is continuous
        if ~isempty(signal.event) && isfield(signal.event,'target')
            targetmask = ~cellfun('isempty',{signal.event.target});
            if any(targetmask)                
                % ... and it has target events: strictly partition in terms of target events
                % this partitioning is perfectly consistent with what would happen if the data were epoched
                % (assuming that target markers with boundary-crossing epochs were already removed in
                %  set_targetmarkers, and assuming that the given epoch bounds are an upper bound of those
                %  used in set_makepos)
                if isempty(epoch_bounds) %#ok<NODEF>
                    disp_once('Warning: Attempting to partition a continuous data set based on target events, but no EpochBounds argument is present. Assuming 1-sample epochs.'); 
                    ival = 0;
                else
                    % make sure that the epoch bounds include the event itself (and add a bit more slack 
                    % to not lose it later even under resampling and rounding errors...)
                    epoch_bounds(1) = min(epoch_bounds(1),-0.25);
                    epoch_bounds(2) = max(epoch_bounds(2),+0.25);
                    ival = round(epoch_bounds(1)*signal.srate) : round(epoch_bounds(2)*signal.srate);
                end
                
                % check for sorting issues
                if ~issorted(idxset)
                    warn_once('Continuous-time data cannot reasonably be used with randomized partitioning; consider using chronological/blockwise partitioning! The index set will be sorted in the following.');
                    idxset = sort(idxset);
                end

                % build a logical mask that marks all retained intervals, starting with all false
                select_mask = false(1,size(signal.data,2));
                deselect_mask = false(1,size(signal.data,2));
                
                % find all target event latencies and the selected latencies
                targlats = round([signal.event(targetmask).latency]);                
                selected = targlats(idxset);
                
                % find the deselected latencies
                deselmask = true(size(targlats));
                deselmask(idxset) = false;
                deselected = targlats(deselmask);
                
                % add to the selected & deselected masks
                select_mask(min(length(select_mask),max(1,vec(bsxfun(@plus,ival,selected'))))) = true;
                deselect_mask(min(length(select_mask),max(1,vec(bsxfun(@plus,ival,deselected'))))) = true;
                
                % warn if the ranges overlap
                if any(select_mask & deselect_mask)
                    disp_once('Note: the epochs for some selected markers are overlapping with epochs of deselected markers. Make sure that you are using long enough safety margins if perfoming cross-validation.'); end
                
                % identify gaps in the select_mask
                gaps = find(diff([false select_mask false]));
                gaps = reshape(gaps(2:end-1),2,[])';
                gaps(:,2) = gaps(:,2)-1;
                
                % find all gaps that do not contain a deselected sample and fuse those...
                for g=1:size(gaps,1)
                    gap = gaps(g,1):gaps(g,2);
                    if ~any(deselect_mask(gap))
                        select_mask(gap) = true; end
                end
                
                % disable the now de-selected target events (settings their target flag back to empty)
                % (note that partitioning is strictly in terms of target markers)
                targetinds = find(targetmask);
                [signal.event(targetinds(deselmask)).target] = deal([]);
                
                % finally, re-extract the intervals and perform the actual partitioning
                retain_intervals = reshape(find(diff([false select_mask false])),2,[])';
                retain_intervals(:,2) = retain_intervals(:,2)-1;
                res = exp_eval(set_selinterval(signal,retain_intervals,'samples'));
                return;
            end
        end
        
        % otherwise partition the raw data based on sample intervals
        
        % create a logical mask from the index set
        mask = false(1,size(signal.data,2)); mask(idxset) = true;
        % find the inclusive, 1-based intervals of non-zeros
        retain_intervals = reshape(find(diff([false mask false])),2,[])';
        retain_intervals(:,2) = retain_intervals(:,2)-1;
        res = exp_eval(set_selinterval(signal,retain_intervals,'samples'));
    end
end

exp_endfun;
