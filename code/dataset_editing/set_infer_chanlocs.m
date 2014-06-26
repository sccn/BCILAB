function data = set_infer_chanlocs(data,disambiguation_rule)
% Infer the chanlocs subfields (positions and type) from labels.
% Data = set_infer_chanlocs(Data,DisambiguationRule)
%
% In:
%   Data   : some EEGLAB data set with a chanlocs field or the chanlocs field itself
%
%   DisambiguationRule : rule to apply to disambiguate montages when a given set of channel labels
%                        match multiple montage files similarly well. Can be:
%                        * 'coverage' : picks the one with maximum coverage
%                        * 'first' : picks the first one in the list of candidates (out of those that 
%                                    have reasonable coverage)
%                        * 'deduce' : uses correlation analysis to deduce the best fit (can be slow)
%                        * filename : uses the montage with the given file name
%                        (default: 'deduce')
%
% Out:
%   Data   : original set with updated chanlocs field, or the updated chanlocs themselves
%
% Notes:
%   parameters cannot be passed by name to this function.
%
% Examples:
%   % given a data set with labeled channels (but no locations), fill in standard channel locations
%   eeg = set_infer_chanlocs(eeg)
%
%   % given a .chanlocs field of a data set with only labels, fill in standard channel locations
%   chanlocs = set_infer_chanlocs(eeg.chanlocs)
%
%   % given a cell array of channel labels, obtain a chanlocs structure with locations filled in
%   chanlocs = set_infer_chanlocs({'C3','Cz','C4'})
%
%   % given a chanlocs file in a format recognized by EEGLAB, obtain a corresponding chanlocs struct
%   chanlocs = set_infer_chanlocs('mylocs.sfp')
%
% See also:
%   io_loadset
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-11-05

% available cap files
caps = {'resources:/caps/Standard-10-5-Cap385.cap','resources:/caps/sccn_LSIE_cap128.cap','resources:/caps/sccn_BEM_coregistered_128_v2.cap','resources:/sccn_BEM_coregistered_256_v1.cap','resources:/caps/sccn_BEM_coregistered_128_v1.cap'};

if nargin < 2
    disambiguation_rule = 'deduce'; end

% obtain the chanlocs
signal = [];
if isfield(data,{'head','parts'})
    data = exp_eval(data); end
if isfield(data,'chanlocs')
    locs = data.chanlocs;
    signal = data;
    if ~isfield(data.chaninfo,'labelscheme')
        data.chaninfo.labelscheme = '10-20'; end
else
    locs = data;
end

% channel labels given as a cell array
if iscellstr(locs)
    locs = struct('labels',locs); end;
if ischar(locs)
    locs = readlocs(env_translatepath(locs)); end
if ~isstruct(locs)
    error('The given data is not in an appropriate format (expected either a file name, a cell array of labels, a .chanlocs struct, or an EEGLAB dataset struct).'); end

% --- find best-matching cap if some locs missing ---

if (~all(isfield(locs,{'X','Y','Z'})) || all(cellfun('isempty',{locs.X}))) && ~isempty(locs)
    
    % for each cap file, try to match channel labels
    for c = 1:length(caps)
        try
            cap = caps{c};
            capdata{c} = hlp_microcache('chanlocs',@io_load,cap,'-mat');
            locdb = capdata{c}.CAP;
            fitlocs{c} = locs;
            switch capdata{c}.MATCH
                case {'exact','nocase_exact'}
                    % exact matching (both cases supported)
                    if ~isempty(strfind(capdata{c}.MATCH,'nocase'))                        
                        [found,idx_in_locdb,idx_in_res] = intersect(lower({locdb.labels}),lower({locs.labels})); %#ok<ASGLU>
                    else
                        [found,idx_in_locdb,idx_in_res] = intersect({locdb.labels},{locs.labels}); %#ok<ASGLU>
                    end
                    for fname = fieldnames(locdb)'
                        [fitlocs{c}(idx_in_res).(fname{1})] = locdb(idx_in_locdb).(fname{1}); end
                    [fitlocs{c}(idx_in_res).type] = deal('EEG');
                case {'substr','nocase_substr'}
                    % sub-string matching (both cases supported)
                    [found,idx_in_res,idx_in_locdb] = hlp_microcache('chanlocs',@do_matching,fitlocs{c},locdb,~isempty(strfind(capdata{c}.MATCH,'nocase'))); %#ok<ASGLU>
                    for fname = fieldnames(locdb)'
                        [fitlocs{c}(idx_in_res).(fname{1})] = locdb(idx_in_locdb).(fname{1}); end %#ok<*AGROW>
                    [fitlocs{c}(idx_in_res).type] = deal('EEG');
            end
            fraction_found(c) = 1-mean(cellfun('isempty',{fitlocs{c}.X}));
        catch e
            fprintf('Not matching channel labels with cap %s due to error: %s',cap,hlp_handleerror(e));
            fraction_found(c) = 0;
        end
    end

    % determine whether we have a clear-cut case
    [sorted_fractions,order] = sort(fraction_found,'descend');
    if max(sorted_fractions) < 0.5
        fprintf('WARNING: only %.0f%% of your channel labels match any known cap design; inferred locations might be wrong. Please consider adding your cap to set_infer_chanlocs.',100*max(sorted_fractions)); end
    if sorted_fractions(1) >= 1.5*sorted_fractions(2)
        % at least twice as many channels matched for this cap as for the next best one: keep it
        bestidx = order(1);
    elseif ~isempty(signal) && size(signal.data,1)==length(locs)
        % not a clear-cut case but we have a signal to work with
        switch disambiguation_rule
            case 'coverage'
                fprintf('Multiple known cap designs match to your data''s channel labels; picking the one with greatest coverage...\n');
                bestidx = order(1);
            case 'first'
                fprintf('Multiple known cap designs match to your data''s channel labels; picking the first one in the list...\n');
                bestidx = 1;
            case 'deduce'
                fprintf('Multiple known cap designs match to your data''s channel labels; determining the best fit...\n');
                quality = hlp_diskcache('montage_quality',@fit_quality,signal,sorted_fractions,order,fitlocs);
                % assess quality of best fit
                if max(quality)<0.3
                    fprintf('WARNING: Did not find a cap design that matches your data. You need to add channel locations for your cap to set_infer_chanlocs to use any function that requires these locations.\n'); 
                    bestidx = NaN;
                else
                    if max(quality)<0.5
                        fprintf('WARNING: The agreement between your data and the available cap designs is low and the chosen locations might be partially wrong. Either your data is too noisy for a reliable determination or your cap design differs from the known ones.\n');  end
                    bestidx = order(argmax(quality));
                end
            otherwise
                % assume that a montage file name was given
                match = ~cellfun('isempty',strfind(caps,disambiguation_rule));
                if nnz(match) == 1
                    bestidx = find(match);
                elseif nnz(match) > 1
                    fprintf('WARNING: multiple montages match your montage name: picking the first one.\n');
                    bestidx = find(match,1);
                else
                    % assume that this is a file name and try to load it
                    bestidx = NaN;
                    tmp = getfield(hlp_microcache('chanlocs',@io_load,disambiguation_rule,'-mat'),'CAP');
                    if ~isequal({tmp.labels},{locs.labels})
                        error('Currently the given cap file must match the available channel labels exactly.'); end
                end
        end
    elseif sorted_fractions(1)>sorted_fractions(2)
        % not a clear-cut case but minor evidence in favor of best cap
        fprintf('NOTE: multiple cap designs are compatible with your channel labels, using the best-matching one; inferred locations might be wrong.\n');
        bestidx = order(1);
    else
        % no clear winner
        fprintf('WARNING: multiple cap designs match your channel labels equally well, using the first cap in the list; inferred locations may be wrong.\n');
        bestidx = order(1);
    end
    
    % use the picked locations
    if ~isnan(bestidx)
        fprintf('  using montage %s.\n',caps{bestidx});
        locs = fitlocs{bestidx};
        if isstruct(data)
            data.chaninfo.labelscheme = capdata{bestidx}.SCHEME;
            data.chaninfo.nosedir = capdata{bestidx}.NOSE;
        end
    end
        
    try
        % optimize the head center
        if  mean(cellfun('isempty',{locs.X})) < 0.75
            locs = hlp_microcache('chanlocs',@pop_chanedit,locs,'eval','chans = pop_chancenter( chans, [],[]);'); end
    catch e
        fprintf('Failed trying to optimized the head center location for your cap montage due to error: %s\n',hlp_handleerror(e));
    end
end    

% --- do post-processing ---

if ~isfield(locs,'labels')
    for k=data.nbchan:-1:1
        locs(k).labels = num2str(k); end
end
if ~isfield(locs,'type')
    locs(1).type = []; end
if ~isfield(locs,'theta')
    locs(1).theta = []; end

% reset 'unknown' channel types to empty
for k=1:length(locs)
    if strcmpi(locs(k).type,'unknown')
        locs(k).type = ''; end
end

% if some channel types are missing, try to look them up
empty_type = cellfun('isempty',{locs.type});
if any(empty_type)
    % assign all successfully retrieved EEG channels
    [locs(~cellfun('isempty',{locs.theta}) & empty_type).type] = deal('EEG');

    % assign the most common prefixes
    prefixes = {'eo','EOG','em','EMG','ec','ECG','ex','EXG','gsr','GSR'};
    for p=1:2:length(prefixes)
        [locs(strncmpi(prefixes{p},{locs.labels},length(prefixes{p})) & empty_type).type] = deal(prefixes{p+1}); end
    
    % also look for occurrence of EO, EM, ECG, EX, and GSR patterns
    patterns = {'EO','EOG','EM','EMG','ECG','ECG','EX','EXG','GSR','GSR','RE','EOG','LE','EOG'};
    for p=1:2:length(patterns)
       [locs(~cellfun('isempty',strfind({locs.labels},patterns{p})) & empty_type).type] = deal(patterns{p+1}); end
    
   % assign type 'EEG' to all numerically labeled channels (default assumption)
    [locs(cellfun(@(x)~isempty(str2num(x)),{locs.labels})).type] = deal('EEG'); %#ok<ST2NM>
   
    % assign 'EEG' to all unknown channels that are preceded and followed by EEG channels...
    filled_in = false;
    for k = find(cellfun('isempty',{locs.type}))
        % check if there is an EEG channel before this one (besides other unknowns)
        eeg_before = false;
        for j=k-1:-1:1
            if ~isempty(locs(j).type)
                if strcmp(locs(j).type,'EEG')
                    eeg_before = true; end
                break;
            end
        end
        % check if there is an EEG channel after this one (besides other unknowns)
        eeg_after = false;
        for j=k+1:length(locs)
            if ~isempty(locs(j).type)
                if strcmp(locs(j).type,'EEG')
                    eeg_after = true; end
                break;
            end
        end        
        if eeg_before && eeg_after
            locs(k).type = 'EEG'; 
            filled_in = true;
        end        
    end
    if filled_in
        disp('Note: some channels could not clearly be recognized as EEG (not in the 10-5 system), but were preceded and followed by EEG channels; assuming they are EEG.'); end
    
    % assign 'unknown' to the rest
    [locs(cellfun('isempty',{locs.type})).type] = deal('unknown');
end

if isfield(data,'chanlocs')
    data.chanlocs = locs;
else
    data = locs;
end


function [found,idx_in_res,idx_in_locdb] = do_matching(fitlocs,locdb,use_lower)
found = {};
idx_in_res = [];
idx_in_locdb = [];
fitlabels = {fitlocs.labels};
dblabels = {locdb.labels};
if use_lower
    fitlabels = lower(fitlabels);
    dblabels = lower(dblabels);
end
for k=1:length(locdb)
    match = ~cellfun('isempty',strfind(fitlabels,dblabels{k}));
    if any(match)
        found{end+1} = locdb(k).labels;
        idx_in_res(end+1) = find(match,1);
        idx_in_locdb(end+1) = k;
    end
end
                    

function quality = fit_quality(signal,sorted_fractions,order,fitlocs)
% first calculate bandpass-filtered signal X
attenuation = 80;
frequencies = min(1,2*[0.5 2 45 50]/signal.srate);        
w = design_kaiser(frequencies(1),frequencies(2),attenuation,true);
B = design_fir(length(w)-1,[0 frequencies 1],[0 0 1 1 0 0],[],w);
try
    % try to process all data at once
    X = filtfilt_fast(B,1,signal.data')';
catch
    % fall back to channel-by-channel processing (if out of memory)
    for c=signal.nbchan:-1:1
        X(c,:) = filtfilt_fast(B,1,signal.data(c,:)')'; end
end

% for each cap design that matches reasonably well...
for k=1:nnz(sorted_fractions(1) < 2*sorted_fractions)
    o = order(k);
    % get matched locations
    tmp = fitlocs{o};
    [x,y,z] = deal({tmp.X},{tmp.Y},{tmp.Z});
    usable_channels = find(~cellfun('isempty',x) & ~cellfun('isempty',y) & ~cellfun('isempty',z));
    positions = [cell2mat(x(usable_channels));cell2mat(y(usable_channels));cell2mat(z(usable_channels))];
    % calculate channel interpolation matrix
    M = hlp_diskcache('montages',@interpolation_matrix,positions);
    % calculated interpolated signal Y
    Y = M*X(usable_channels,:);
    % calculated median of windowed correlation between X and Y
    window_len = 1*signal.srate;
    wnd = 0:window_len-1;
    offsets = round(1:window_len:size(Y,2)-window_len);
    W = length(offsets);
    corrs = [];
    for o=W:-1:1
        XX = X(usable_channels,offsets(o)+wnd)'; XX = bsxfun(@minus,XX,mean(XX));
        YY = Y(:,offsets(o)+wnd)'; YY = bsxfun(@minus,YY,mean(YY));
        corrs(:,o) = sum(XX.*YY)./(sqrt(sum(XX.^2)).*sqrt(sum(YY.^2)));
    end
    chancorr = median(corrs,2);
    % use average self-correlation across channels as quality index for this cap design
    quality(k) = mean(chancorr);
end


function M = interpolation_matrix(positions)
% interpolation_matrix_version<1.0> -- for the cache
fprintf('  calculating interpolation matrix for your montage; this is a one-time process');
range = 1:length(positions);
M = zeros(length(positions));
for k=range
    M(k,range~=k) = real(sphericalSplineInterpolate(positions(:,range~=k),positions(:,k)))'; 
    fprintf('.');
end
fprintf('\n');

