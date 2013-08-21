function [data,coords] = set_infer_chanlocs(data)
% Infer the chanlocs subfields (positions and type) from labels.
% [Data] = set_infer_chanlocs(Data)
%
% In:
%   Data   : some EEGLAB data set with a chanlocs field or the chanlocs field itself
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

% obtain the chanlocs
if isfield(data,{'head','parts'})
    data = exp_eval(data); end
if isfield(data,'chanlocs')
    locs = data.chanlocs;
    if ~isfield(data.chaninfo,'labelscheme')
        data.chaninfo.labelscheme = '10-20'; end
else
    locs = data;
end

% channel labels given as a cell array
if iscell(locs)
    locs = struct('labels',locs); end;

if ischar(locs)
    locs = readlocs(env_translatepath(locs)); end

% if no channel locations present, try to look them up
if (~all(isfield(locs,{'X','Y','Z'})) || all(cellfun('isempty',{locs.X}))) && ~isempty(locs)
    % go through a variety of caps and retain the best (i.e. largest) match
    caps = {'resources:/Standard-10-5-Cap385.sfp','resources:/sccn_LSIE_cap128.xyz','resources:/sccn_BEM_coregistered_128_v2.xyz','resources:/sccn_BEM_coregistered_256_v1.xyz','resources:/sccn_BEM_coregistered_128_v1.xyz'};
    schemes = {'10-20','sccn_128_v3','sccn_128_v2','sccn_256_v1','sccn_128_v1'};
    matching = {'nocase','substr','substr','substr'};
    nosedir = {'+X','+Y','+Y','+Y'};
    for c = 1:length(caps)
        try
            cap = caps{c};
            fitlocs{c} = locs;
            switch matching{c}
                case 'nocase'
                    % not case sensitive
                    fitlocs{c} = pop_chanedit(locs,'lookup',env_translatepath(cap));
                case 'case'
                    % case sensitive matching
                    fitlocs = locs;
                    locdb = readlocs(env_translatepath(cap));
                    [found,idx_in_locdb,idx_in_res] = intersect({locdb.labels},{locs.labels}); %#ok<ASGLU>
                    for fname = fieldnames(locdb)'
                        [fitlocs{c}(idx_in_res).(fname{1})] = locdb(idx_in_locdb).(fname{1}); end
                    [fitlocs{c}(idx_in_res).type] = deal('EEG');
                case 'substr'
                    % case-sensitive sub-string matching
                    locdb = readlocs(env_translatepath(cap));
                    found = {};
                    idx_in_res = [];
                    idx_in_locdb = [];
                    for k=1:length(locdb)
                        match = ~cellfun('isempty',strfind({fitlocs{c}.labels},locdb(k).labels));
                        if any(match)
                            found{end+1} = locdb(k).labels;
                            idx_in_res(end+1) = find(match,1);
                            idx_in_locdb(end+1) = k;
                        end
                    end
                    for fname = fieldnames(locdb)'
                        [fitlocs{c}(idx_in_res).(fname{1})] = locdb(idx_in_locdb).(fname{1}); end
                    [fitlocs{c}(idx_in_res).type] = deal('EEG');
            end
            fraction_found(c) = 1-mean(cellfun('isempty',{fitlocs{c}.X}));
            if fraction_found(c) > 0.75
                break; end
        catch
        end
    end
    try 
        bestidx = argmax(fraction_found);
        locs = fitlocs{bestidx};
        if isstruct(data)
            data.chaninfo.labelscheme = schemes{bestidx};
            data.chaninfo.nosedir = nosedir{bestidx};
        end
        % optimize the head center
        if  mean(cellfun('isempty',{locs.X})) < 0.75
            locs = pop_chanedit(locs,'eval','chans = pop_chancenter( chans, [],[]);'); end
    catch
    end
end    

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
    [locs(cellfun(@(x)~isempty(str2num(x)),{locs.labels})).type] = deal('EEG');
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
