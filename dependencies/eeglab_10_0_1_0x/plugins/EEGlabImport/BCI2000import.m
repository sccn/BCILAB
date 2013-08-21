function EEG = BCI2000import(filename,pickevents,mergeposition,concatruns,maxevents)
% Import a BCI2000 data set
%
% In:
%   Filename : optional file name to open (file dialog will open otherwise)
%
%   PickEvents : let the user pick the events to import (default: true)
%                (if false, all events will be imported)
%
%   MergePosition : merge the "position" id of a BCI2000 event into
%                   the event type, as _# (default: false)
%
%   ConcatRuns : whether to concatenate and return all runs related to a single selected
%                session file (default: false)
%
%   MaxEvents : maximum number of events for a State type to be accepted into the data set
%               (default: 5000)
%
% Out:
%   EEG : EEGLAB data set
%
% (C) 2000-2009, BCI2000 Project
% http://www.bci2000.org
%
% Adapted by Christian Kothe (Swartz Center for Computational Neuroscience, UCSD)

if ~exist('maxevents','var') maxevents = 5000; end

warning off MATLAB:intMathOverflow

addpath(genpath('..\..\functions'));
EEG = [];
%% load the data
if ~exist('filename','var') || isempty(filename)
    % use the GUI
    [fname, fpath] = uigetfile('*.dat','multiselect','on');
    if fpath == 0
        return;
    end
    if ~iscell(fname)
        tmp = fname;
        clear fname;
        fname{1} = tmp; clear tmp;
    end
    for i=1:length(fname)
        fname{i} = [fpath,fname{i}];
    end
else
    % do it directly
    fname = filename;
    if ~iscell(fname)
        fname = {fname}; 
    end
end

% concatenate runs?
if exist('concatruns','var') && concatruns
    % just one file picked?
    if length(fname) == 1
        [p,n,x] = fileparts(fname{1});
        fname = {};
        ids = [];
        % correct filename format?
        if length(n) > 2 && (n(end-2) == 'R') && ~isempty(str2num(n(end-1:end)))
            d = dir(p);
            % list all files that begin appropriately
            cands = {d(strncmp({d.name},n,length(n)-2)).name};
            for c=1:length(cands)
                [pp,nn,xx] = fileparts(cands{c});
                % check if the remainder is actually numeric, and if the extension 
                % is identical to the one of the picked file
                id = str2num(nn(length(n)-1:end));
                if ~isempty(id) && strcmp(xx,x)
                    fname{end+1} = [p filesep cands{c}];
                    ids(end+1) = id;
                end
            end
        end
        % sort them in ascending order
        [ids,order] = sort(ids);
        fname = fname(order);
        if any(diff(ids) > 1)
            disp('Warning: At least one of the runs is missing.'); 
        end
    end
end


files = struct('name',fname);
startPos = 0;
evCount = 1;
signal = [];
states = [];
[signal, states, parms] = load_bcidat(files.name);
signal = single(signal);

stateNames = fieldnames(states);
commonStates = {'TargetCode','ResultCode','Dwelling','Feedback',...
    'IntertrialInterval'};
sel = [];
for i=1:length(stateNames)
    a = strFindCell(commonStates,stateNames{i},'exact');
    if a > 0
        sel = [sel, i];
    end
end

if ~exist('pickevents','var') || pickevents
    [sel, ok] = listdlg('liststring',stateNames,...
        'selectionmode','multiple','promptstring','Select BCI Events','initialvalue',sel);
    if ok == 0
        data = [];
        return;
    end
else
    sel = 1:length(stateNames);
end

selectedEvents = stateNames(sel);
for k=1:length(stateNames)
    if ismember(stateNames{k}, selectedEvents)
        continue;
    end
    states = rmfield(states, stateNames{k});
end
states = compressData(states);
stateNames = fieldnames(states);

EEG = eeg_emptyset;
if ~isempty(parms.SamplingRate.NumericValue)
    EEG.srate = parms.SamplingRate.NumericValue;
elseif ~isempty(parms.SamplingRate.Value)
    EEG.srate = str2num(strrep(parms.SamplingRate.Value{1},'Hz',''));
end
EEG.nbchan = size(signal,2);
EEG.pnts = size(signal,1);
EEG.trials = 1;
EEG.data = signal';
try
    EEG.chanlocs = struct('labels',parms.ChannelNames.Value(1:EEG.nbchan));
    if length(parms.ChannelNames.Value) > EEG.nbchan
        disp('Warning: There are more channel labels in the BCI2000 file than there were recorded channels...'); end
catch,end
EEG = eeg_checkset(EEG);

evCount = 1;
for i=1:length(stateNames)    
    st = getfield(states,stateNames{i});
    if length(st.latency) <= maxevents
        for ev=1:length(st.latency)
            EEG.event(evCount).latency = st.latency(ev);
            EEG.event(evCount).position = st.value(ev);
            EEG.event(evCount).type = stateNames{i};
            evCount = evCount+1;
        end
    end
end

EEG = eeg_checkset(EEG, 'eventconsistency');
EEG = pop_editeventvals(EEG, 'sort', {'latency',[0]});

if exist('mergeposition','var') && mergeposition
    for k=1:length(EEG.event)
        tmp = EEG.event(k).position;
        if ~isempty(tmp)            
            if isnumeric(tmp)
                tmp = num2str(tmp);
            end
            EEG.event(k).type = [EEG.event(k).type '_' tmp];
        end
    end
end

%------------------------------
%-----------------------------------------------------------------------
function pos = strFindCell(str, toFind, option)

pos = 0;
if ~exist('option')
    option = 'all';
end

if ~iscellstr(str)
    error('Input to strFindCell must be a cell array of strings.');
    pos = 0;
    return;
end

switch option
    case 'all'
        count = 1;
        for i=1:length(str)
            if ~isempty(strfind(str{i}, toFind))
                pos(count) = i;
                count = count+1;
            end
        end

    case 'exact'
        for i=1:length(str)
            if (strcmp(str{i}, toFind))
                pos = i;
                return;
            end
        end
end

% -------------------------------------------------------------
function data = compressData(data, name)

compressStates = {'TargetCode','ResultCode','IntertrialInterval','Feedback','Dwelling','StimulusCode'};

if isstruct(data)
    stateNames = fieldnames(data);
    for i=1:length(stateNames)
        if ismember(stateNames{i}, compressStates)
            value = getfield(data, stateNames{i});
            if isstruct(value) continue; end
            data = setfield(data, stateNames{i}, compressData(value, stateNames{i}));
            clear value;
        else
            % use some trickery to see if we should try to compress it or
            % not
            value = getfield(data,stateNames{i});
            if isstruct(value) continue; end
            
            % assume that the type of data that we will use as events is
            % less than 5% of the total length of the data
            U = unique(value);
            if length(U) > 1 && length(U) < length(value)*.05
                data = setfield(data, stateNames{i}, compressData(value, stateNames{i}));
            else
                disp(['Warning: ', stateNames{i}, ' is not a valid state for EEGlab, and is not being imported.']);
                data = rmfield(data, stateNames{i});
            end
        end
    end
else
    % do a run-length encode into events
    lat = find([1; diff(data(:))] ~= 0);
    dur = diff([lat;numel(data)+1])';
    data = struct('latency',lat,'value',data(lat),'duration',dur');
end
