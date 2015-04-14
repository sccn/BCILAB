% pop_WearableSensing() - import WearableSensing files into EEGLAB
%
%       This function is intended to import .CSV files produced by
%       Wearable Sensing's DSI streaming software into EEGLAB. This 
%	is version 1.07.
%
% Usage:
%   >> OUTEEG = pop_WearableSensing; % pop up window
%   >> OUTEEG = pop_WearableSensing( filename, channels, type);
%
% Inputs:
%   filename - [string] file name
%
% Optional inputs:
%   'channels'   - [integer array] list of channel indices
%
%   'blockrange' - [min max] integer range of data blocks to import, in seconds.
%                  Entering [0 3] will import the first three blocks of data.
%                  Default is empty -> import all data blocks.
%
%  'importevent' - ['on'|'off'] import events. Default if 'on'.'
%
%  'lowpass'     - [integer]  This does a lowpass filter on the data. We
%                   recommend [70] hz.  If undefined, does not lowpass filter
%                   the data.
%
% Outputs:
%   OUTEEG   - EEGLAB data structure
%
% Author: Steven Pillen, steven@wearablesensing.com
% Copyright (C) Wearable Sensing, 2015.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
function EEG = pop_WearableSensing(filename, varargin);
EEG = [];


if nargin < 1
    % ask user
    [filename, filepath] = uigetfile('*.csv', 'Choose a data file -- pop_WearableSensing()'); %%% this is incorrect in original version!!!!!!!!!!!!!!
    drawnow;
    
    if filename == 0 return; end;
    filename = [filepath filename];
    
    
    % open file to get infos
    % ----------------------
    
    CSV = importdata2(filename);
    
    
    
    %specify data (we want to see the total time length here)
    % Identify montage
    for CheckChannel = 1:length(CSV.textdata)
        if any(strfind(CSV.textdata{CheckChannel}, 'Reference'))
            disp('Reference identified.')
            Reference = CSV.textdata{13}(23:24);
            if isletter(CSV.textdata{13}(25)) || ismember(CSV.textdata{13}(25),{'1' '2' '3' '4' '5' '6' '7' '8' '9' '0'})
                Reference = [Reference, CSV.textdata{13}(25)];
            end
        end
        
        if any(strfind(CSV.textdata{CheckChannel}, 'Sensor Position'))
            while  double(CSV.textdata{CheckChannel}(end)) == 10||double(CSV.textdata{CheckChannel}(end)) == 13
                CSV.textdata{CheckChannel}(end) = [];
            end
            tmp = regexp(CSV.textdata{CheckChannel},'([^ ,]*)','tokens');
            a = cat(1,tmp{:})';
            Montage = [a(4:end) {'Trigger'}];
            Data = CSV.data(:, 2:length(Montage)+1);
            
        else
            MontageStart = strfind(CSV.colheaders, 'ime');
            MontageStart = find(not(cellfun('isempty', MontageStart)));
            MontageStart = MontageStart(1)+1;
            
            MontageEnd = strfind(CSV.colheaders, 'rigger');
            MontageEnd = find(not(cellfun('isempty', MontageEnd)));
            MontageEnd = MontageEnd(1);
            Montage = CSV.colheaders(MontageStart:MontageEnd);
            Data = CSV.data(:, MontageStart:MontageEnd);
            
        end
    end
    if ~exist('Reference')
            Reference = 'Pz';
            disp('No reference identified.  Defaulting to Pz.')
    end
    %Here, we insert the reference channel.
    a = size(Data);
    Data(:,a(2)+1) = Data(:,a(2));
    Data(:,a(2)) = Data(:,a(2))*0;
    
    a = size(Montage);
    Montage(:,a(2)+1) = Montage(:,a(2));
    Montage{a(2)} = Reference;
    
    % %         Montage = []
    % %         for montagestep = 3:length(a)
    % %             Montage = [Montage, a{montagestep}, ' '];
    % %         end
    %     end
    
    
    
    
    
    %Sampling rate
    
    for FreqFinder = 1:length(CSV.textdata)
        
        if any(strfind(CSV.textdata{FreqFinder}, 'Sample'))
            if any(char(regexp(CSV.textdata{FreqFinder},'\d\d\d', 'match')))
                sratepart = char(regexp(CSV.textdata{FreqFinder},'\d\d\d', 'match'));
                Srate = str2double(sratepart);
                break
            elseif any(char(regexp(CSV.textdata{FreqFinder},'\d\d', 'match')))
                sratepart = char(regexp(CSV.textdata{FreqFinder},'\d\d', 'match'));
                Srate = str2double(sratepart);
                break
            end
            
            
        end
    end
    
    [Duration, Channels] = size(Data);
    DataDuration = (Duration/Srate);
    
    
    %     disp('Reading data file header...');
    %     CSV = sopen(filename, 'r', [], 'OVERFLOWDETECTION:OFF');
    %     if ~isfield(dat, 'NRec')
    %         error('Unsuported data format');
    %     end;
    
    %DataDuration = 1; %This is a dummy value for the data's time duration.  I will want to modify this into something dynamic from the data input.
    
    uilist = { { 'style' 'text' 'String' 'Channels to include (default is all):' } ...
        { 'style' 'edit' 'string' '' } ...
        { 'style' 'text' 'String' [ 'Data range (in seconds) to read (default all [0 ' int2str(DataDuration) '])' ] } ...
        { 'style' 'edit' 'string' '' } ...
        { 'style' 'text' 'String' 'Import Event Triggers (check for "yes")' 'value' 1} ...
        { 'style' 'checkbox' 'string' '' 'value' 1 } {} ...
        { 'style' 'text' 'String' 'Lowpass FIR Filter* (For no filter, leave empty)'}...
        { 'style' 'edit' 'string' '70' } {}...
        { 'style' 'text' 'String' ['*If included, the lowpass filter must be a single value lower than the sampling rate ( ', ...
        num2str(Srate), ' )'] } ...
        { 'style' 'text' 'String' ' It is recommended that you do lowpass the data at 70 Hz'}};
    geom = { [3 1] [3 1] [3 0.35 0.5]  [3 1] 1 1 1};
    
    GUIinput = inputgui( geom, uilist, 'pophelp(''pop_WearableSensing'')', ...
        'Load data using Wearable Sensing .csv -- pop_WearableSensing()');
    if length(GUIinput) == 0 return; end;
    
    % decode GUI params
    % -----------------
    options = {};
    if ~isempty(GUIinput{1}), options = { options{:} 'channels'   eval( [ '[' GUIinput{1} ']' ] ) }; end;
    if ~isempty(GUIinput{2}), options = { options{:} 'blockrange' eval( [ '[' GUIinput{2} ']' ] ) }; end;
    if ~GUIinput{3},          options = { options{:} 'importevent'  'off'  }; end;
    if ~isempty(GUIinput{4}), options = { options{:} 'lowpass'   eval( [ '[' GUIinput{4} ']' ] ) }; end;
else %if thee are user input parametersr
    options = varargin;
    CSV = importdata2(filename);
    
    MontageStart = find(ismember(CSV.colheaders,'Time'))-1; %should be the first channel
    MontageEnd = find(ismember(CSV.colheaders,'Trigger')); %should come after the last channel.  Trigger has event values
    Montage = CSV.colheaders(MontageStart:MontageEnd);
    Data = CSV.data(:, MontageStart:MontageEnd);
end

% Make sure the input parameters make sense
% -----------------------

g = finputcheck( options, { 'blockrange'   'integer' [0 Inf]    [];
    'channels'     'integer' [0 Inf]    [];
    'lowpass'     'integer'  [0 Srate]    [];
    'importevent'  'string'  { 'on';'off' } 'on';
    }, 'pop_WearableSensing');
if ischar(g), error(g); end;

% import data
% -----------

EEG = eeg_emptyset;


%   if ~isempty(GUIinput{1}), options = { options{:} 'channels'   eval( [ '[' GUIinput{1} ']' ] ) }; end;
%     if ~isempty(GUIinput{2}), options = { options{:} 'blockrange' eval( [ '[' GUIinput{2} ']' ] ) }; end;
%     if ~isempty(GUIinput{3}), options = { options{:} 'lowpass'   eval( [ '[' GUIinput{3} ']' ] ) }; end;
%     if ~result{3},          options = { options{:} 'importevent'  'off'  }; end;

%channels, sampling rate, the datablock itself
%import events
firststep = 0;
if strcmp(g.importevent, 'on')
    events = [];
    for eventsearch = 1:length(Data(:,find(ismember(Montage,'Trigger'))))
        firststep = firststep + 1;
        
        if firststep == 1 %this is for the instance where the first value is nonzero
            continue
        end
        
        if Data(eventsearch,find(ismember(Montage,'Trigger'))) ~= 0 % if we have a nonzero trigger,
            if Data(eventsearch-1,find(ismember(Montage,'Trigger'))) == 0 %and if it is preceded by zero,
                events = [events, eventsearch]; %this means it is an event.  It is added to this list.
            end
        end
    end
    if ~isempty(events)
        for event = 1:length(events)
            EEG.event(event).type = [ num2str(Data(events(event),find(ismember(Montage,'Trigger')))) ]; %put 'trigger ', at the beginning of the brackets to label it
            EEG.event(event).latency = events(event);
        end
    end
    
end



%import channels

%Eliminate extraneous "channels" from data, montage
% 
% Data(:,find(ismember(Montage,'S17'))) = [];
% Montage(:,find(ismember(Montage,'S17'))) = [];
% 
% Data(:,find(ismember(Montage,'S18'))) = [];
% Montage(:,find(ismember(Montage,'S18'))) = [];
% 
 Data(:,find(ismember(Montage,'CM'))) = [];
 Montage(:,find(ismember(Montage,'CM'))) = [];
% 
% Data(:,find(ismember(Montage,'S21'))) = [];
% Montage(:,find(ismember(Montage,'S21'))) = [];

Data(:,find(ismember(Montage,'Trigger'))) = [];
Montage(:,find(ismember(Montage,'Trigger'))) = [];

%Eliminate the channels specified in the user input from the data and montage
if ~isempty(g.channels)
    for channel = length(Montage):-1:1
        if ~any(ismember(g.channels, channel))
            Data(:,channel) = [];
            Montage(:,channel) = [];
        end
    end
end

%limit the time range to what is specified

if ~isempty(g.blockrange)
    interval(1) = blockrange(1) * dat.SampleRate(1) + 1;
    interval(2) = blockrange(2) * dat.SampleRate(1);
    Data = Data(interval(1):interval(2),:);
end


%Now let's check the lowpass
filter = 0;
if ~isempty(g.lowpass)
    if length(g.lowpass) ~= 1 || g.lowpass < 1
        error('Lowpass value not valid. Make sure it is a single positive value, lower than the Sampling Rate.')
    else
	filter = 1;
        %This is where the fir filter script could be placed.  That will
        %come later.
    end
end

%populate sample rate, the data, and the length of the data
EEG.srate = Srate;
EEG.data = Data';
EEG.nbchan = size(EEG.data, 1);
EEG.xmax = DataDuration;
EEG.trials = 1;
EEG.pnts = size(EEG.data,2);
EEG = eeg_checkset(EEG);



%populate the Channel labels
channelindex = 1;
for label = Montage
    EEG.chanlocs(channelindex).labels = char(label);
    EEG.chanlocs(channelindex).ref = '';
    EEG.chanlocs(channelindex).theta = [];
    EEG.chanlocs(channelindex).radius= [];
    EEG.chanlocs(channelindex).X = [];
    EEG.chanlocs(channelindex).Y = [];
    EEG.chanlocs(channelindex).Z = [];
    EEG.chanlocs(channelindex).sph_theta = [];
    EEG.chanlocs(channelindex).sph_phi = [];
    EEG.chanlocs(channelindex).sph_radius = [];
    EEG.chanlocs(channelindex).type = '';
    EEG.chanlocs(channelindex).urchan = [];
    channelindex = channelindex+1;
end
EEG = eeg_checkset(EEG,'makeur');   % Make EEG.urevent field
EEG = eeg_checkset(EEG);
%Mark the referenced value
for i = 1:EEG.nbchan
    EEG.chanlocs(i).ref = Reference;
end
if filter == 1
	EEG = pop_eegfiltnew(EEG, [], g.lowpass);
end
eeglab redraw;

end

