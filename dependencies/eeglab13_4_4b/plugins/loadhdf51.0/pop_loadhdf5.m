% pop_loadhdf5() - Read g-tec 24bit EEG file from g-recorder
%
% Usage:
%   >> EEG = pop_loadhdf5;                                      % an interactive uigetfile window
%   >> EEG = pop_loadhdf5( filename, filepath );                % no pop-up window 
%   >> EEG = pop_loadhdf5( filename, filepath, rejectchannels);
%
% Optional input:
%   filename       - Filename of the *.hdf5 file
%   filepath       - File path of the *-hdf5 file
%   rejectchan     - Channel number to reject (remove) from the data 
%
% Outputs:
%   EEG            - EEGLAB data structure
%
% Author: Simon Lind Kappel, Aarhus University, 2015


function [EEG, command] = pop_loadhdf5(varargin)

command = '';
EEG = [];

%% Create input argument structure
if ~isempty(varargin)
    try 
        g = struct(varargin{:});
    catch
        disp('Setevent: wrong syntax in function arguments'); 
        return; 
    end
else
    g = [];
end;

%% History
command = sprintf('[EEG, command] = pop_loadhdf5(%s);', vararg2str(varargin));

%% Create new EEG struct
EEG = eeg_emptyset;

%% test the presence of variables
if ~isfield(g, 'filename')
    [g.filename, g.filepath, FilterIndex] = uigetfile({'*.hdf5','*.hdf5 (gRecorder)'},'Select a recording from the g-tec amplifier');
    if FilterIndex < 1, return, end
end

if ~isfield(g, 'filepath')
    g.filepath = pwd;
end

if (~strcmp(g.filename(end-4:end),'.hdf5'))
    g.filename = strcat(g.filename,'.hdf5');
end

EEG.filename = g.filename;
EEG.filepath = g.filepath;

Filename = fullfile(g.filepath,g.filename);

%% Read the hdf5 file and extract the EEG data and other options

%Read info about all attributes in the hdf5 file.
fileinfo = hdf5info(Filename);

%Read Raw EEG data and scale from µV to V
[EEG.data, ~] = hdf5read(Filename,'/RawData/Samples');
EEG.data = single(EEG.data);

%Get the XML data in the AcquisitionTaskDescription attribute
[AcquisitionTaskDescription, ~] = hdf5read(Filename,'/RawData/AcquisitionTaskDescription');

%Extract samplerate
FsTemp = ReadXmlAttribute(AcquisitionTaskDescription.Data, 'SampleRate');
EEG.srate = str2double(FsTemp{1});

%Extract Channel number
ChNrTemp = ReadXmlAttribute(AcquisitionTaskDescription.Data, 'PhysicalChannelNumber');
EEG.gtec.chnumber = str2double(ChNrTemp);

%Extract Channel number
DeviceNumberTemp = ReadXmlAttribute(AcquisitionTaskDescription.Data, 'DeviceNumber');
EEG.gtec.devnumber = str2double(DeviceNumberTemp);

%Get names of the EEG channels
ChannelNames = ReadXmlAttribute(AcquisitionTaskDescription.Data, 'ChannelName');

%Find the array index of the '/AsynchronData' item.
AsynchronDataIdx = 0;
for n=1:length(fileinfo.GroupHierarchy.Groups)
    if strcmp(fileinfo.GroupHierarchy.Groups(n).Name,'/AsynchronData')
        AsynchronDataIdx = n;
    end
end

%% Extract the trigger data from the hdf5 file and create the event struct for EEGlab

%The trigger data is stored in '/AsynchronData' [ Groups(1) ]. Check if any data is
%stores in this directory.
if (AsynchronDataIdx > 0)
    if (length(fileinfo.GroupHierarchy.Groups(AsynchronDataIdx).Datasets) > 1)
        
        %Get the XML data in the '/AsynchronData/AsynchronSignalTypes' attribute
        [TriggerChannelInfo, ~] = hdf5read(Filename,'/AsynchronData/AsynchronSignalTypes');
        TriggerName = ReadXmlAttribute(TriggerChannelInfo.Data, 'Name');
        
        %TrigChannelNum = ReadXmlAttribute(TriggerChannelInfo.Data, 'ChannelNumber');
        TrigID = str2double(ReadXmlAttribute(TriggerChannelInfo.Data, 'ID'));
        
        [TrigTime, ~] = hdf5read(Filename,'/AsynchronData/Time');
        [TrigTypeID, ~] = hdf5read(Filename,'/AsynchronData/TypeID');
        
        if length(TrigTime) ~= length(TrigTypeID)
           warning(sprintf('Trigtype is missing for the last %i triggers',length(TrigTime)-length(TrigTypeID)));
           TrigTime = TrigTime(1:length(TrigTypeID));
        end
        
        EEG.event = struct('type',{},'position',[],'latency',[],'urevent',[]); %Columns: 1=type, 2=position, 3=latency, 4=urevent
        [~,TrigSort] = sort(TrigTime);
        
        TrigTime = double(TrigTime);
        
        for n=1:length(TrigTime);
            EEG.event(n).type = TriggerName{TrigID == TrigTypeID(TrigSort(n))};
            EEG.event(n).position = 1;
            EEG.event(n).latency = TrigTime(TrigSort(n));
            EEG.event(n).urevent = n;
        end
    end
end

%% handle input arguments
rejectchans = [];

tmpfields = fieldnames(g);
for curfield = tmpfields'
    switch lower(curfield{1})
        case {'filename' } % do nothing now
        case {'filepath' } % do nothing now
        case 'rejectchan'   
            rejectchans = getfield(g, {1}, curfield{1});            
        otherwise, error(['pop_editset() error: unrecognized field ''' curfield{1} '''']);
    end;
end;

if ~isfield(lower(curfield{1}),'rejectchan' )
    res = inputgui( 'geometry', { [1 1] }, ...
        'geomvert', [1], 'uilist', { ...
        { 'style', 'text', 'string', [ 'Reject channels:' 1 1 1 ] }, ...
        { 'style', 'edit', 'string', ''}});
    if ~isempty(res)
        try
            rejectchans = eval(res{1});
        catch
            warning('The reject channel string was not formatted correctly');
        end
    end
end

if ~isempty(rejectchans) && isempty(find(isnan(rejectchans) == 1,1))
    range_chan = 1:size(EEG.data,1);
    range_chan(rejectchans) = [];
    EEG.data = EEG.data(range_chan,:);
    ChannelNames = ChannelNames(range_chan);
    EEG.gtec.chnumber = EEG.gtec.chnumber(range_chan);
    EEG.gtec.devnumber = EEG.gtec.devnumber(range_chan);
end

%% Asign field in the EEG struct
EEG.chanlocs        = struct('labels', cellstr(ChannelNames));
EEG.nbchan          = size(EEG.data,1);
EEG.pnts            = size(EEG.data,2);
EEG.trials          = 1;
EEG.xmin            = 0;
EEG.xmax            = EEG.xmin + (EEG.pnts-1)/EEG.srate;
EEG.setname 		= EEG.filename;

EEG = eeg_checkset(EEG);

%% Extract data from Xml string
function [AttributeValue] = ReadXmlAttribute(DataString, AttributeName)
AttributeLocated = 0;
AttributeValueIndex = 1;
temp = '';
for i=length(AttributeName)+3:length(DataString)-length(AttributeName)-3
    
    % Start of string has been located
    if strcmp(DataString(i-length(AttributeName)-2:i-1), ['<' AttributeName '>'])
        AttributeLocated = 1;
        temp = '';
    end
    
    % End of string has been located
    if strcmp(DataString(i:i+length(AttributeName)+2), ['</' AttributeName '>'])
        AttributeLocated = 0;
        AttributeValue{AttributeValueIndex} = temp;
        AttributeValueIndex = AttributeValueIndex + 1;
    end
    
    if (AttributeLocated)
        temp(AttributeLocated) = DataString(i);
        AttributeLocated = AttributeLocated +1;
    end  
end

if (~exist('AttributeValue','var'))
    AttributeValue{1} = '';
end
