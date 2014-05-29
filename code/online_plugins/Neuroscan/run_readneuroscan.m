function run_readneuroscan(varargin)
% Receive real-time data from Neuroscan Scan Recorder
% run_readneuroscan(MatlabStream,InputHost,Port)
%
%
% In:
%   MatlabStream : name of the stream to create in the MATLAB environment (default: 'laststream')
%     
%   InputHost: Source TCP hostname. Can be a computer name, URL, or IP
%              address (default: '127.0.0.1')
%   
%   Port : the port on which to connect to the TCP host (default: 4000)
%     
%   UpdateFrequency : update frequency, in Hz (default: 25)
%     
%     Examples:
%     % open an input stream that named 'newstream' on port 5000, using all
%     % other defaults
%     run_readneuroscan('MatlabStream', 'newstream', 'Port', 5000);
%
% Author: Visual Attention and Cognition Lab, Dan Roberts, and Nick Peñaranda, George Mason University, Spring 2014
%         Released under the GPLv3, see COPYING.txt
%         Based on the BrainVision BCILAB plug-in by Hal Greenwald

% declare the name of this component (shown in the menu)
declare_properties('name','Neuroscan Recorder');

% read options
opts = arg_define(varargin, ...
    arg({'new_stream','MatlabStream'}, 'laststream',[],'New Stream to create. This is the name of the stream within the MATLAB environment.'), ...
    arg({'src_hostname','InputHost'}, '127.0.0.1',[],'Source TCP hostname. Can be a computer name, URL, or IP address.'), ...
    arg({'src_port','Port'}, 4000,[],'TCP Host Port.'), ...
    arg({'update_freq','UpdateFrequency'}, 25,[],'Update frequency (Hz). New data is polled at this rate, in Hz.'),...
    arg({'chan_labels', 'ChannelLabels'}, 'standard', [], 'Enter custom channel labels if applicable')...
    );

% open a connection
h = ns_open(opts.src_hostname, opts.src_port);

if strcmpi(opts.chan_labels, 'standard')
    if h.numChan == 40 % assume standard NuAmps 40-channel montage
        h.channelNames = {'HEOL', 'HEOR', 'FP1', 'FP2', 'VEOU', 'VEOL', 'F7', 'F3', 'FZ',...
            'F4', 'F8', 'FT7', 'FC3', 'FCZ', 'FC4', 'FT8', 'T3', 'C3', 'CZ', 'C4', 'T4', ...
            'TP7', 'CP3', 'CPZ', 'CP4', 'TP8', 'A1', 'T5', 'P3', 'PZ', 'P4', 'T6', ...
            'A2', 'O1', 'OZ', 'O2', 'FT9', 'FT10', 'PO1', 'PO2'};
    elseif h.numChan == 68 % assume standard SynAmps2 montage,
        % which includes 64 EEG channels + VEOG, HEOG, EKG, EMG
        h.channelNames = {'FP1','FPZ','FP2','AF3','AF4','F7','F5','F3',...
            'F1','FZ','F2','F4','F6','F8','FT7','FC5','FC3','FC1','FCZ',...
            'FC2','FC4','FC6','FT8','T7','C5','C3','C1','CZ','C2','C4',...
            'C6','T8','M1','TP7','CP5','CP3','CP1','CPZ','CP2','CP4',...
            'CP6','TP8','M2','P7','P5','P3','P1','PZ','P2','P4','P6',...
            'P8','PO7','PO5','PO3','POZ','PO4','PO6','PO8','CB1','O1',...
            'OZ','O2','CB2','VEO','HEO','EKG','EMG'};        
    else
       error('unknown default channel labels for this montage'); 
    end
else % custom channel labels
    opts.chan_labels = evalin('base', opts.chan_labels);
    if length(opts.chan_labels) ~= h.numChan
        errorMsg = sprintf('the number of custom channel labels provided (%i) does not match the number of channels in the data stream (%i)', length(opts.chan_labels), h.numChan);
        ns_close(h);
        error(errorMsg); %#ok<SPERR>
    end
end

h.name = opts.new_stream;

%Create and initialize online stream
onl_newstream(opts.new_stream, 'srate', h.srate, 'chanlocs', h.channelNames, 'data', zeros(length(h.channelNames),0,0),'xmin',toc(uint64(0)));

% start background acquisition
onl_read_background(opts.new_stream,@()ns_read(h), opts.update_freq);