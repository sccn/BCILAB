function bci_Initialize(in_signal_dim, out_signal_dim)
% Set processing state according to the bci_Parameters assigned by BCI2000.
% bci_Initialize(In-Signal-Dim,Out-Signal-Dim)
%
% Assumes that bci_Preflight completed successfully.
%
% In:
%  In-Signal-Dim  : a two-element vector of the form [#Channels, #Samples] which describes the input signal's dimensions
%
% Out:
%  Out-Signal-Dim : a two-element vector of the form [#Channels, #Samples] which describes the input signal's dimensions
%
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2010-09-08

% these variables are shared with Initialize and Process
global model meta_data in_samples out_samples next_outputs;

quicklog('bcilab-log.txt','called with in_signal_dim = [%i, %i] and out_signal_dim = [%i, %i]',in_signal_dim(1),in_signal_dim(2),out_signal_dim(1),out_signal_dim(2));

if isempty(meta_data)
    quicklog('bcilab-log.txt','Meta-data field is empty (but was properly set in bci_Preflight - aborting');
    error('Meta-data is empty, apparently bci_Preflight was not called properly...'); 
end

% create the stream
try
    quicklog('bcilab-log.txt','Now creating a new stream...');
    quicklog('bcilab-log.txt','srate will be %.1f',meta_data.srate);
    quicklog('bcilab-log.txt','chanlocs will be %s',hlp_tostring({meta_data.chanlocs.labels}));
    onl_newstream('bci2000_stream',meta_data);
    quicklog('bcilab-log.txt','Success...');
catch le
    quicklog('bcilab-log.txt','Stream creation failed; error: %s',le.message);
    quicklog('bcilab-log.txt','Full meta-data properties: %s',hlp_tostring(meta_data));
    rethrow(le);
end

% load a predictor
try
    quicklog('bcilab-log.txt','Now loading the predictive model...');
    onl_newpredictor('bci2000_model',model,{'bci2000_stream'});
catch le
    quicklog('bcilab-log.txt','Model creation failed; error: %s',le.message);
    rethrow(le);
end
    
% set the consumed/computed sample numbers to zero
in_samples = 0;  % number of samples consumed & emitted so far
out_samples = 0; % number of samples computed so far (computation happens at a usually lower rate)

% this should be set in bci_StartRun
next_outputs = {zeros(length(meta_data.chanlocs),1)};

quicklog('bcilab-log.txt','Now ready for processing; exiting...');
