function run_pipetcpip(varargin)
% Run BCILAB's real-time engine as a processing node via TCP/IP
% run_pipetcpip(Arguments...)
%
% This plugin reads raw data from some source host/port and writes predictions to some output
% host/port. This example code is implemented using MATLAB's Instrument Control Toolbox. 
%
% The example plugin reads raw signal data from a hypothetical TCP/IP server, formatted as one
% sample per line (and each line containing a space-separated list of floating-point values, in
% string format). For simplicity, the sampling rate and channel names of the input data stream are
% not communicated over the TCP link, and are therefore to be manually specified in user parameters.
% The plugin processes the incoming data and produces outputs at a fixed rate, which are forwarded
% to some destination TCP server. The output sampling rate is usually much lower than the raw-data
% rate (e.g. matching the screen refresh rate), and the output data form can be specified in an
% optional parameter. The message format is one line per output, as space-separated list of
% floating-point values.
%
%
% In:
%   InputHost : hostname (or IP string) of the source machine (default: '127.0.0.1')
%
%   InputPort : listening port of the source service (default: 12345)
%
%   InputMetadata : Meta-data of the input stream. This is a struct or cell array of name-value 
%                   pairs with meta-data fields to use. The mandatory fields are 'srate' and
%                   'chanlocs', where chanlocs is either a channel locations struct array, or a cell
%                   array of channel names, or the number of channels (in which case a cell array of
%                   the form {'A1','A2', ..., 'A32','B1', ...} is created). Optionally, the field
%                   'datasource' can be set to point to a dataset on disk or in a MATLAB workspace
%                   variable.
%
%   Model : a file, struct or workspace variable name that contains a predictive model, as 
%           previously computed by bci_train (default: 'lastmodel')
%
%   OutputHost : hostname (or IP string) of the destination machine (default: '127.0.0.1')
%
%   OutputPort : listening port of the destination service (default: 12346)
%
%   OutputSamplingRate : rate, in Hz, at which the output stream shall be sampled (default: 20)
%                        make sure that this rate is low enough so that BCILAB can process in real 
%                        time (otherwise it would lag)
%
%   OutputForm : format of the data sent to the output stream, can be one of the following:
%                'expectation': the expected value (= posterior mean) of the outputs; can be 
%                               multi-dimensional but is usually 1d (default) this mode is
%                               appropriate for simple applications that expect a smooth control
%                               signal, or for applications that expect a regression output
%                'distribution' : parameters of the output distribution; for discrete 
%                                 distributions, this is one probability value for each target
%                                 (adding up to 1) this mode is appropriate for more advanced
%                                 applications that use the full output distribution (e.g., for
%                                 decision-theoretical processing) (default)
%                'mode' : the most likely output value (currently only supported for discrete 
%                         distributions) this mode is appropriate for simple applications that
%                         take a non-probabilistic classifier decision (e.g., as from a Support
%                         Vector Machine)
%   
% Example:
%   run_pipetcpip('Model','lastmodel', 'InputPort',2050, 'OuputHost','192.168.1.10','OutputPort',...
%       2051, 'InputMetadata',{'srate',256,'chanlocs',{'C3','Cz','C4'}})
%
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2011-01-18

% declare the name of this component (shown in the menu)
declare_properties('name','TCP (Instrument Control Toolbox)');

% define arguments...
arg_define(varargin, ...
    arg({'in_hostname','InputHost'}, '127.0.0.1', [],'Source TCP hostname. Can be a computer name, URL, or IP address.'), ...
    arg({'in_port','InputPort'}, 12345, [],'Source TCP port. Depends on the source application, usually > 1024.'), ...
    arg_sub({'in_metadata','InputMetadata'},{},@utl_parse_metadata, 'Meta-data of the input stream. These are fields as they appear in EEGLAB data sets; only sampling rate and channel labels are mandatory.'), ...
    arg({'pred_model','Model'}, 'lastmodel', [], 'Predictive model. As obtained via bci_train or the Model Calibration dialog.','type','expression'), ...
    arg({'out_hostname','OutputHost'}, '127.0.0.1',[],'Destination TCP hostname. Can be a computer name, URL, or IP address.'), ...
    arg({'out_port','OutputPort'}, 12346, [],'Destination TCP port. Depends on the destination application, usually > 1024.'), ...
    arg({'out_srate','OutputSamplingRate'}, 20,[],'Output sampling rate. This is the rate at which estimate should be computed. If this value is too high, the BCI will start to lag behind.'), ...
    arg({'out_format','OutputForm'}, 'distribution',{'expectation','distribution','mode'},'Form of the produced output values. Can be the expected value (posterior mean) of the target variable, or the distribution over possible target values (probabilities for each outcome, or parametric distribution), or the mode (most likely value) of the target variable.'));

% parse meta-data specification
meta = utl_parse_metadata(in_metadata);

% set up BCILAB stream
onl_newstream('stream_tcpip','srate',meta.srate,'chanlocs',meta.chanlocs);

% load the given predictor
onl_newpredictor('predictor_tcpip',pred_model,'stream_tcpip');

% connect to source
src = tcpip(in_hostname,in_port);
fopen(src);

% connect to destination
dst = tcpip(out_hostname,out_port);
fopen(out);

t = 0; % current stream clock
while 1
    % get a sample and append it to the BCILAB stream
    sample = str2num(fgetl(src))';
    onl_append('stream_tcpip',sample);
    t = t + 1/meta.srate;
    % if it is time to produce an output sample...
    if t > 1/out_srate
        t = t - 1/out_srate;
        % compute the output
        result = onl_predict('predictor_tcpip',out_format);
        % and send a message
        fprintf(dst,'%s/n',sprintf('%.3f ',result));
    end
end
