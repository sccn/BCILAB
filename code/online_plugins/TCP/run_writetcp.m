function run_writetcp(varargin)
% Send real-time BCI outputs via TCP to some remote destination.
% run_writetcp(Model,SourceStream,OutputHost,OutputPort,OutputForm,MessageFormat,UpdateFrequency,PredictorName,VerboseOutput)
%
% This function runs in the background and processes data from some MATLAB stream (created with some
% other background input plugin, e.g., the BioSemi reader). The processed data is periodically 
% forwarded to some destination via TCP, in a user-specified message format.
%
% In:
%   Model : predictive model to use (see onl_newpredictor) (default: 'lastmodel')
%
%   SourceStream : real-time stream name to read from (in MATLAB workspace) (default: 'laststream')
%
%   OutputHost : destination host name to send results to (computer name, URL or IP address)
%                (default: '127.0.0.1')
%
%   OutputPort : destination port to which the results are sent (default: 12345)
%
%   OutputForm : output data form, see onl_predict (default: 'distribution')
%
%   MessageFormat : Formatting function or format string to use for each emitted samples 
%                   (default: @(D)[sprintf('%.3f ',D) '/n'], one line per output, in text format)
%                   if this is a string, it needs to be a formatting string as accepted by sprintf,
%                   written in extra quotes. If it is a function, it needs to accept a data item 
%                   (see onl_predict) an return a string.
%
%   UpdateFrequency : update frequency (default: 10)
%
%   PredictAt : Predict at markers. If nonempty, this is a cell array of online target markers relative 
%               to which predictions shall be made. If empty, predictions are always made on the most recently 
%               added sample. (default: {})
%
%   PredictorName : name for new predictor, in the workspace (default: 'lastpredictor')
%
%   Verbose : whether to display verbose outputs (e.g. connection failure) (default: false)

%
% Examples:
%   % open a new BCILAB processing stream, using the previously learned predictive model 'mymodel',
%   % and reading from a previously opened input stream named 'mystream'; send outputs to some IP:port
%   run_writetcp('mymodel','mystream','192.168.1.5',12467)
%
%   % as before, but pass arguments by name
%   run_writetcp('Model','mymodel','SourceStream','mystream','OutputHost','192.168.1.5','OutputPort',12467)
%
%   % as before, but transmit the data in a raw binary encoding
%   run_writetcp('Model','mymodel','SourceStream','mystream','OutputHost','192.168.1.5','OutputPort',12467,'MessageFormat',@(x)typecast(x,'uint8'))
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-11-19

% declare the name of this component (shown in the menu)
declare_properties('name','TCP');

% define arguments
opts = arg_define(varargin, ...
    arg({'pred_model','Model'}, 'lastmodel', [], 'Predictive model. As obtained via bci_train or the Model Calibration dialog.','type','expression'), ...
    arg({'in_stream','SourceStream'}, 'laststream',[],'Input Matlab stream. This is the stream that shall be analyzed and processed.'), ...
    arg({'out_hostname','OutputHost'}, '127.0.0.1',[],'Destination TCP hostname. Can be a computer name, URL, or IP address.'), ...
    arg({'out_port','OutputPort'}, 12346, [],'Destination TCP port. Depends on the destination application, usually > 1024.'), ...
    arg({'out_form','OutputForm'},'distribution',{'expectation','distribution','mode'},'Output form. Can be the expected value (posterior mean) of the target variable, or the distribution over possible target values (probabilities for each outcome, or parametric distribution), or the mode (most likely value) of the target variable.'), ...
    arg({'msg_format','MessageFormat'},'@(D)[sprintf(''%.3f '',D) ''/n'']',[],'Message Format. Either a formatting function (Data to string) or an fprintf-style format string (in quotes).','type','expression'), ...
    arg({'update_freq','UpdateFrequency'},10,[],'Update frequency. This is the rate at which the graphics are updated.'), ...
    arg({'predict_at','PredictAt'}, {},[],'Predict at markers. If nonempty, this is a cell array of online target markers relative to which predictions shall be made. If empty, predictions are always made on the most recently added sample.','type','expression'), ...
    arg({'pred_name','PredictorName'}, 'lastpredictor',[],'Name of new predictor. This is the workspace variable name under which a predictor will be created.'), ...
    arg({'verbose_output','Verbose'}, false,[],'Verbose output. Whether to display verbose outputs (e.g., connection failure).'));

% convert format strings to formatting functions
if ischar(opts.msg_format)
    opts.msg_format = @(D) sprintf(opts.msg_format,D); end

% get a fresh connection id
id = sprintf('c%.0f',fresh_id('tcpsocket'));

% start background writer job
onl_write_background( ...
    'ResultWriter',@(y)send_message(y,opts.verbose_output,id,opts.out_hostname,opts.out_port,opts.msg_format),...
    'MatlabStream',opts.in_stream, ...
    'Model',opts.pred_model, ...
    'OutputFormat',opts.out_form, ...
    'UpdateFrequency',opts.update_freq, ...
    'PredictorName',opts.pred_name, ...
    'PredictAt',opts.predict_at, ...
    'Verbose',opts.verbose_output, ...
    'StartDelay',0,...
    'EmptyResultValue',[]);

disp('Now writing...');


% background message sending function
function send_message(yy,verbose_output,id,host,port,formatter)
persistent conns;
try
    % try to connect if necessary
    if ~isfield(conns,id)
        conns.(id) = connect(host,port); end
    % for each prediction...
    for k=1:size(yy,1)        
        y = yy(k,:);
        % send it off
        try
            strm = conns.(id).strm;
            strm.writeBytes(char(formatter(y)));
            strm.flush();
        catch e1
            if strcmp(e1.identifier, 'MATLAB:Java:GenericException')
                % failed to send: try to re-connect...
                conns.(id) = connect(host,port);
            else
                rethrow(e1);
            end
        end
    end
catch e2
    if strcmp(e2.identifier, 'MATLAB:Java:GenericException')
        if verbose_output
            fprintf('Could not connect to %s:%.0f\n',host,port); end
    else
        env_handleerror(e2);
    end
end

function newconn = connect(host,port)
import java.io.*
import java.net.*
import java.lang.*
newconn.sock = Socket();
newconn.sock.setTcpNoDelay(1);
newconn.sock.setSoTimeout(5);
newconn.sock.connect(InetSocketAddress(host,port),1);
newconn.strm = DataOutputStream(newconn.sock.getOutputStream());
disp(['Connection to ' host ' established.']);
