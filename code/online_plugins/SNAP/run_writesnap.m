function run_writesnap(varargin)
% Send real-time BCI outputs to a SNAP-based stimulus presentation system.
% run_writetcp(Model,SourceStream,OutputHost,OutputPort,OutputForm,MessageFormat,UpdateFrequency,PredictorName,ConnTimeout)
%
% This function runs in the background and processes data from some MATLAB stream (created with some
% other background input plugin, e.g., the BioSemi reader). The processed data is periodically 
% forwarded to a program running the SNAP framework.
%
% In:
%   Model : predictive model to use (see onl_newpredictor) (default: 'lastmodel')
%
%   SourceStream : real-time stream name to read from (in MATLAB workspace) (default: 'laststream')
%
%   OutputHost : destination host name to send results to (computer name, URL or IP address)
%                (default: 'localhost')
%
%   OutputPort : destination port to which the results are sent (default: 12345)
%
%   OutputForm : output data form, see onl_predict (default: 'distribution')
%
%   TargetVariable : name of the variable in the SNAP module which should carry the BCI predictions
%                    (default: 'bci')
%
%   UpdateFrequency : update frequency (default: 10)
%
%   PredictorName : name for new predictor, in the workspace (default: 'lastpredictor')
%
%   Verbose : whether to display verbose outputs (e.g. connection failure) (default: false)
%
% Examples:
%   % open a new BCILAB processing stream, using the previously learned predictive model 'mymodel',
%   % and reading from a previously opened input stream named 'mystream'; send outputs to some IP:port
%   run_writesnap('mymodel','mystream','192.168.1.5')
%
%   % as before, but pass arguments by name
%   run_writetcp('Model','mymodel','SourceStream','mystream','OutputHost','192.168.1.5')
%
%   % as before, but transmit the data in a raw binary encoding
%   run_writetcp('Model','mymodel','SourceStream','mystream','OutputHost','192.168.1.5','TargetVariable','alpha')
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-11-19

% declare the name of this component (shown in the menu)
declare_properties('name','SNAP');

% define arguments
arg_define(varargin, ...
    arg({'pred_model','Model'}, 'lastmodel', [], 'Predictive model. As obtained via bci_train or the Model Calibration dialog.','type','expression'), ...
    arg({'in_stream','SourceStream'}, 'laststream',[],'Input Matlab stream. This is the stream that shall be analyzed and processed.'), ...
    arg({'out_hostname','OutputHost','Host'}, 'localhost',[],'Destination TCP hostname. Can be a computer name, URL, or IP address.'), ...
    arg({'out_port','OutputPort','Port'}, 7897, [],'Destination TCP port. Matches the SNAP default.'), ...
    arg({'out_form','OutputForm','Form'},'expectation',{'expectation','distribution','mode'},'Output form. Can be the expected value (posterior mean) of the target variable, or the distribution over possible target values (probabilities for each outcome, or parametric distribution), or the mode (most likely value) of the target variable.'), ...
    arg({'target_variable','TargetVariable','Target'},'bci',[],'Target Variable. This is the variable name in the SNAP module that should receive the BCI value'), ...
    arg({'update_freq','UpdateFrequency'},10,[],'Update frequency. This is the rate at which the output is updated.'), ...
    arg({'pred_name','PredictorName'}, 'lastpredictor',[],'Name of new predictor. This is the workspace variable name under which a predictor will be created.'), ...
    arg({'verbose_output','Verbose'}, false,[],'Verbose output. Whether to display verbose outputs (e.g., connection failure).'));

% convert format strings to formatting functions
if ~isvarname(target_variable)
    disp('Note: The given target variable is likely not a valid variable name.'); end

% get a fresh connection id
id = sprintf('c%.0f',fresh_id('snapsocket'));

% start background writer job
onl_write_background(@(y)send_message(y,verbose_output,id,out_hostname,out_port,target_variable),in_stream,pred_model,out_form,update_freq,0,pred_name);

% background message sending function
function send_message(y,verbose_output,id,host,port,varname)
persistent conns;
try
    % try to connect if necessary
    if ~isfield(conns,id)
        conns.(id) = connect(host,port); end
    try
        strm = conns.(id).strm;
        if isscalar(y)
            strm.writeBytes(char([sprintf('setup %s=%.5f',varname,y) 10]));
        elseif ~isempty(y)
            strm.writeBytes(char(['setup ' varname '=(' sprintf('%.5f,',y) ')' 10]));
        end
        strm.flush();
    catch e1
        if strcmp(e1.identifier, 'MATLAB:Java:GenericException')
            % failed to send: try to re-connect...
            conns.(id) = connect(host,port);
        else
            rethrow(e1);
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
