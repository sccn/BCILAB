function run_writesnap(varargin)
% Send real-time BCI outputs to a SNAP-based stimulus presentation system.
% run_writetcp(Model,SourceStream,OutputHost,OutputPort,OutputForm,MessageFormat,UpdateFrequency,PredictorName,ConnTimeout)
%
% This function runs in the background and processes data from some MATLAB stream (created with some
% other background input plugin, e.g., the BioSemi reader). The processed data is periodically 
% forwarded to a program running the SNAP framework.
%
% In:
%   Model : A model data structure (as obtained from bci_train) based on which the predictor shall be 
%           created; typically this is a model struct, but for convenience it can be a file name, 
%           variable name in the base workspace, or a cell array of {file name, variable name} to 
%           refer to a variable inside a .mat file. The model is not modified by this function.
%           (default: 'lastmodel')
%
%   SourceStreamNames : Optional names of stream data structures in the MATLAB base workspace to
%                       consider as possible data sources (previously created with onl_newstream); 
%                       if a stream contains all channels that are needed by the predictor, or 
%                       alternatively has the right number and type of channels it will be considered 
%                       as a potential source stream unless ambiguous. (default: 'laststream')
%
%   OutputHost : destination host name to send results to (computer name, URL or IP address)
%                (default: 'localhost')
%
%   OutputPort : destination port to which the results are sent (default: 7897)
%
%   OutputForm : output data form, see onl_predict (default: 'distribution')
%
%   TargetVariable : name of the variable in the SNAP module which should carry the BCI predictions
%                    (default: 'bci')
%
%   UpdateFrequency : The rate at which new outputs will be computed. (default: 10)
%
%   PredictAt : Predict at markers. If nonempty, this is a cell array of online target markers relative 
%               to which predictions shall be made. If empty, predictions are always made on the most recently 
%               added sample. (default: {})
%
%   PredictorName : Name of the predictor to be created; a variable of this name will be created in 
%                   the MATLAB base workspace to hold the predictor's state. If a variable with this
%                   name already exists it will be overridden. (default: 'lastpredictor')
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
opts = arg_define(varargin, ...
    arg({'pred_model','Model'}, 'lastmodel', [], 'Predictive model. As obtained via bci_train or the Model Calibration dialog.','type','expression'), ...
    arg({'in_stream','SourceStreamNames','SourceStream'}, 'laststream',[],'Input Matlab stream name(s). Optional names of stream data structures in the MATLAB base workspace to consider as possible data sources (previously created with onl_newstream); if a stream contains all channels that are needed by the predictor, or alternatively has the right number and type of channels it will be considered as a potential source stream unless ambiguous.','typecheck',false,'shapecheck',false), ...
    arg({'out_hostname','OutputHost','Host'}, 'localhost',[],'Destination TCP hostname. Can be a computer name, URL, or IP address.'), ...
    arg({'out_port','OutputPort','Port'}, 7897, [],'Destination TCP port. Matches the SNAP default.'), ...
    arg({'out_form','OutputForm','Form'},'expectation',{'expectation','distribution','mode'},'Output form. Can be the expected value (posterior mean) of the target variable, or the distribution over possible target values (probabilities for each outcome, or parametric distribution), or the mode (most likely value) of the target variable.'), ...
    arg({'target_variable','TargetVariable','Target'},'bci',[],'Target Variable. This is the variable name in the SNAP module that should receive the BCI value'), ...
    arg({'update_freq','UpdateFrequency'},10,[0 Inf],'Update frequency. The rate at which new outputs will be computed.'), ...
    arg({'predict_at','PredictAt'}, {},[],'Predict at markers. If nonempty, this is a cell array of online target markers relative to which predictions shall be made. If empty, predictions are always made on the most recently added sample.','type','expression'), ...
    arg({'pred_name','PredictorName'}, 'lastpredictor',[],'Name of new predictor. A variable of this name will be created in the MATLAB base workspace to hold the predictor''s state. If a variable with this name already exists it will be overridden.'), ...
    arg({'verbose_output','Verbose'}, false,[],'Verbose output. Whether to display verbose outputs (e.g., connection failure).'));

% convert format strings to formatting functions
if ~isvarname(opts.target_variable)
    disp('Note: The given target variable is likely not a valid variable name.'); end

% get a fresh connection id
id = sprintf('c%.0f',fresh_id('snapsocket'));

% start background writer job
onl_write_background( ...
    'ResultWriter',@(y)send_message(y,opts.verbose_output,id,opts.out_hostname,opts.out_port,opts.target_variable),...
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
function send_message(yy,verbose_output,id,host,port,varname)
persistent conns;
try
    % try to connect if necessary
    if ~isfield(conns,id)
        conns.(id) = connect(host,port); end
    % for each prediction...
    for k=1:size(yy,1)
        y = yy(k,:);
        try
            % send it off
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
    end
catch e2
    if strcmp(e2.identifier, 'MATLAB:Java:GenericException')
        if verbose_output
            fprintf('Could not connect to %s:%.0f\n',host,port); end
    else
        hlp_handleerror(e2);
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
