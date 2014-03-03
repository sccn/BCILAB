function run_writeosc(varargin)
% Send real-time BCI outputs via OSC to some remote destination.
% run_writeosc(Model,SourceStream,OutputHost,OutputPort,OutputPath,OutputForm,MessageFormatter,UpdateFrequency,PredictorName,ConnTimeout)
%
% This function runs in the background and processes data from some MATLAB stream (created with some
% other background input plugin, e.g., the BioSemi reader). The processed data is periodically 
% forwarded to some destination via OSC, in a user-specified message format.
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
%   OutputHost : Destination host name to send results to (computer name, URL or IP address)
%                (default: '127.0.0.1')
%
%   OutputPort : Destination port to which the results are sent (default: 12345)
%
%   OutputPath : Output OSC address path (default: '/bci')
%
%   OutputForm : Output data form, see onl_predict (default: 'distribution')
%
%   MessageFormatter : Formatting function to convert each BCI output into an OSC message (which is
%                      a cell array of scalars of certain types) 
%                      (default: '@(D) num2cell(single(D(:))))')
%
%   UpdateFrequency : The rate at which new outputs will be computed. (default: 10)
%
%   PredictorName : Name of the predictor to be created; a variable of this name will be created in 
%                   the MATLAB base workspace to hold the predictor's state. If a variable with this
%                   name already exists it will be overridden. (default: 'lastpredictor')
%
% Examples:
%   % using previously learned predictive model 'mymodel', process data from the input stream 'mystream'
%   % and send the result to OSC destination mymac:22050, under the path /bci
%   run_writeosc('mymodel','mystream','mymac',22050)
%
%   % as before, but pass arguments by name
%   run_writeosc('Model','mymodel','SourceStream','mystream','OutputHost','mymac','OutputPort',22050)
%
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-11-19

% declare the name of this component (shown in the menu)
declare_properties('name','OSC');

% define arguments
opts = arg_define(varargin, ...
    arg({'pred_model','Model'}, 'lastmodel', [], 'Predictive model. As obtained via bci_train or the Model Calibration dialog.','type','expression'), ...
    arg({'in_stream','SourceStreamNames','SourceStream'}, 'laststream',[],'Input Matlab stream name(s). Optional names of stream data structures in the MATLAB base workspace to consider as possible data sources (previously created with onl_newstream); if a stream contains all channels that are needed by the predictor, or alternatively has the right number and type of channels it will be considered as a potential source stream unless ambiguous.'), ...
    arg({'out_hostname','OutputHost'}, '127.0.0.1',[],'Destination TCP hostname. Can be a computer name, URL, or IP address.'), ...
    arg({'out_port','OutputPort'}, 12346, [],'Destination TCP port. Depends on the destination application, usually > 1024.'), ...
    arg({'out_path','OutputPath'},'/bci',[],'Output OSC address. This is the local name of the OSC method that should be invoked (begins with a /)'), ...
    arg({'out_form','OutputForm'},'distribution',{'expectation','distribution','mode'},'Output form. Can be the expected value (posterior mean) of the target variable, or the distribution over possible target values (probabilities for each outcome, or parametric distribution), or the mode (most likely value) of the target variable.'), ...
    arg({'msg_format','MessageFormatter'},'@(D) num2cell(single(D(:))))',[],'Message Formatting function. Converts a BCI output (usually scalar or row vector, in some complex cases a matrix) into a cell array of OSC-compatible data types (float32,int32, ...).','type','expression'), ...
    arg({'update_freq','UpdateFrequency'},10,[0 Inf],'Update frequency. The rate at which new outputs will be computed.'), ...
    arg({'predict_at','PredictAt'}, {},[],'Predict at markers. If nonempty, this is a cell array of online target markers relative to which predictions shall be made. If empty, predictions are always made on the most recently added sample.','type','expression'), ...
    arg({'pred_name','PredictorName'}, 'lastpredictor',[],'Name of new predictor. A variable of this name will be created in the MATLAB base workspace to hold the predictor''s state. If a variable with this name already exists it will be overridden.'), ...
    arg({'verbose_output','Verbose'}, false,[],'Verbose output. Whether to display verbose outputs (e.g., connection failure).'));

% convert format strings to formatting functions
if ischar(opts.msg_format)
    opts.msg_format = @(D) sprintf(opts.msg_format,D); end

% connect to remote address via OSC
if ~exist('osc_new_address','file')
    try
        build_osc;
    catch e
        error('The OSC library has not been built for your platform yet; see dependencies/OSC* for more info.'); 
    end
end
conn = osc_new_address(opts.out_hostname,opts.out_port);
deleter = onCleanup(@()osc_free_address(conn)); % if the last reference to this is dropped, the connection is closed (on MATLAB 2008a+)

% start background writer job
onl_write_background( ...
    'ResultWriter',@(y)send_message(y,conn,opts.out_path,opts.msg_format,deleter),...
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
function send_message(yy,conn,opath,formatter,deleter)
for k=1:size(yy,1)
    y = yy(k,:);
    msg = struct('path',opath,'data',{formatter(y)});
    if osc_send(conn,msg) == 0
        error('OSC transmission failed.'); end
end