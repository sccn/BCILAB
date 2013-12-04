function run_writelsl(varargin)
% Provide real-time BCI outputs to the lab streaming layer.
% run_writetcp(Model,SourceStream)
%
% This function runs in the background and processes data from some MATLAB stream (created with some
% other background input plugin, e.g., the BioSemi reader). The resulting estimates are offered to 
% the lab streaming layer for use by other programs.
%
% In:
%   Model : predictive model to use (see onl_newpredictor) (default: 'lastmodel')
%
%   SourceStream : real-time stream name to read from in MATLAB workspace (default: 'laststream')
%
%   LabStreamName : name of the stream in the lab streaming layer (default: 'BCI')
%
%   ChannelNames : labels for each output channel (if this is for classification and the output form 
%                  is set to distribution, the channels correspond to the labels of the respective classes 
%                  (default: {'class1','class2', ..., 'classN'} depending on the # of classes)
%
%   OutputForm : output data form, see onl_predict (default: 'distribution')
%
%   UpdateFrequency : update frequency (default: 10)
%
%   PredictAt : Predict at markers. If nonempty, this is a cell array of online target markers relative 
%               to which predictions shall be made. If empty, predictions are always made on the most 
%               recently added sample. (default: {})
%
%   Verbose : Verbose output. If false, the console output of the online pipeline will be suppressed.
%             (default: false)
%
%   PredictorName : name for new predictor, in the workspace (default: 'lastpredictor')
%
% Notes:
%   This code is currently untested.
%
% Examples:
%   % open a new BCILAB processing stream, using the previously learned predictive model 'mymodel',
%   % and reading from a previously opened input stream named 'mystream'. Name the stream 'BCI' and 
%   % name the channels according to the classes whose probabilities they carry (here for some hypothetical
%   % error perception classifier)
%   run_writelsl('mymodel','mystream','BCI',{'error','no_error'})
%
%   % as before, but pass arguments by name
%   run_writelsl('Model','mymodel','SourceStream','mystream','ChannelNames',{'error','no_error'})
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2012-03-21

persistent lib;

% declare the name of this component (shown in the menu)
declare_properties('name','Lab streaming layer');

% define arguments
opts = arg_define(varargin, ...
    arg({'pred_model','Model'}, 'lastmodel', [], 'Predictive model. As obtained via bci_train or the Model Calibration dialog.','type','expression'), ...
    arg({'in_stream','SourceStream'}, 'laststream',[],'Input Matlab stream. This is the stream that shall be analyzed and processed.'), ...
    arg({'out_stream','LabStreamName','Target'},'bci',[],'Name of the lab stream. This is the name under which the stream is provided to the lab streaming layer.'), ...
    arg({'channel_names','ChannelNames'},{'class1','class2'},[],'Output channel labels. These are the labels of the stream''s channels. In a typical classification setting each channel carries the probability for one of the possible classes.'), ...
    arg({'out_form','OutputForm','Form'},'expectation',{'expectation','distribution','mode'},'Output form. Can be the expected value (posterior mean) of the target variable, or the distribution over possible target values (probabilities for each outcome, or parametric distribution), or the mode (most likely value) of the target variable.'), ...
    arg({'update_freq','UpdateFrequency'},10,[],'Update frequency. This is the rate at which the output is updated.'), ...
    arg({'predict_at','PredictAt'}, {},[],'Predict at markers. If nonempty, this is a cell array of online target markers relative to which predictions shall be made. If empty, predictions are always made on the most recently added sample.','type','expression'), ...
    arg({'verbose','Verbose'}, false,[],'Verbose output. If false, the console output of the online pipeline will be suppressed.'), ...
    arg({'source_id','SourceID'}, 'input_data',{'input_data','model'},'Use as source ID. This is the data that determines the source ID of the stream (if the stream is restarted, readers will continue reading from it if it has the same source ID). Can be input_data (use a hash of dataset ID + target markers used for training) or model (use all model parameters).'), ...
    arg({'pred_name','PredictorName'}, 'lastpredictor',[],'Name of new predictor. This is the workspace variable name under which a predictor will be created.'));

% load the model
model = utl_loadmodel(opts.pred_model);

% check if channel labels make sense for the model
if strcmp(opts.out_form,'distribution')
    if isfield(model,'classes') && ~isempty(model.classes)
        if length(opts.channel_names) ~= length(model.classes)
            disp('The number of classes provided by the model does not match the number of provided channel names; falling back to default names.');
            opts.channel_names = cellfun(@(k)['class' num2str(k)],num2cell(1:length(model.classes),1),'UniformOutput',false);
        end
    end
else
    if isfield(model,'classes') && ~isempty(model.classes)
        if length(opts.channel_names) ~= 1
            disp('A classification model will produce just one channel if the output is not in distribution form, but a different number of channels was given. Falling back to the default channel label.');
            opts.channel_names = {'class'};
        end        
    end
end

% try to calculate a UID for the stream
try
    if strcmp(opts.source_id,'model')
        uid = hlp_cryptohash({rmfield(model,'timestamp'),predict_at,in_stream,out_stream});
    elseif strcmp(opts.source_id,'input_data')
        uid = hlp_cryptohash({model.source_data,predict_at,in_stream,out_stream});
    else
        error('Unsupported SourceID option: %s',hlp_tostring(opts.source_id));
    end
catch e
    disp('Could not generate a unique ID for the predictive model; the BCI stream will not be recovered automatically after the provider system had a crash.');
    hlp_handleerror(e);
    uid = '';
end

% instantiate the library
disp('Opening the library...');
if isempty(lib)
    lib = lsl_loadlib(env_translatepath('dependencies:/liblsl-Matlab/bin')); end

% describe the stream
disp('Creating a new streaminfo...');
info = lsl_streaminfo(lib,opts.out_stream,'MentalState',length(opts.channel_names),opts.update_freq,'cf_float32',uid);
% ... including some meta-data
desc = info.desc();
channels = desc.append_child('channels');
for c=1:length(opts.channel_names)
    newchn = channels.append_child('channel');
    newchn.append_child_value('name',opts.channel_names{c});
    newchn.append_child_value('type',opts.out_form);
end

% create an outlet
outlet = lsl_outlet(info);

% start background writer job
onl_write_background( ...
    'ResultWriter',@(y)send_samples(outlet,y),...
    'MatlabStream',opts.in_stream, ...
    'Model',opts.pred_model, ...
    'OutputFormat',opts.out_form, ...
    'UpdateFrequency',opts.update_freq, ...
    'PredictorName',opts.pred_name, ...
    'PredictAt',opts.predict_at, ...
    'Verbose',opts.verbose, ...
    'StartDelay',0,...
    'EmptyResultValue',[]);

disp('Now writing...');


function send_samples(outlet,y)
if ~isempty(y)
    outlet.push_chunk(y'); end

