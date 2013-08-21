function y = onl_predict(name,outfmt,suppress_output)
% Query a predictor given the current contents of the stream(s) referenced by it.
% Result = onl_predict(Name,Format)
%
% After a predictive model has been loaded successfully into the online system (which involves 
% opening and linking it to the necessary data streams), it can be "queried", i.e. its outputs can
% be requested, at any time and any rate, using this function.
%
% In:
%   Name : name of a predictor (under which is was previously created with onl_newpredictor)
%
%   Format : the desired form of the prediction (see also ult_formatprediction), can be:       
%            * 'raw': the raw prediction, as defined by ml_predict (default)
%            * 'expectation': the output is the expected value (i.e., posterior mean) of the
%                             quantity to be predicted; can be multi-dimensional [1xD], but D
%                             equals in most cases 1
%            * 'distribution': the output is the probability distribution (discrete or
%                              continuous) of the quantity to be predicted usually, this is a
%                              discrete distribution - one probability value for every possible
%                              target outcome [1xV] it can also be the parameters of a
%                              parametric distribution (e.g., mean, variance) - yielding one
%                              value for each parameter [DxP]
%            * 'mode': the mode [1xD], or most likely output value (only supported for discrete
%                      probability distributions)
%
%   SuppressOutput : whether to suppress console output (default: true)
%
% Out:
%   Result : Predictions of the selected model(s) w.r.t. to the most recent data.
%
% Example:
%   % obtain a prediction from a previoussly loaded model
%   output = onl_predict('mypredictor')
%
% See also:
%   onl_newpredictor, onl_newstream, onl_append
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-03

if nargin < 2
    outfmt = 'raw'; end
if nargin < 3
    suppress_output = true; end


try
    if ~suppress_output
        % run do_predict() with the regular online scope set (expression system disabled, and is_online set to 1)
        y = hlp_scope({'disable_expressions',1,'is_online',1},@do_predict,name,outfmt);
    else
        % run do_predict() with console output suppressed, and using the regular online scope
        % (expression system disabled, and is_online set to 1)
        [output,y] = evalc('hlp_scope({''disable_expressions'',1,''is_online'',1},@do_predict,name,outfmt)'); %#ok<ASGLU>
    end
catch e
    if ~exist('name','var')
        error('BCILAB:onl_predict:noname','The name of the predictor to use must be specified'); end
    if ~strcmp(e.identifier,'BCILAB:set_makepos:not_enough_data')
        disp('onl_predict() encountered an error; Traceback: ');
        env_handleerror(e); 
    end
    y = NaN;
end



function y = do_predict(name,outfmt)
% This function does the actual work -- depending on the extra options on onl_predict it is called 
% in different ways

% get the predictor from the workspace
try    
    pred = evalin('base',name);
catch
    error(['A predictor with name ' name ' does not exist in the workspace.']);
end

% make sure that the prediction function has the right format
if ischar(pred.tracking.prediction_function)
    % prediction function given as a string
    if strncmp(pred.tracking.prediction_function,'Paradigm',8)
        % class reference: instantiate
        instance = eval(pred.tracking.prediction_function); %#ok<NASGU>
        pred.tracking.prediction_function = eval('@instance.predict');
    else
        % some other function
        pred.tracking.prediction_function = str2func(pred.tracking.prediction_function);
    end
end

% get chunks of processed data for each input signal of the prediction function
for p=length(pred.pipelines):-1:1
    [buffers{p},pred.pipelines{p}] = onl_filtered(pred.pipelines{p}, pred.tracking.prediction_window(p),false,false); end


% check if we have enough data in each chunk
if ~all(pred.tracking.prediction_window==cellfun(@(buf)size(buf.data,2),buffers) | pred.tracking.prediction_window==0)
    % not enough data yet
    y = NaN;
else
    
    % ensure that it has a .stateful field
    if ~isfield(pred,'stateful')
        pred.stateful = is_stateful(pred.tracking.prediction_function,[],[]); end
    % ensure that .arg_direct is set to 1 (this way the prediction functions doesn't start
    % re-parsing its arguments since we assume that the model is specified completely)
    pred.arg_direct = 1;
    
    % invoke the prediction function appropriately
    if pred.stateful
        [y,pred] = pred.tracking.prediction_function(struct('streams',{buffers}),pred);
    else
        y = pred.tracking.prediction_function(struct('streams',{buffers}),pred);
    end
    
    % format the results
    y = utl_formatprediction(y,outfmt);
    
end

% write back the updated predictor
assignin('base',name,pred);
