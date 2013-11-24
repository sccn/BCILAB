function y = onl_predict(name,outfmt,suppress_console_output,empty_result_value)
% Query a predictor given the current contents of the stream(s) referenced by it.
% Result = onl_predict(Name,Format,SuppressOutput,EmptyResultValue)
%
% After a predictive model has been loaded successfully into the online system using
% onl_newpredictor, it can be "queried", i.e. its outputs can be requested, at any time and any
% rate, using onl_predict.
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
%   EmptyResultValue : The value to return when no result is available (e.g., not enough data yet)
%                      (default: NaN)
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

    % handle inputs
    if nargin < 4
        empty_result_value = NaN;
        if nargin < 2
            outfmt = 'raw'; end
        if nargin < 3
            suppress_console_output = true; end
    end

    % run do_predict() with the online scope set (expression system disabled, and is_online set to 1)
    try
        if suppress_console_output
            [output,y] = evalc('hlp_scope({''disable_expressions'',1,''is_online'',1},@do_predict)'); %#ok<ASGLU>
        else
            y = hlp_scope({'disable_expressions',1,'is_online',1},@do_predict);
        end
    catch e
        hlp_handleerror(e);
        y = empty_result_value;
    end

    
    function y = do_predict()
        % get predictor from base workspace
        pred = evalin('base',name);
        % get new data from each input pipeline of the prediction function
        for p=length(pred.pipelines):-1:1
            [buffers{p},pred.pipelines{p}] = onl_filtered(pred.pipelines{p},0,false,false); 
            empty(p) = isempty(buffers{p}.data);
        end
        if ~any(empty)
            % invoke the prediction function
            if pred.stateful
                [y,pred] = pred.tracking.prediction_function(struct('streams',{buffers}),pred);
            else
                y = pred.tracking.prediction_function(struct('streams',{buffers}),pred);
            end
            % format the results
            y = utl_formatprediction(y,outfmt);
        else
            y = empty_result_value;
        end
        % write back the updated predictor
        assignin('base',name,pred);
    end
end
