function prediction = utl_formatprediction(prediction,out_format)
% Utility for post-processing the outputs of onl_predict into some form.
% Formatted-Prediction = utl_formatprediction(Raw-Prediction,Format)
%
% In:
%   Prediction : a prediction as made by ml_predict / onl_predict / bci_predict (for one or more target values)
%
%   Format     : format of the output prediction (in the descriptions, N is the number of predictions); can be one of:
%                 * 'expectation': the output is the expected value (i.e., posterior mean) of the quantity to be predicted; can be multi-dimensional [NxD]
%                 * 'distribution': the output is the probability distribution (discrete or continuous) of the quantity to be predicted 
%                                   usually, this is a discrete distribution - one probability value for every possible target outcome [NxV]
%                                   it can also be the parameters of a parametric distribution (e.g., mean, variance) - yielding one value for each parameter [NxP]
%                 * 'mode': the mode [Nx1], or most likely output value (only supported for discrete probability distributions)
%                 * 'raw': the raw prediction, as defined by ml_predict
% 
% Out:
%   Formatted-Prediction : The formatted output, as described in Format
%
% See also:
%   ml_predict
%
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2010-09-08

% raw format: nothing to do
if strcmp(out_format,'raw')
    return; end

% we have valid outputs...
if iscell(prediction)
    switch out_format
        case 'expectation'
            % expected value is requested
            if strcmp(prediction{1},'disc')
                % discrete distribution
                prediction = prediction{2}*prediction{3};
            elseif ischar(prediction{1}) && ~strcmp(prediction{1},'struct')
                % standard distribution (e.g. Normal, Binomial, ...)
                args = cell(1,size(prediction{2},2));
                for k = 1:size(prediction{2},2)
                    args{k} = prediction{2}(:,k); end
                prediction = feval([prediction{1} 'stat'],args{:});
            else
                % either a user-specified distribution or a structured prediction
                error('Unsupported predicted distribution');
            end
        case 'mode'
            % the mode of the distribution is requested
            if strcmp(prediction{1},'disc')
                % discrete distribution: most likely class
                prediction = prediction{3}(argmax(prediction{2}'));
            else
                % anything else: we have no way to compute this
                error('Unsupported predicted distribution');
            end
        case 'distribution'
            % we just return the distribution's parameters; the type of the distribution is meta-data
            prediction = prediction{2};
        otherwise
            error('Unknown output format');
    end
elseif isnumeric(prediction)
    % we have a point estimate; there is nothing to be transformed anyway...
    if ~isnan(prediction) && strcmp(out_format,'distribution')
        error('A distribution is not available from this predictor'); end
else
    error('Unknown prediction output');
end

prediction = double(prediction);
