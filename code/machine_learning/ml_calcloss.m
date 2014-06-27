function [measure,stats] = ml_calcloss(type,T,P)
% Calculate the loss for a set of predictions, given knowledge about the target values.
% [Loss,Stats] = ml_calcloss(Type, Target, Prediction)
%
% The performance of predictive models on data [1], or the support of hypotheses given data can be estimated with help of this function, which compares
% a set of target values with predicted values (which are supposed to match the targets) and computes the 'loss' (or cost, or degree of failure). 
% Since both targets and predictions can assume categorical, continuous or multivariate values, or probability distributions over such values, 
% there is a variety of different loss functions to chose from. Depending on the data types of the targets and predictions, a reasonable choice
% of default loss function can often be auto-selected by the function (if Type is specified as []), for example mis-classification rate for categorical
% variables, mean square error for continuous variables, etc.. 
%
% Ultimately, the loss function should reflect most closely the actual loss incurred by a wrong prediction, in the context where the method is supposed to
% be used. For example, if a predictive model that shall predict either +1 or -1 can assign (gradual) probabilities to the two possible outcomes, the most 
% informative loss measure would be the misfit between the actual and predicted probabiltities of those outcomes. If an application maps the output to
% a one-dimensional control signal (e.g. turning a robot gradually left or right), this loss measure would be appropriate (or mean-square error could also 
% be used). If the application, however, makes an exact binary decision based on the probabilities (e.g. playing a sound in the +1 case and not playing
% it otherwise), then the most true-to-the-fact loss function would be the mis-classification rate, and consequently, a model optimized under this 
% loss can be assumed to work better in practice, even though the loss does not capture all information present in the target & prediction values. If the
% relative cost of wrongfully chosing one class (e.g., play sound) is different from that of wrongfully chosing the other class (play no sound), and is 
% not fixed/known in advance, an even better loss would be the (negative) area under the ROC (Receiver-Operating Characteristic) curve.
% 
% In:
%   Type       : the type of loss to compute:
%                * 'auto' / empty: use 'mcr','mse','nll','kld', depending on supplied target & prediction data formats
%                   * misclassification rate is used for discrete probability distributions,
%                   * mean squared error is used for regression outputs
%                   * negative log-likelihood is used for continuous probability distributions
%                   * kullback-leibler divergence is used when both targets and predictions are distributions
%                * 'kld': Kullback-Leibler divergence D_KL(Target||Prediction)
%                * 'nll': negative log-likelihood of target under prediction or vice versa (depending on what is specified as a distribution)
%                * 'mcr'/'err': misclassification rate, for discrete outcomes
%                * 'mae'/'l1': mean absolute error
%                * 'mse'/'l2': mean square error
%                * 'smse': standardized mean square error
%                * 'max'/'linf': maximum absolute error
%                * 'sign': mis-classification rate on the sign of continuous outcomes
%                * 'rms': root mean square error
%                * 'bias': directed mean bias
%                * 'medse': median square error
%                * 'auc': negative area under ROC, for 2-class outcomes
%                * 'cond_entropy': conditional entropy of Target given Prediction, for discrete outcomes
%                * 'cross_entropy': cross-entropy of Target under Prediction
%                * 'f_measure': negative F-score, harmonic mean between precision and recall, for 2-class outcomes
%                * cell array of strings: compute loss for each specified type and return an array of results
%
%   Target     : target variable; can be in any of the formats produced by ml_predict.
%
%   Prediction : predicted variable; can be in any of the formats produced by ml_predict.
%                format is expected to be consistent with (though not necessarily identical to) the Target format, 
%                and the number of samples needs to match that of the target variable
%                
% Out:
%   Loss       : the scalar loss, aggregated over all samples, or a vector of loss measures, if a cell array was specified for the type
%   Stats      : struct of additional statistics (depending on type), or a struct array if a cell array was specified for the type
%
% Note:
%   Simple of the implemented losses have not been tested for all applicable combinations of target/prediction formats. For standard applications,
%   there should be no errors, however. When new losses are implemented, they should be added to this file.
%
% References:
%   [1] MacKay, D. J. C. "Information theory, inference, and learning algorithms." 
%       Cambridge University Press, 2003.
%
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2010-04-22



% the default metric; if type is empty,
max_nll_loss = 1000; % to prevent negative log-likelihoods from blowing up

if iscell(type) && all(cellfun(@ischar,type))
    % call recursively for each type and aggregate statistics
    for t=1:length(type)
        [measure(t),stats{t}] = ml_calcloss(type{t},T,P); end
    measure = measure(1);
    stats = hlp_aggregatestructs(stats,'replace');
else
    % add some stats
    if is_discrete(T)
        stats.classes = T{3};
        stats.class_ratio = mean(T{2});
    elseif is_discrete(P) && is_point(T)
        stats.classes = P{3};
        for c=1:length(stats.classes)
            stats.class_ratio(c) = mean(T==stats.classes(c)); end
    end
    
    % auto-determine measure type
    if isempty(type) || strcmp(type,'auto')
        if is_1d_regression(T) && is_1d_regression(P)
            % we have regression outputs in both cases -> mean squared error
            type = 'mse';
        elseif is_1d_regression(T) && is_discrete(P)
            % the targets are regression values, the estimates are discrete posteriors -> misclassification rate
            type = 'mcr';
        elseif (is_1d_regression(T) && is_1d_distributed(P)) ...
                || (is_1d_distributed(T) && is_1d_regression(P))
            % one is a distribution, the other is a point estimate -> negative log-likelihood
            type = 'nll';
        elseif is_1d_distributed(T) && is_1d_distributed(P)
            % both are 1d distributions (hopefully of compatible types -> kullback-leibler divergence)
            type = 'kld';
        elseif is_discrete(T) && is_discrete(P)
            % both are discrete distributions -> kullback-leibler divergence
            type = 'kld';
        elseif is_structured(T) && is_structured(P)
            % both are structured distributions -> misclassification rate
            type = 'mcr';
        else
            error('an appropriate measure cannot be auto-determined from the data; please specify.');
        end
    end
    
    switch lower(type)
        % compute distribution measures
        case 'kld'
            % kullback-leibler divergence D_KL(targ||pred)
            if is_discrete(T) && is_discrete(P)
                measure = mean(kl_divergence_discrete(T,P));
            elseif is_discrete(P) && is_point(T) && length(unique(T)) <= length(P{3})
                % convert T into a discrete probability distribution
                tmp = zeros(max(T),length(T)); tmp(T+max(T).*(0:(length(T)-1))') = 1; 
                T = {'disc' tmp' P{3}};
                measure = mean(kl_divergence_discrete(T,P));
            elseif is_gaussian(T) && is_gaussian(P)
                measure = mean(kl_divergence_gaussian(T,P));
            elseif is_1d_distributed(T) && is_1d_distributed(P)
                % use the cdf method
                error('not yet implemented');
            else
                error('cannot compute the kullback-leibler divergence for the given training & prediction data format.');
            end
        case 'nll'
            % negative log-likelihood
            if is_point(T) && is_distributed(P)
                measure = sum(neg_loglike(T,P,max_nll_loss));
            elseif is_distributed(T) && is_point(P)
                measure = sum(neg_loglike(P,T,max_nll_loss));
            else
                error('cannot compute log-likelihood for the given training and prediction data format.');
            end
        case {'mcr','err'}
            % misclassification rate
            if is_structured(T) && is_structured(P)
                measure = mean(cellfun(@isequal,T,P));
            else
                % we can enumerate classes & compute a per-class error
                if isnumeric(T) && isnumeric(P) && isequal(size(T),size(P)) && size(T,2) == 1
                    classes = unique([T;P]);
                elseif (is_1d_regression(T) || is_discrete(T)) && (is_1d_regression(P) || is_discrete(P))
                    if is_discrete(T)
                        classes = T{3};
                        T = T{3}(hlp_getresult(2,@max,T{2}')');
                    end
                    if is_discrete(P)
                        classes = P{3};
                        P = P{3}(hlp_getresult(2,@max,P{2}')');
                    end
                else
                    error('target and/or prediction formats not suited for classification');
                end
                try
                    corrects = T == P;
                catch
                    disp('Data unaligned; working around this; this is a severe bug.');
                    len = max(size(T,1),size(P,1));
                    T = T(round((1:len)*size(T,1)/len));
                    P = P(round((1:len)*size(P,1)/len));
                    corrects = T == P;
                end
                    
                measure = 1-mean(corrects);
                % also add per-class error
                stats.classes = classes;
                stats.per_class = nan(size(classes,1),1);
                for c=1:length(classes)
                    if any(T==classes(c))
                        stats.per_class(c) = 1-mean(corrects(T==c)); end
                end
                if length(classes) == 2
                    stats.TP = sum(T==classes(2) & P==classes(2))/sum(T==classes(2));
                    stats.TN = sum(T==classes(1) & P==classes(1))/sum(T==classes(1));
                    stats.FP = sum(T==classes(1) & P==classes(2))/sum(T==classes(1)); % fixed FP/FN rates, reported by Matt Peterson.
                    stats.FN = sum(T==classes(2) & P==classes(1))/sum(T==classes(2));
                end
            end
            
        otherwise
            % compute sample-based measures
            Tx = expected_value(T);
            Px = expected_value(P);
            switch lower(type)
                case {'mae','l1'}
                    measure = mean(abs(Px-Tx));
                case {'mse','l2'}
                    measure = mean((Px-Tx).^2);
                case 'sign'
                    measure = mean(sign(Px) ~= sign(Tx));
                case {'smse'}
                    measure = mean((Px-Tx).^2) ./ var(Tx);
                case {'max','linf'}
                    measure = max(abs(Px-Tx));
                case 'rms'
                    measure = sqrt(mean((Px-Tx).^2));
                case 'bias'
                    measure = mean(Px-Tx);
                case 'medse'
                    measure = median((Px-Tx).^2);
                case 'auc'
                    % derive binary T, continuous P, and the class identities
                    if is_discrete(P)
                        classes = P{3};
                        P = P{2}(:,2);
                        if is_discrete(T)
                            T = T{2}(:,2);
                        elseif ~is_1d_regression(T)
                            error('AUC cannot be computed for non-class targets');
                        end
                    elseif is_discrete(T) && is_1d_regression(P)
                        classes = T{3};
                        T = T{2}(:,2);
                    elseif is_1d_regression(T) && is_1d_regression(P)
                        classes = unique(T);
                    end
                    % compute the AUC
                    if length(classes) ~= 2
                        error('AUC can only be computed for 2-class outcomes'); end
                    if length(unique(T)) > 2
                        error('AUC can only be computed for nonprobabilistic targets.'); end
                    T = T==classes(2);
                    nT = sum(T);
                    ranks = tied_ranks(P);
                    measure = -(sum(ranks(T)) - (nT^2 + nT)/2) / (nT*(length(T)-nT));
                case 'f_measure'
                    classes = [];
                    if is_discrete(T)
                        classes = T{3};
                        T = T{2}*T{3};
                    elseif is_discrete(P) 
                        classes = P{3};
                        P = P{2}*P{3};
                    elseif is_point(T)
                        classes = unique(T);
                    end
                    if length(classes) ~= 2
                        error('computation of the f measure requires 2-class outcomes'); end
                    T = T==classes(2);
                    P = P==classes(2);
                    stats.TP = sum((P==1)&(T==1));
                    stats.TN = sum((P==0)&(T==0));
                    stats.FP = sum((P==1)&(T==0));
                    stats.FN = sum((P==0)&(T==1));
                    measure = -2*stats.TP/(2*stats.TP+stats.FP+stats.FN);
                case 'cross_entropy'
                    if is_discrete(P) && is_1d_regression(T)
                        % convert T into a discrete distribution
                        T = {'disc' repmat(T,1,length(P{3})) == repmat(P{3}',size(T,1),1) P{3}};
                    elseif is_discrete(T) && is_1d_regression(P)
                        % convert P into a discrete distribution
                        if length(unique(P)) > length(T{3})
                            error('P must contain class predictions, not gradual predictions'); end
                        P = {'disc' repmat(P,1,length(T{3})) == repmat(T{3}',size(P,1),1) T{3}};
                    elseif ~(is_discrete(T) || is_discrete(P))
                        error('cross-entropy can currently only be computed for class outcomes.');
                    end
                    measure = mean(kl_divergence_discrete(T,P));
                case 'cond_entropy'
                    error('not yet implemented.');
                otherwise
                    error('unknown measure specified.');
            end
    end
end
stats.measure = type;
stats.(type) = measure;
end



% expectation of distribution D
function V = expected_value(D)
if ~is_distributed(D)
    V = D;
elseif is_discrete(D)
    V = D{2}*D{3};
elseif is_1d_distributed(D)
    params = mat2cell(D{2}, size(D{2},1), ones(1,size(D{2},2)));
    V = feval([D{1},'stat'], params{:});
elseif is_Nd_distributed(D)
    means = cellfun(@(x)x{1}',D{2},'UniformOutput',false);
    V = vertcat(means{:});
else
    error('unsupported distribution format; cannot compute expected value');
end
end



% negative log-likelihood of X under D
function V = neg_loglike(X,D,maxval)
if ~exist('maxval','var')
    maxval = Inf; end
if is_discrete(D)
    match = repmat(X,1,size(D{3},1)) == repmat(D{3}',size(X,1),1);
    V = -log(D{2}(match));
elseif is_1d_distributed(D)
    params = mat2cell(D{2}, size(D{2},1), ones(1,size(D{2},2)));
    V = -log(feval([D{1},'pdf'],X,params{:}));
elseif is_Nd_distributed(D)
    V = zeros(size(X,1),1);
    for k=1:size(X,1)
        V(k) = mvnpdf(X(k,:),D{2}{k}{1},D{2}{k}{2}); end
else
    error('the second argument must be a probability distribution');
end
V(isinf(V)) = maxval;
end




% compute the KL divergence D_KL(P,Q) for sets of discrete probability distributions P, Q
function res = kl_divergence_discrete(P,Q)
if ~is_discrete(P) || ~is_discrete(Q)
    error('computation of the discrete KL divergence requires discrete distributions.'); end
% extract probability mass functions (and renormalize for good measure)
P = P{2}./repmat(sum(P{2},2),1,size(P{2},2));
Q = Q{2}./repmat(sum(Q{2},2),1,size(Q{2},2));
if ~isequal(size(P),size(Q))
    error('P and Q must have the same number of classes/bins and samples'); end
if any(~isfinite(P(:))) || any(~isfinite(Q(:)))
    error('the inputs contain non-finite values!'); end
res = P.*log(P./Q);
res(isnan(res)) = 0;
res = sum(res,2);
end



% compute the KL divergence D_KL(P,Q) for sets of gaussian probability distributions P,Q
function res = kl_divergence_gaussian(P,Q)
if ~is_gaussian(P) || ~is_gaussian(Q)
    error('P and Q must contain gaussian distributions'); end
% re-represent univariate distributions as a multivariate ones
P = quickif(strcmp(P{1},'norm'), cellfun(@(x){x(1),x(2)},mat2cell(P{2},ones(size(P{2},1),1),2),'UniformOutput',false), P{2});
Q = quickif(strcmp(Q{1},'norm'), cellfun(@(x){x(1),x(2)},mat2cell(Q{2},ones(size(Q{2},1),1),2),'UniformOutput',false), Q{2});
if length(P) ~= length(Q)
    error('P and Q must have the same number of distributions.'); end
res = zeros(length(P),1);
for i=1:length(P)
    [mu0,sig0] = P{i}{:};
    [mu1,sig1] = Q{i}{:};
    res(i) = (logdet(sig1)-logdet(sig0) + trace(sig1\sig0) + (mu1-mu0)' * (sig1\(mu1-mu0)) - length(mu0)) / 2;
end
end
