function pred = ml_predictgmm(trials, model)
% Prediction function for Gaussian Mixture Models.
% Prediction = ml_predictgmm(Trials, Model)
%
% In:
%   Trials  : the data a matrix, as in ml_predict
%
%   Model   : predictive model as produced by ml_trainlogreg
%
% Out:
%   Prediction  : discrete probability distribution, formatted as
%                 {'disc' [NxC] [Cx1]}, with element #2 being the per-class probability and 
%                 element #3 the original target values per class
%                 thus, the expected target values are Prediction{2}*Prediction{3}
%
% Examples:
%   targets might look like this: [-1 -1 1 -1 1 -1 -1 1 -1 -1 1 -1 -1 1 ...]'
%
%   model = ml_traingmm(data,targets)
%   p = ml_predictgmm(data, model); expectation = p{2}*p{3};
%   now expectation might look like this: [-0.6 -0.9 +0.4 -0.7 +0.8 -0.1 +0.5 +1.0 -0.9 +1.0 -1.0 -1.0 +1.0 ...]'
%
% See also:
%   ml_traingmm
%
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2010-04-05

% scale data
trials = hlp_applyscaling(trials,model.sc_info);

% predict probabilities
switch model.variant
    case {'avdp','vdp','bj','bjrnd','cdp','csb','vb'}
        for c=1:length(model.classes)
            model.class{c}.test_data = trials';
            % obtain log probabilities under each component
            tmp = vdpgm(model.data{c}',model.class{c});
            pdfmat(:,c) = exp(tmp.predictive_posterior');
        end
        % normalize the log probabilities, incorporate the prior, and turn them into true conditional probabilities
        probs = gmmb_normalize(pdfmat);
    case {'em','fj','gem'} 
        pdfmat = gmmb_pdf(trials, model.class);
        probs = gmmb_normalize(gmmb_weightprior(pdfmat, model.class));
end
pred = {'disc', probs, model.classes};
