classdef ParadigmRCSP < ParadigmDataflowSimplified
    % Advanced paradigm for oscillatory processes via the Regularized Common Spatial Patterns (CSP) algorithm(s).
    %
    % RCSP is an extension of the standard Common Spatial Patterns (CSP) algorithm [1] as defined in [2].
    % Implements a variety of improved CSP variants.
    %
    %
    % References:
    %  [1] Ramoser, H., Gerking, M., and Pfurtscheller, G. "Optimal spatial filtering of single trial EEG during imagined hand movement."
    %      IEEE Trans. Rehab. Eng 8 (2000), 446, 441.
    %  [2] Lotte, F., Guan, G., "Regularizing Common Spatial Patterns to Improve BCI Designs: Unified Theory and New Algorithms."
    %      IEEE Trans Biomed Eng 58, 2 (2011), 355-362.
    %
    % Name:
    %   Regularized Common Spatial Patterns
    %
    %                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
    %                                2011-09-09
    
    methods
      
        function defaults = preprocessing_defaults(self)
            defaults = {'FIRFilter',{'Frequencies',[6 8 28 32],'Type','minimum-phase'}, 'EpochExtraction',[0.5 3.5], 'Resampling',100};
        end
                
        function model = feature_adapt(self,varargin)
            args = arg_define(varargin, ...
                arg_norep('signal'), ...
                arg({'patterns','PatternPairs'},3,uint32([1 1 64 1000]),'Number of CSP patterns (times two).','cat','Feature Extraction','type','expression','shape','row'), ...
                arg({'alpha','ObjectiveRegularizer'},search(10.^(-10:-1)),[0 Inf],'Objective-function regularizer. This parameter regularizes the objective function according to a particular target matrix K.'), ...
                arg({'beta','CovariancePrior'},search(0:0.1:0.9),[0 1],'Covariance prior. This is the parameter beta that blends between the raw covariance (if 0) and a user-supplied "generic" covariance matrix for the respective class (e.g., pooled over subjects).'), ...
                arg({'gamma','CovarianceShrinkage'},search(0:0.1:0.9),[0 1],'Covariance shrinkage. This is the shrinkage parameter gamma that blends between the empirical covariance (if 0) and the identity matrix (if 1).'), ...
                arg({'priorcovs','PriorCovariances'},[],[],'Prior Covariance matrices. This is one matrix per class; needs to have the same number of channels as the data to be analyzed.'), ...
                arg({'objtarget','ObjectiveTarget'},[],[],'Objective-function target matrix. This matrix controls the behavior of the penalty parameter alpha (matrix K in Lotte''s framework ([2])). If this is empty, the identity matrix will be assumed (which yields Tikhonov-regularized CSP, or TRCSP).'));

            signal = args.signal;
            if signal.nbchan == 1
                error('CSP does intrinsically not support single-channel data (it is a spatial filter).'); end
            if signal.nbchan < args.patterns
                error('CSP prefers to work on at least as many channels as you request output patterns. Please reduce the number of pattern pairs.'); end
            if isempty(args.objtarget)
                args.objtarget = eye(size(signal.data,1)); end
            
            for k=1:2
                trials{k} = exp_eval_optimized(set_picktrials(signal,'rank',k));
                if strcmp(args.gamma,'auto')
                    if args.beta ~= 0
                        error('The ''auto'' covariance shrinkage is not compatible with a beta parameter ~= 0.'); end
                    % get raw covariance matrices
                    covar{k} = hlp_diskcache('featuremodels',@cov_shrink,reshape(trials{k}.data,size(trials{k}.data,1),[])');
                    covar{k}(~isfinite(covar{k})) = 0;
                else
                    % get raw covariance matrices
                    covar{k} = cov(reshape(trials{k}.data,size(trials{k}.data,1),[])');
                    covar{k}(~isfinite(covar{k})) = 0;
                    % incorporate the pior covariance, if available...
                    if ~isempty(args.priorcovs)
                        covar{k} = (1-args.beta)*covar{k} + args.beta*args.priorcovs{k}; end
                    % implement shrinkage towards identity
                    covar{k} = (1-args.gamma)*covar{k} + args.gamma*eye(size(signal.data,1));
                end
            end
            
            for k=1:2
                % implement objective-function regularization
                M{k} = covar{k}/(covar{3-k} + args.alpha*args.objtarget);
                try
                    [V{k},D{k}] = eig(M{k});
                    % invert to get forward projection
                    P{k} = inv(V{k});
                    if ~(all(isfinite(P{k})) && all(isfinite(V{k})) && all(isfinite(D{k})))
                        error('Divergence.'); end
                catch e
                    % keep going, this particular result will be weeded out by the parameter search
                    V{k} = randn(size(M{k}));
                    D{k} = randn(size(M{k}));
                    P{k} = randn(size(M{k}));
                end
            end
            
            % wrap up
            model.filters = [P{2}(1:args.patterns,:); P{1}(args.patterns:-1:1,:)]';
            model.patterns = [V{2}(:,1:args.patterns) V{1}(:,args.patterns:-1:1)]';
            model.chanlocs = signal.chanlocs;
        end
        
        function features = feature_extract(self,signal,featuremodel)
            features = zeros(size(signal.data,3),size(featuremodel.filters,2));
            for t=1:size(signal.data,3)
                features(t,:) = log(var(signal.data(:,:,t)'*featuremodel.filters)); end
        end
        
        function visualize_model(self,varargin) %#ok<*INUSD>
            args = arg_define([0 3],varargin, ...
                arg_norep({'myparent','Parent'},[],[],'Parent figure.'), ...
                arg_norep({'featuremodel','FeatureModel'},[],[],'Feature model. This is the part of the model that describes the feature extraction.'), ...
                arg_norep({'predictivemodel','PredictiveModel'},[],[],'Predictive model. This is the part of the model that describes the predictive mapping.'), ...
                arg({'patterns','PlotPatterns'},true,[],'Plot patterns instead of filters. Whether to plot spatial patterns (forward projections) rather than spatial filters.'), ...
                arg({'paper','PaperFigure'},false,[],'Use paper-style font sizes. Whether to generate a plot with font sizes etc. adjusted for paper.'));
            arg_toworkspace(args);
            
            % number of pairs, and index of pattern per subplot
            np = size(featuremodel.patterns,1)/2; 
            idx = [1:np 2*np:-1:np+1];
            % for each CSP pattern...
            for p=1:np*2
                subplot(2,np,p,'Parent',myparent);
                if args.patterns
                    topoplot(featuremodel.patterns(idx(p),:),featuremodel.chanlocs);
                else
                    topoplot(featuremodel.filters(:,idx(p)),featuremodel.chanlocs);
                end
                t = title(['CSP Pattern ' num2str(idx(p))]);
                if args.paper
                    set(t,'FontUnits','normalized');
                    set(t,'FontSize',0.1);                    
                end
            end
        end
        
        function layout = dialog_layout_defaults(self)
            layout = {'SignalProcessing.Resampling.SamplingRate', 'SignalProcessing.FIRFilter.Frequencies', ...
                'SignalProcessing.FIRFilter.Type', 'SignalProcessing.EpochExtraction', '', ...
                'Prediction.FeatureExtraction.PatternPairs', 'Prediction.FeatureExtraction.ObjectiveRegularizer', ...
                'Prediction.FeatureExtraction.CovariancePrior', 'Prediction.FeatureExtraction.CovarianceShrinkage','', ...
                'Prediction.MachineLearning.Learner'};
        end
        
        function tf = needs_voting(self)
            tf = true;
        end
        
    end
end
