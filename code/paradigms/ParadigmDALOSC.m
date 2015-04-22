classdef ParadigmDALOSC < ParadigmDataflowSimplified
    % Advanced paradigm for oscillatory processes via the Dual-Augmented Lagrange method.
    %
    % The DAL-HF paradigm is an easy-to-use and generic method for operating on oscillatory processes.
    % It is a special case of the general framework described in [1] (using only the "second-order"
    % detector); the general variant, which can use multiple bands and time windows, as well as slow
    % cortical potentials under a joint optimality criterion, is further explained in para_dal. DAL is
    % the name of the optimization method [2], and not an accepted or recognized name for BCI paradigms
    % using it (but is used here for the lack of a better name).
    %
    % Like CSP (para_csp), DAL-HF can be used for a wide range of oscillatory processes and should
    % theoretically out-perform it in many cases, at the cost of much increased computation times (due
    % to the need for regularization). Like for CSP, the frequency band and time window on which DAL-LF
    % shall operate needs to be determined a priori. The frequency band is usually chosen to be the
    % alpha and/or beta band for the majority of tasks, but under some circumstances, it can be chosen
    % heuristically, as in [3].  The learner parameter is either chosen to be 'dal', to run without
    % regularization for quick results, or {'dal',search(2.^(-6:0.5:10)} to run with full
    % regularization, which takes much longer but is the proper way to apply the method.
    %
    % Example: Consider a user of an eye tracking software to control a mouse cursor. The user wants to
    % execute a control thought (e.g., sensorimotor) in order to issue selections (i.e., a 'click'). We
    % assume that the control signal is reflected in oscillatory processes. either in the alpha or the
    % beta band, and that it takes approx. 2 seconds to think it. For a possible calibration design,
    % see, e.g., [4], in which the data set contains events in which the user is either in the 'search'
    % condition or in the 'select' condition. A predictive model can then be obtained with DAL-HF as
    % follows:
    %
    %   calib = io_loadset('data sets/john/selection.eeg')
    %   myapproach = {'DALOSC', 'SignalProcessing', {'EpochExtraction',[0 2]}};
    %   [loss,model,stats] = bci_train('Data',calib, 'Approach',myapproach, 'TargetMarkers',{'search','select'});
    %
    %
    % References:
    %  [1] Ryota Tomioka and Klaus-Robert Mueller, "A regularized discriminative framework for EEG analysis with application to brain-computer interface",
    %      Neuroimage, 49 (1) pp. 415-432, 2010.
    %  [2] Ryota Tomioka & Masashi Sugiyama, "Dual Augmented Lagrangian Method for Efficient Sparse Reconstruction",
    %      IEEE Signal Proccesing Letters, 16 (12) pp. 1067-1070, 2009.
    %  [3] Benjamin Blankertz, Ryota Tomioka, Steven Lemm, Motoaki Kawanabe, and Klaus-Robert Mueller.
    %      "Optimizing spatial filters for robust EEG single-trial analysis."
    %      IEEE Signal Process Mag, 25(1):41-56, January 2008
    %  [4] Zander T.O., Gaertner M., Kothe C., Vilimek R. "Combining Eye Gaze Input with a Brain-Computer Interface for Touchless Human-Computer Interaction"
    %      International Journal of Human-Computer Interaction, in press.
    %
    % Name:
    %   High-Frequency DAL
    %
    %                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
    %                           2010-06-25
    
    methods
      
        function defaults = preprocessing_defaults(self)
            defaults = {'FIRFilter',{'Frequencies',[6 8 28 32],'Type','minimum-phase'}, 'EpochExtraction',[0.5 3.5], 'Resampling',100};
        end

        function defaults = machine_learning_defaults(self)
            defaults = {'dal', 'Lambdas',2.^(4:-0.25:-3), 'NumFolds',5,'FoldMargin',1};
        end
        
        function model = feature_adapt(self,varargin)
            args = arg_define(varargin, ...
                arg_norep('signal'), ...
                arg({'normalizers','NormalizationExponents'},[-0.5,-0.5],[],'Normalization exponents [lhs, rhs]. Two-element array of powers for the left-hand-side and right-hand-side normalization matrices that are applied to the data from the region.','cat','Feature Extraction','guru',true), ...
                arg({'shrinkage_rescalers','ShrinkageRescalers'},true,[],'Use shrinkage to estimate rescalers. If enabled, the rescaling matrices will be estimated using shrinkage.'), ...
                arg({'normalize','Normalize'},true,[],'Perform data normalization.'));
            
            if args.normalize
                X = num2cell(args.signal.data,[1 2]);
                for t=1:length(X)
                    X{t} = cov(X{t}'); end
                % calc rescaling matrices
                if args.shrinkage_rescalers
                    model.P = {hlp_diskcache('featuremodels',@cov_shrink,cat(2,X{:})')^args.normalizers(1),hlp_diskcache('featuremodels',@cov_shrink,cat(1,X{:}))^args.normalizers(2)};
                else
                    model.P = {cov(cat(2,X{:})')^args.normalizers(1),cov(cat(1,X{:}))^args.normalizers(2)};
                end
                % make sure that the data is properly normalized so that our lambda range does not have 
                % to be tuned to fit the data scale
                Y = zeros(size(args.signal.data,1),size(args.signal.data,1),size(args.signal.data,3));
                for t=1:size(args.signal.data,3)
                    Y(:,:,t) = model.P{1}*cov(args.signal.data(:,:,t)')*model.P{2}; end                    
                model.scale = 1/sum(sum(sqrt(var(Y,[],3)))/(size(Y,1)*size(Y,2)));
            end
            model.chanlocs = args.signal.chanlocs;
            model.cov = cov(args.signal.data(:,:)');
            model.normalize = args.normalize;
        end
        
        function features = feature_extract(self,signal,featuremodel)
            features = zeros(size(signal.data,1),size(signal.data,1),size(signal.data,3));
            if featuremodel.normalize
                for t=1:size(signal.data,3)
                    features(:,:,t) = featuremodel.scale*(featuremodel.P{1}*cov(signal.data(:,:,t)')*featuremodel.P{2}); end
            else
                for t=1:size(signal.data,3)
                    features(:,:,t) = cov(signal.data(:,:,t)'); end
            end
        end

        function visualize_model(self,varargin) %#ok<*INUSD>
            args = arg_define([0 3],varargin, ...
                arg_norep({'myparent','Parent'},[],[],'Parent figure.'), ...
                arg_norep({'fmodel','FeatureModel'},[],[],'Feature model. This is the part of the model that describes the feature extraction.'), ...
                arg_norep({'pmodel','PredictiveModel'},[],[],'Predictive model. This is the part of the model that describes the predictive mapping.'), ...
                arg({'maxcomps','MaxComponents'},Inf,[],'Maximum components to plot. Maximum number of components to plot (if too many).'), ...
                arg({'patterns','PlotPatterns'},true,[],'Plot patterns instead of filters. Whether to plot spatial patterns (forward projections) rather than spatial filters.'), ...
                arg({'paper','PaperFigure'},false,[],'Use paper-style font sizes. Whether to generate a plot with font sizes etc. adjusted for paper.'));
            arg_toworkspace(args);
            
            if isempty(myparent)
                myparent = figure('Name','Per-window weights'); end
            % get the spatial preprocessing matrix.
            P = fmodel.P{1};
            % obtain & reshape the model
            M = reshape(pmodel.model.w,size(P,2),[]);
            % do an SVD to get spatial and temporal filters
            [U,S,V] = svd(M);
            % display the model contents
            N = min(rank(M),args.maxcomps);
            weights = P*U(1:N,:)';
            if args.patterns
                weights = fmodel.cov*weights; end
            topoplot_grid(weights,fmodel.chanlocs);
        end
        
        function layout = dialog_layout_defaults(self)
            layout = {'SignalProcessing.Resampling.SamplingRate', 'SignalProcessing.FIRFilter.Frequencies', ...
                'SignalProcessing.FIRFilter.Type', 'SignalProcessing.EpochExtraction', '', ...
                'Prediction.MachineLearning.Learner.Lambdas','Prediction.MachineLearning.Learner.LossFunction',...
                'Prediction.MachineLearning.Learner.Regularizer'};
        end
        
    end
end

