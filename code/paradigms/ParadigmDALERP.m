classdef ParadigmDALERP < ParadigmDataflowSimplified
    % Advanced paradigm for slow cortical potentials via the Dual-Augmented Lagrange method.
    %
    % The DAL-LF paradigm is, like para_windowmeans, a general method for operating on slow cortical
    % potentials. It is a special case of a more general framework described in [1] (using only its
    % "first-order" detector); the general variant, which can in addition capture oscillatory processes,
    % is further explained in para_dal. DAL is the name of the optimization method, and not an accepted
    % or recognized name for BCI paradigms using it (but is used here for the lack of a better name).
    %
    % The paradigm does not make a clear distinction between signal processing, feature extraction and
    % machine learning, unlike most others, but instead is a jointly optimized mapping from raw signal
    % (epoch) to probabilistic prediction, using an efficient regularized optimization method (further
    % detailed in [2]). The method usually out-performs the windowed means paradigm, and in addition
    % does not require any user parameters aside from the epoch limits and lowpass filtering band, and
    % is therefore one of the most useful BCI paradigms. The major drawback is the required computation
    % time (and in some cases, the required memory -- which can be ameliorated by reducing the sampling
    % rate of the data) due to the need for regularization. For this reason, it is a good strategy to
    % first run the paradigm without regularization to get a ball-park estimate of the attainable
    % accuracy, and only run the complete regularization when it makes sense.
    %
    % Just like the windowed means paradigm, DAL-LF is applicable to a wide range of event-related and
    % non-event-related scenarios, some of which are listed in para_windowmeans.
    %
    % Example: Consider the goal of predicting whether a person perceives a fixated on-screen item as
    % being unexpected (and/or erroneous, non-rewarding) or not. A calibration data set for this task
    % could be annotated with an event for every gaze fixation made by the user (obtained from an eye
    % tracker) while reading short on-screen text fragments which are either semantically correct or
    % incorrect. The two event types which identify the conditions sare 'corr' and 'err'. From the
    % literature [4,5], it can be assumed that these events should be accompanied by a characteristic
    % slow cortical potential in the EEG, which allows to infer the condition. The 'learner' parameter
    % will be specified as the default (relatively fine-grained) search over possible DAL regularization
    % parameter values.
    %
    %   calib = io_loadset('data sets/john/reading-errors.eeg')
    %   myapproach = {'DALERP', 'SignalProcessing',{'EpochExtraction',[0 0.8]}};
    %   [loss,model,stats] = bci_train('Data',calib, 'Approach',myapproach, 'TargetMarkers',{'corr','err'});
    %
    %
    % References:
    %  [1] Ryota Tomioka and Klaus-Robert Mueller, "A regularized discriminative framework for EEG analysis with application to brain-computer interface",
    %      Neuroimage, 49 (1) pp. 415-432, 2010.
    %  [2] Ryota Tomioka & Masashi Sugiyama, "Dual Augmented Lagrangian Method for Efficient Sparse Reconstruction",
    %      IEEE Signal Proccesing Letters, 16 (12) pp. 1067-1070, 2009.
    %  [3] Marcel van Gerven, Ali Bahramisharif, Tom Heskes and Ole Jensen, "Selecting features for BCI control based on a covert spatial attention paradigm."
    %      Neural Networks 22 (9), 1271-1277, 2009
    %  [4] Gehring, W.J., Coles, M.G.H., Meyer, D.E., Donchin, E.
    %      "The error-related negativity: an event-related brain potential accompanying errors."
    %      Psychophysiology 27, 34-41, 1990
    %  [5] Oliveira, F.T.P., McDonald, J.J., Goodman, D.
    %      "Performance monitoring in the anterior cingulate is not all error related: expectancy deviation and the representation of action-outcome associations"
    %      Journal of Cognitive Neuroscience. 19(12), 1994-2004, 2007
    %
    % Name:
    %   Low-Frequency DAL
    %
    %                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
    %                           2010-06-25
    
    methods
        
        function defaults = preprocessing_defaults(self)
            defaults = {'FIRFilter',{[0.1 0.5],'highpass'},'EpochExtraction',[-1.5 1.5],'Resampling',60,'SpectralSelection',[0.1 15]};
        end
        
        function defaults = machine_learning_defaults(self)
            defaults = {'dal', 'Lambdas',2.^(10:-1.5:-5), 'NumFolds',5,'FoldMargin',1};
        end
        
        function model = feature_adapt(self,varargin)
            args = arg_define(varargin, ...
                arg_norep('signal'), ...
                arg({'normalizers','NormalizationExponents'},[-0.25,-0.25],[],'Normalization exponents [lhs, rhs]. Two-element array of powers for the left-hand-side and right-hand-side normalization matrices that are applied to the data from the region.','guru',true,'cat','Feature Extraction'), ...
                arg({'shrinkage_cov','ShrinkageCov'},true,[],'Use shrinkage covariance. This is slower but works better in the case of few trials.'), ...
                arg({'apply_to','ApplyTo'},'channels',{'channels','sources','components','full CSD'},'Apply classifier to. Allows to select the type of time series to apply this model to.','cat','Feature Extraction'));
            
            switch args.apply_to
                case 'channels'
                    X = num2cell(args.signal.data,[1 2]);
                case 'sources'
                    X = num2cell(args.signal.srcpot,[1 2]);
                case 'full CSD'
                    X = num2cell(args.signal.srcpot_all,[1 2]);
                    % compute diagonal covariance matrices right away, since it's hopeless to try to get a full cov
                    model.P = {diag(var(cat(2,X{:})'))^args.normalizers(1),diag(var(cat(1,X{:})))^args.normalizers(2)};
                case 'components'
                    if isempty(args.signal.icaact) && ~isempty(args.signal.icaweights)
                        args.signal.icaact = reshape((args.signal.icaweights*args.signal.icasphere)*args.signal.data(args.signal.icachansind,:), [], args.signal.pnts, args.signal.trials); end
                     X = num2cell(args.signal.icaact,[1 2]);
            end
            model.chanlocs = args.signal.chanlocs;
            if ~isfield(model,'P')
                if args.shrinkage_cov
                    model.P = {hlp_diskcache('featuremodels',@cov_shrink,cat(2,X{:})')^args.normalizers(1),hlp_diskcache('featuremodels',@cov_shrink,cat(1,X{:}))^args.normalizers(2)}; 
                else
                    model.P = {cov(cat(2,X{:})')^args.normalizers(1),cov(cat(1,X{:}))^args.normalizers(2)}; 
                end
            end
            model.times = args.signal.xmin + (0:args.signal.pnts-1)/args.signal.srate;
            model.apply_to = args.apply_to;
            
            % store some extra info
            model.cov = cov(args.signal.data(:,:)');            
            global tracking;
            tracking.inspection.dal_model = model;
        end
        
        function features = feature_extract(self,signal,featuremodel)
            switch featuremodel.apply_to
                case 'channels'
                    features = signal.data;
                case 'sources'
                    features = signal.srcpot;
                case 'full CSD'
                    features = signal.srcpot_all;
                case 'components'
                    features = signal.icaact;
                otherwise
                    error('Unsupported type of time series selected as source data: %s',hlp_tostring(featuremodel.apply_to));
            end
            for t=1:size(features,3)
                features(:,:,t) = featuremodel.P{1}*features(:,:,t)*featuremodel.P{2}; end
        end
        
        function visualize_model(self,varargin) %#ok<*INUSD>
            args = arg_define([0 3],varargin, ...
                arg_norep({'myparent','Parent'},[],[],'Parent figure.'), ...
                arg_norep({'fmodel','FeatureModel'},[],[],'Feature model. This is the part of the model that describes the feature extraction.'), ...
                arg_norep({'pmodel','PredictiveModel'},[],[],'Predictive model. This is the part of the model that describes the predictive mapping.'), ...
                arg({'maxcomps','MaxComponents'},Inf,[],'Maximum components to plot. Maximum number of components to plot (if too many).'), ...
                arg({'regcurve','PlotRegcurve'},true,[],'Plot regularization curve. Whether to plot the regularization curve.'), ...
                arg({'paper','PaperFigure'},false,[],'Use paper-style font sizes. Whether to generate a plot with font sizes etc. adjusted for paper.'), ...
                arg({'patterns','PlotPatterns'},true,[],'Plot patterns instead of filters. Whether to plot spatial patterns (forward projections) rather than spatial filters.'));                
            arg_toworkspace(args);
            
            % no parent? --> create new figure
            if isempty(myparent)
                myparent = figure('Name','Per-window weights'); end
            % get the spatial preprocessing matrix.
            P = fmodel.P{1};
            Q = fmodel.P{2};
            % obtain & reshape the model
            M = reshape(pmodel.model.w,size(P,2),[]);
            % do an SVD to get spatial and temporal filters
            [U,S,V] = svd(M);
            % display the model contents
            N = min(rank(M),args.maxcomps) + double(args.regcurve);
            px = ceil(sqrt(N));
            py = ceil(N/px);
            lim = -Inf;
            for x=1:N
                lim = max([lim;abs(inv(Q)*V(:,x)*S(x,x))]); end
            for x=1:N
                col = mod(x-1,px);
                row = floor((x-1) / px);
                idx = 1 + col + 2*row*px;
                if x < N || (x==N && ~args.regcurve)
                    subplot(2*py,px,idx,'Parent',myparent);
                    if args.patterns
                        topoplot(fmodel.cov*P*U(x,:)',fmodel.chanlocs);
                    else
                        topoplot(P*U(x,:)',fmodel.chanlocs);
                    end
                    t = title(sprintf('Component %.0f',x));
                    camzoom(1.2);
                    subplot(2*py,px,idx+px,'Parent',myparent);
                    p1 = plot(fmodel.times,inv(Q)*V(:,x)*S(x,x),'black');
                    ylim([-lim lim]);
                    hold; p2 = plot(fmodel.times,zeros(length(Q),1),'black--');
                    l1 = xlabel('Time in ms');
                    l2 = ylabel('Weight');
                elseif args.regcurve
                    subplot(2*py,px,idx+px,'Parent',myparent);
                    t = title('Regularization curve');
                    p1 = plot(mean(pmodel.model.losses)); p2=[];
                    l1 = xlabel('Regularization parameter #');
                    l2 = ylabel('Prediction loss');
                end
                if args.paper
                    set([p1,p2],'LineWidth',3);
                    set([l1,l2,t],'FontUnits','normalized');
                    set([l1,l2,t],'FontSize',0.1);
                    set(gca,'FontUnits','normalized');
                    set(gca,'FontSize',0.1);
                end
            end
        end
        
        function layout = dialog_layout_defaults(self)
            layout = {'SignalProcessing.Resampling.SamplingRate', 'SignalProcessing.FIRFilter.Frequencies',...
                'SignalProcessing.EpochExtraction', ...
                'SignalProcessing.SpectralSelection.FrequencySpecification', '', ...
                'Prediction.MachineLearning.Learner.Lambdas','Prediction.MachineLearning.Learner.LossFunction',...
                'Prediction.MachineLearning.Learner.Regularizer'};
        end
        
    end
end

