classdef ParadigmFBCSP < ParadigmDataflowSimplified
    % Paradigm for complex oscillatory processes using the filter-bank CSP algorithm.
    % Result = para_multiband_csp(Input-Data, Operation-Mode, Options...)
    %
    % Filter-bank CSP [1,2] is a simple extension of the basic CSP method (see ParadigmCSP), in which for
    % each of several time/frequency filters a set of CSP filters is learned, followed by log-variance
    % feature extraction, concatenation of all features (over all chosen spectral filters) and
    % subsequent machine learning. It is not a general replacement for CSP due to the problem of
    % overfitting, but is very useful whenever oscillatory processes in different frequency bands (and
    % with different spatial topographies) are jointly active, and their concerted behavior must be
    % taken into account for a given prediction task. Filter-bank CSP can also be used to capture
    % oscillations in multiple time windows, instead of frequency windows (for example for the detection
    % of complex event-related dynamics).
    %
    % Since the dimensionality of the feature space is larger than in CSP, and since complex
    % interactions may be present, a more complex classifier than the default LDA may be necessary to
    % learn an appropriate model. On the other hand, more flexibility amplifies the risk of overfitting
    % (especially with only little calibration data), so that the performance should always be compared
    % to standard CSP (and Spec-CSP). Another reason is that complex (relevant) interactions between
    % different frequency bands are seemingly rarely observed in practice. The most important
    % user-configurable parameters are the selection regions in time and frequency and the learner
    % component.
    %
    % Typical applications would be those in which either complex event-related oscillatory dynamics
    % happen (for example when reacting to a particular stimulus) and/or where non-trivial interactions
    % between frequency bands (e.g. alpha/theta) are relevant, such as, for example, in workload
    % measurements.
    %
    % Example: Consider a calibration data set in which a subject is maintaining and updating
    % different number of items in his/her working memory at different times, e.g. while performing
    % the n-Back task [2]. Events with types 'n1','n2','n3' indicate challenge stimuli in which the
    % respective number of items is being processed by the person. The goal is to be able to predict
    % the working-memory load of the person following the presentation of such a memory-related
    % challenge. An epoch of 3 seconds relative to each challenge is selected, and three different
    % regions are chosen, two of them over the entire interval, covering the theta and alpha ryhthm,
    % respectively, and one region that is restricted to a window around the time of heaviest
    % cognitive processing. The three regions are specified as a cell array of flt_select
    % parameters.
    %
    %   data = io_loadset('data sets/mary/nback.eeg')
    %   myapproach = {'FBCSP' 'SignalProcessing',{'EpochExtraction',[-0.5 2.5]}, ...
    %       'Prediction', {'FeatureExtraction',{'FreqWindows',[4 6; 7 15; 7 15],'TimeWindows',[-0.5 2.5; -0.5 2.5; 0.25 1.25]}, ...
    %                      'MachineLearning',{'Learner','logreg'}}}
    %   [loss,model,stats] = bci_train('Data',data, 'Approach','ParadigmFBCSP, 'TargetMarkers',{'n1','n2','n3'})
    %
    % References;
    %   [1] Quadrianto Novi, Cuntai Guan, Tran Huy Dat, and Ping Xue, "Sub-band Common Spatial Pattern (SBCSP) for Brain-Computer Interface"
    %       Proceedings of the 3rd International IEEE EMBS Conference on Neural Engineering Kohala Coast, Hawaii, USA, May 2-5, 2007
    %   [2] Kai K. Ang, Zhang Y. Chin, Haihong Zhang, Cuntai Guan, "Filter Bank Common Spatial Pattern (FBCSP) in Brain-Computer Interface"
    %       In 2008 IEEE International Joint Conference on Neural Networks (IEEE World Congress on Computational Intelligence) (June 2008), pp. 2390-2397.
    %   [3] Owen, A. M., McMillan, K. M., Laird,A. R. & Bullmore, E. "N-back working memory paradigm: A meta-analysis of normative functional neuroimaging studies."
    %       Human Brain Mapping, 25, 46-59, 2005
    %
    % Name:
    %   Filter-Bank CSP
    %
    %                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
    %                           2010-04-29
    
    methods
      
        function defaults = preprocessing_defaults(self)
            % define the default pre-processing parameters of this paradigm
            defaults = {'EpochExtraction',[0.5 3.5],'Resampling',200};
        end
                
        function model = feature_adapt(self,varargin)
            % adapt a feature representation using the CSP algorithm
            args = arg_define(varargin, ...
                arg_norep('signal'), ...
                arg({'patterns','PatternPairs'},3,uint32([1 1 64 10000]),'CSP patterns per band (times two).','cat','Feature Extraction'), ...
                arg({'shrinkage_cov','ShrinkageCovariance'},false,[],'Shrinkage covariance estimator. Whether to use shrinkage to estimate the covariance matrices.'), ...
                arg({'robust_cov','RobustCovariance'},false,[],'Robust covariance estimation. Whether to use robust cov estimation.'), ...
                arg({'freqwnds','FreqWindows'},[0.5 3; 4 7; 8 12; 13 30; 31 42],[0 0.5 200 1000],'Frequency bands of interest. Matrix containing one row for the start and end of each frequency band from which CSP patterns shall be computed. Values in Hz.','cat','Feature Extraction'), ...
                arg({'timewnds','TimeWindows'},[],[],'Time windows of interest. Matrix containing one row for the start and end of each time window from which CSP patterns shall be computed. Values in seconds. If both this and the freqwnds parameter are non-empty, they should have the same number of rows.','cat','Feature Extraction'), ...
                arg({'winfunc','WindowFunction'},'rect',{'barthann','bartlett','blackman','blackmanharris','bohman','cheb','flattop','gauss','hamming','hann','kaiser','nuttall','parzen','rect','taylor','triang','tukey'},'Type of window function. Typical choices are rect (rectangular), hann, gauss, blackman and kaiser.'),...
                arg({'winparam','WindowParameter','param'},[],[],'Parameter of the window function. This is mandatory for cheb, kaiser and tukey and optional for some others.','shape','scalar'));
            
            if args.signal.nbchan == 1
                error('Multi-band CSP does intrinsically not support single-channel data (it is a spatial filter).'); end
            if args.signal.nbchan < args.patterns
                error('Multi-band CSP prefers to work on at least as many channels as you request output patterns. Please reduce the number of pattern pairs.'); end
            if ~isempty(args.freqwnds) && ~isempty(args.timewnds) && size(args.freqwnds,1) ~= size(args.timewnds,1)
                error('If both time and frequency windows are specified, both arrays must have the same number of rows (together they define the windows in time and frequency).'); end
            if isempty(args.timewnds)
                args.timewnds = zeros(size(args.freqwnds,1),0); end
            if isempty(args.freqwnds)
                args.freqwnds = zeros(size(args.timewnds,1),0); end
            
            filters = [];
            patterns = [];
            % for each window
            for w = 1:max(size(args.freqwnds,1),size(args.timewnds,1))
                % pre-parse arguments for flt_window and flt_spectrum (for fast subsequent online use)
                time_args{w} = arg_report('vals',@flt_window,{'time',{args.timewnds(w,:),args.winfunc,args.winparam}});
                freq_args{w} = arg_report('vals',@flt_spectrum,{'freq',args.freqwnds(w,:)});
                for k=1:2
                    % filter trial subrange in time and frequency
                    data = exp_eval_optimized(flt_spectrum('signal',flt_window('signal',set_picktrials(args.signal,'rank',k),time_args{w}),freq_args{w}));
                    if args.robust_cov
                        covar{k} = hlp_diskcache('featuremodels',@cov_blockgeom,reshape(data.data,size(data.data,1),[])',max([data.nbchan*2,data.srate*2,size(data,3)]));
                    else
                        if args.shrinkage_cov
                            covar{k} = hlp_diskcache('featuremodels',@cov_shrink,reshape(data.data,size(data.data,1),[])');
                        else
                            covar{k} = cov(reshape(data.data,size(data.data,1),[])');
                        end
                    end
                    covar{k}(~isfinite(covar{k})) = 0;
                end
                [V,D] = eig(covar{1},covar{1}+covar{2}); %#ok<NASGU>
                P = inv(V);                                
                filters = [filters V(:,[1:args.patterns end-args.patterns+1:end])];
                patterns = [patterns P([1:args.patterns end-args.patterns+1:end],:)'];
            end
            model = struct('filters',{filters},'patterns',{patterns},'time_args',{time_args},'freq_args',{freq_args},'chanlocs',{args.signal.chanlocs});
        end
        
        function features = feature_extract(self,signal,featuremodel)
            W = length(featuremodel.freq_args);
            F = size(featuremodel.filters,2);
            T = size(signal.data,3);
            features = zeros(T,F);
            for w = 1:W
                % filter data in time & frequency
                data = exp_eval_optimized(flt_spectrum('signal',flt_window('signal',signal,featuremodel.time_args{w}),featuremodel.freq_args{w}));
                inds = (w-1)*(F/W)+(1:(F/W));
                for t=1:T
                    features(t,inds) = log(var(data.data(:,:,t)' * featuremodel.filters(:,inds))); end
            end
        end
        
        function visualize_model(self,varargin) %#ok<*INUSD>
            args = arg_define([0 3],varargin, ...
                arg_norep({'myparent','Parent'},[],[],'Parent figure.'), ...
                arg_norep({'featuremodel','FeatureModel'},[],[],'Feature model. This is the part of the model that describes the feature extraction.'), ...
                arg_norep({'predictivemodel','PredictiveModel'},[],[],'Predictive model. This is the part of the model that describes the predictive mapping.'), ...
                arg({'patterns','PlotPatterns'},true,[],'Plot patterns instead of filters. Whether to plot spatial patterns (forward projections) rather than spatial filters.'), ...
                arg({'weight_scaled','WeightScaled'},false,[],'Scaled by weight. Whether to scale the patterns by weight.'));
            arg_toworkspace(args);
            
            % find the relevant components
            scores = predictivemodel.model.w;
            scores = sqrt(abs(scores));
            % optionally remove the bias if included in w
            if length(scores) == size(featuremodel.patterns,2)+1
                scores = scores(1:end-1); end 
            % frequency labels
            % titles = repmat({'delta','theta','alpha','beta','gamma'},8,1); titles = titles(:);
            % extract relevant patterns
            patterns = featuremodel.patterns(:,find(scores)); %#ok<FNDSB>
            filters = featuremodel.filters(:,find(scores)); %#ok<FNDSB>
            % plot them
            if args.weight_scaled
                if args.patterns
                    topoplot_grid(patterns,featuremodel.chanlocs,'scales',scores(find(scores))/max(scores)*1);
                else
                    topoplot_grid(filters,featuremodel.chanlocs,'scales',scores(find(scores))/max(scores)*1);
                end
            else
                if args.patterns
                    topoplot_grid(patterns,featuremodel.chanlocs);
                else
                    topoplot_grid(filters,featuremodel.chanlocs);
                end
            end
            % figure;
        end
        
        function layout = dialog_layout_defaults(self)
            % define the default configuration dialog layout
            layout = {'SignalProcessing.Resampling.SamplingRate', 'SignalProcessing.EpochExtraction', '', ...
                'Prediction.FeatureExtraction.FreqWindows', 'Prediction.FeatureExtraction.TimeWindows', ...
                'Prediction.FeatureExtraction.WindowFunction', '', 'Prediction.FeatureExtraction.PatternPairs', '', ...
                'Prediction.MachineLearning.Learner'};
        end
        
        function tf = needs_voting(self)
            tf = true;
        end
    end
end

