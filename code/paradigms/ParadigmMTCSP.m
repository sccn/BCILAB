classdef ParadigmMTCSP < ParadigmDataflowSimplified
    % Experimental paradigm for all-frequency Common Spatial Patterns.
    %
    % The basic idea is to calculate CSP for each covariance matrix in the cross-spectrum, 
    % and to use multi-taper spectral estimation to ensure an optimal tradeoff between
    % spectral precision and estimation noise. The default classifier is sparse logistic
    % regression with elastic-net penalty.
    %
    % This paradigm also implements a second approach in which the cross-spectrum is not spatially
    % filtered but directly submitted to the classifier (Disciplined Cross-Spectral Regression).
    %
    % Name:
    %   Multi-Taper CSP
    %
    %                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
    %                           2013-04-26
    
    methods
      
        function defaults = preprocessing_defaults(self)
            % define the default pre-processing parameters of this paradigm
            defaults = {'FIRFilter',{[0.5 1],'highpass'}, 'EpochExtraction',[0.5 3.5],'Resampling',200};
        end
        
        function defaults = machine_learning_defaults(self)
            defaults = {'logreg','variant',{'lars','ElasticMixing',0.25}};
        end
                
        function model = feature_adapt(self,varargin)
            % adapt a feature representation using the CSP algorithm
            args = arg_define(varargin, ...
                arg_norep('signal'), ...
                arg({'timewnds','TimeWindows'},[],[],'Time windows of interest. Matrix containing one row for the start and end of each time window for which CSP patterns shall be computed. Values in seconds. If both this and the freqwnds parameter are non-empty, they should have the same number of rows.'), ...
                arg({'winfunc','WindowFunction'},'rect',{'bartlett','barthann','blackman','blackmanharris','flattop','gauss','hamming','hann','kaiser','lanczos','nuttall','rect','triang'},'Type of window function. Typical choices are rect (rectangular), hann, gauss, blackman and kaiser.'),...
                arg({'winparam','WindowParameter','param'},[],[],'Parameter of the window function. This is mandatory for cheb, kaiser and tukey and optional for some others.','shape','scalar'), ...
                arg_sub({'spectral_estimation','SpectralEstimation'},{},@utl_calc_crossspec, 'Spectral estimation parameters.'), ...
                arg({'patterns','PatternPairs'},3,uint32([1 1 64 10000]),'CSP patterns per frequency (times two).'), ...
                arg({'whycsp','SkipCSP','WhyCSP'},false,[],'Classify cross-spectrum directly. This results in much higher-dimensional features, but can be approached with appropriately regularized classifiers (see ml_trainproximal).'), ...
                arg({'normalize_spectrum','NormalizeSpectrum'},false,[],'Normalize the spectrum. Recommended if using sophisticated regularized classifiers.'), ...
                arg({'logtransform','LogTransform'},false,[],'Log-transform output. Log-transformed spectra are more likely to be separable by a linear classifier.'), ...
                arg({'vectorize_features','VectorizeFeatures'},true,[],'Vectorize the features. For compatibility with basic classifiers.'));
            
            if args.signal.nbchan == 1
                error('Multi-taper CSP does intrinsically not support single-channel data (it is a spatial filter).'); end
            if args.signal.nbchan < args.patterns
                error('Multi-taper CSP prefers to work on at least as many channels as you request output patterns. Please reduce the number of pattern pairs.'); end
            if isempty(args.timewnds)
                args.timewnds = struct(); end
            [C,dummy,T] = size(args.signal.data); %#ok<NASGU,ASGLU>
            if args.whycsp
                % shortcut
                for w = size(args.timewnds,1):-1:1
                    time_args{w} = arg_report('vals',@flt_window,{'time',{args.timewnds(w,:),args.winfunc,args.winparam}}); end
                model = struct('spec_args',{args.spectral_estimation}, 'time_args',{time_args},'chanlocs',{args.signal.chanlocs},'vectorize_features',{args.vectorize_features},'whycsp',{args.whycsp},'normalize_spectrum',{args.normalize_spectrum});
            else
                covar = {}; mean_covar = {}; weighted_covar = {};
                % for each time window...
                for w = size(args.timewnds,1):-1:1
                    if length(unique([args.signal.event.target]))>2
                        % SPoC version
                        time_args{w} = arg_report('vals',@flt_window,{'time',{args.timewnds(w,:),args.winfunc,args.winparam}});
                        % calc weighted and average cross-spectra
                        [mean_covar{w},weighted_covar{w}] = hlp_diskcache('featuremodels',@utl_calc_crossspec,args.spectral_estimation,'sum_weights',[args.signal.epoch.target],'signal',exp_eval_optimized(flt_window('signal',args.signal,time_args{w})));
                        mean_covar{w}(~isfinite(mean_covar{w})) = 0; weighted_covar{w}(~isfinite(weighted_covar{w})) = 0;
                        % calculate spatial filters for each frequency
                        for f=size(mean_covar{w},1):-1:1
                            [V,D] = eig(reshape(weighted_covar{w}(f,:,:),C,C),reshape(mean_covar{w,1}(f,:,:),C,C)); %#ok<NASGU>
                            P = inv(V);
                            % retain k best filters/patterns at both ends of the eigenvalue spectrum
                            filters(w,f,:,:) = real(V(:,[1:args.patterns end-args.patterns+1:end]));
                            patterns(w,f,:,:) = real(P([1:args.patterns end-args.patterns+1:end],:))';
                        end
                    else
                        % CSP version
                        for k=1:2
                            subset = exp_eval_optimized(set_picktrials(args.signal,'rank',k));
                            % pre-parse arguments for flt_window and flt_spectrum (for fast subsequent online use)
                            time_args{w} = arg_report('vals',@flt_window,{'time',{args.timewnds(w,:),args.winfunc,args.winparam}});
                            % calc cross-spectrum for the given windowed data subset
                            covar{w,k} = hlp_diskcache('featuremodels',@utl_calc_crossspec,args.spectral_estimation,'signal',exp_eval_optimized(flt_window('signal',subset,time_args{w})));
                            covar{w,k}(~isfinite(covar{w,k})) = 0;
                        end
                        % solve a CSP instance for each frequency
                        for f=size(covar{w,1},1):-1:1
                            [V,D] = eig(reshape(covar{w,1}(f,:,:),C,C),reshape(covar{w,1}(f,:,:),C,C)+reshape(covar{w,2}(f,:,:),C,C)); %#ok<NASGU>
                            P = inv(V);
                            filters(w,f,:,:) = real(V(:,[1:args.patterns end-args.patterns+1:end]));
                            patterns(w,f,:,:) = real(P([1:args.patterns end-args.patterns+1:end],:))';
                        end                
                    end
                end
                model = struct('filters',{filters},'patterns',{patterns},'time_args',{time_args},'spec_args',{args.spectral_estimation}, 'covar',{covar}, 'mean_covar',{mean_covar}, 'weighted_covar',{weighted_covar}, 'chanlocs',{args.signal.chanlocs},'vectorize_features',{args.vectorize_features},'whycsp',{args.whycsp},'normalize_spectrum',{args.normalize_spectrum},'logtransform',{args.logtransform});
            end
            global tracking; %#ok<TLEV>
            tracking.inspection.signal = args.signal;
            tracking.inspection.chanlocs = args.signal.chanlocs;
        end
                
        function features = feature_extract(self,signal,featuremodel)
            for w = length(featuremodel.time_args):-1:1
                % extract time window
                wnd = exp_eval_optimized(flt_window('signal',signal,featuremodel.time_args{w}));
                % extract cross-spectral features (note: gigantic!)
                if featuremodel.whycsp
                    % W x F x C x C x T
                    features(w,:,:,:,:) = utl_calc_crossspec(featuremodel.spec_args,'signal',wnd,'feature_filters',false);
                    if featuremodel.normalize_spectrum
                        nfreqs = size(features,2);
                        freqs = featuremodel.spec_args.freqwnd;
                        freqs = freqs(1):(freqs(2)-freqs(1))/(nfreqs-1):freqs(2);
                        features = bsxfun(@times,features,max(1,freqs));
                    end
                else
                    % F x W x P x T
                    if onl_isonline
                        features(:,w,:,:) = utl_calc_crossspec(featuremodel.spec_args,'signal',wnd,'feature_filters',squeeze(featuremodel.filters(w,:,:,:)));
                    else
                        features(:,w,:,:) = hlp_diskcache('features',@utl_calc_crossspec,featuremodel.spec_args,'signal',wnd,'feature_filters',squeeze(featuremodel.filters(w,:,:,:)));
                    end
                    if featuremodel.normalize_spectrum
                        nfreqs = size(features,1);
                        freqs = featuremodel.spec_args.freqwnd;
                        freqs = freqs(1):(freqs(2)-freqs(1))/(nfreqs-1):freqs(2);
                        features = bsxfun(@times,features,max(1,1./freqs'));
                    end                    
                    if featuremodel.logtransform
                        features = log(features); end
                end
            end
            % apply minimal conditioning to features
            features = real(features);            
            features(~isfinite(features)) = 0;
            % do final vectorization if desired
            if featuremodel.vectorize_features
                features = reshape(features,[],signal.trials)'; end
        end
        
        function visualize_model(self,parent,featuremodel,predictivemodel,varargin) %#ok<*INUSD>
            % no visualization yet
        end
        
        function layout = dialog_layout_defaults(self)
            % define the default configuration dialog layout
            layout = {'SignalProcessing.Resampling.SamplingRate', 'SignalProcessing.EpochExtraction', '', ...
                'Prediction.FeatureExtraction.TimeWindows','Prediction.FeatureExtraction.WindowFunction', '', ...
                'Prediction.FeatureExtraction.SpectralEstimation.FrequencyRange', ...
                'Prediction.FeatureExtraction.SpectralEstimation.TimeBandwidth', ...
                'Prediction.FeatureExtraction.SpectralEstimation.SubsampleSpectrum', ...
                'Prediction.FeatureExtraction.SpectralEstimation.RobustEstimation', '' ...
                'Prediction.FeatureExtraction.PatternPairs', 'Prediction.FeatureExtraction.VectorizeFeatures', '', ...
                'Prediction.MachineLearning.Learner'};
        end
        
        function tf = needs_voting(self)
            tf = false;
        end
    end
end

