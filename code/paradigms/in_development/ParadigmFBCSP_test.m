classdef ParadigmFBCSP_test < ParadigmDataflowSimplified
    % Paradigm for complex oscillatory processes using the filter-bank CSP algorithm.
    % Result = para_multiband_csp(Input-Data, Operation-Mode, Options...)
    %
    % This is a test method which uses FIR/IIR filters computed on the continuous signal
    % instead of on the segments. It doesn't work online as-is (needs a new continuous filter for
    % that).
    %
    % Name:
    %   Filter-Bank CSP (test version)
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
                arg({'filter_type','FilterType'},'fft',{'fft','fir','iir-butter','iir-cheb1','iir-cheb2','iir-ellip'}, 'Filter type to use. Only FFT is online ready at this point as this is a test version.'), ...
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
                % set up transition bands around the frequency window
                f = args.freqwnds(w,:);
                f_trans = [max(0.1,f(1)-1),f(1),f(2),f(2)+1];
                switch args.filter_type
                    case 'fft'
                        freq_args{w} = arg_report('vals',@flt_spectrum,{'freq',args.freqwnds(w,:)});
                    case 'fir'
                        freq_args{w} = arg_report('vals',@flt_fir,{'Frequencies',f_trans});
                    case 'iir-butter'
                        freq_args{w} = arg_report('vals',@flt_iir,{'Frequencies',f_trans,'Type','butterworth'});
                    case 'iir-cheb1'
                        freq_args{w} = arg_report('vals',@flt_iir,{'Frequencies',f_trans,'Type','chebychev1'});
                    case 'iir-cheb2'
                        freq_args{w} = arg_report('vals',@flt_iir,{'Frequencies',f_trans,'Type','chebychev2'});
                    case 'iir-ellip'
                        freq_args{w} = arg_report('vals',@flt_iir,{'Frequencies',f_trans,'Type','elliptic'});
                    otherwise
                        error('Unsupported filter types: %s',hlp_tostring(args.filter_type,100));
                end
                for k=1:2
                    data = self.filter_data(args.signal,time_args{w},freq_args{w},args.filter_type,k);
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
            model = struct('filters',{filters},'patterns',{patterns},'time_args',{time_args},'freq_args',{freq_args},'ftype',{args.filter_type},'chanlocs',{args.signal.chanlocs});
        end
        
        function features = feature_extract(self,signal,featuremodel)
            W = length(featuremodel.freq_args);
            F = size(featuremodel.filters,2);
            T = size(signal.data,3);
            features = zeros(T,F);
            for w = 1:W
                % filter data in time & frequency
                data = self.filter_data(signal,featuremodel.time_args{w},featuremodel.freq_args{w},featuremodel.ftype);
                inds = (w-1)*(F/W)+(1:(F/W));
                for t=1:T
                    features(t,inds) = log(var(data.data(:,:,t)' * featuremodel.filters(:,inds))); end
            end
        end
        
        function data = filter_data(self,signal,time_args,freq_args,filter_type,picktrials)
            % filter trial subrange in time and frequency
            if strcmp(filter_type,'fft')
                data = exp_eval_optimized(flt_spectrum('signal',flt_window('signal',set_picktrials(signal,'rank',k),time_arg),freq_args));
            else
                % otherwise we hack the desired freq filter into the filter chain prior to set_makepos
                signal = signal.tracking.expression;
                [epo, pos] = utl_find_filter(signal,'set_makepos');
                signal_idx = find(strcmp(epo.parts,'signal'))+1;
                if strncmp(filter_type,'fir',3)
                    epo.parts{signal_idx} = flt_fir('signal',epo.parts{signal_idx},freq_args);
                elseif strncmp(filter_type,'iir',3)
                    epo.parts{signal_idx} = flt_iir('signal',epo.parts{signal_idx},freq_args);
                else
                    error('Unsupported filter type: %s',hlp_tostring(filter_type,100));
                end
                signal = subsasgn(signal,pos,epo);
                % ... and apply the time window etc at the end
                if exist('picktrials','var')
                    data = exp_eval_optimized(flt_window('signal',set_picktrials(signal,'rank',picktrials),time_args));
                else
                    data = exp_eval_optimized(flt_window('signal',signal,time_args));
                end
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

