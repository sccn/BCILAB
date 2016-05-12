classdef ParadigmFBCSPMod < ParadigmDataflowSimplified
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
                arg({'freqwnds','FreqWindows'},[0.5 3; 4 7; 8 12; 13 30; 31 42],[0 0.5 200 1000],'Frequency bands of interest. Matrix containing one row for the start and end of each frequency band from which CSP patterns shall be computed. Values in Hz.','cat','Feature Extraction'), ...
                arg({'timewnds','TimeWindows'},[],[],'Time windows of interest. Matrix containing one row for the start and end of each time window from which CSP patterns shall be computed. Values in seconds. If both this and the freqwnds parameter are non-empty, they should have the same number of rows.','cat','Feature Extraction'), ...
                arg({'winfunc','WindowFunction'},'rect',{'barthann','bartlett','blackman','blackmanharris','bohman','cheb','flattop','gauss','hamming','hann','kaiser','nuttall','parzen','rect','taylor','triang','tukey'},'Type of window function. Typical choices are rect (rectangular), hann, gauss, blackman and kaiser.'),...
                arg({'winparam','WindowParameter','param'},[],[],'Parameter of the window function. This is mandatory for cheb, kaiser and tukey and optional for some others.','shape','scalar'),...
                arg({'nfft','NFFT'}, [], [],'Size of the FFT used in spectrum calculation. Default value is the greater of 256 or the next power of 2 greater than the length of the signal.' ),...
                arg({'winlen','WinLen'},100, [10, 1000], 'Divide the signal into sections of this length for Welch spectrum  calculation.'),...
                arg({'numoverlap','NumOverlap'}, [], [10, 1000], 'Number of overlap samples from section to next for Welch spectrum calculation.'));
            
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
            
            [signal, nof, freqwnds, timewnds, winfunc, winparam, nfft, winlen, numoverlap] = deal(args.signal, args.patterns, args.freqwnds,...
                args.timewnds, args.winfunc, args.winparam, args.nfft, args.winlen, args.numoverlap);
            
            [C,S,dum] = size(signal.data); 
            Fs = signal.srate;
            
            if isempty(numoverlap)
                numoverlap = floor(0.5*winlen); end
            
            if (winlen > S) || (numoverlap > winlen)
                error(' In Welch method, the length of the window should be smaller than the signal length, and the number of overlap should be smaller than the window length.');
            else
                win = window_func(winfunc,winlen,winparam);
                nwin = floor((S-numoverlap)/(winlen-numoverlap));
            end
            
            % The innfft is used internally to design and apply CSP filters
            if isempty(nfft)
                innfft = 2^(nextpow2(signal.pnts));
            else
                innfft = max(nfft, 2^(nextpow2(signal.pnts)))
                if innfft > nfft
                    error(' The chosen length of FFT is too short. '); end
            end
            
            allfreqs = 0:Fs/innfft:Fs;
            allfreqs = allfreqs(1:innfft);
            
            for c=1:2
                % compute the per-class epoched data X and its Fourier transform (along time), Xfft
                X{c} = exp_eval_optimized(set_picktrials(signal,'rank',c));
                [C,S,T] = size(X{c}.data);
                Xfft{c} = zeros(C,T,nfft,nwin);
                
                signal_idx = bsxfun(@plus,[1:winlen]',[0:nwin-1]*(winlen-numoverlap));
                Xdata_seg = repmat(X{c}.data,1,1,1,nwin);  %Xdata_seg -> C,S,T,nwin
                Xdata_seg2 = permute(Xdata_seg,[1 3 2  4]); %Xdata_seg2 ->C,T,S,nwin
                Xdata_seg3 = reshape(Xdata_seg2(:,:,signal_idx),C,T,[],nwin); % Xdata_seg3 -> C,T,winlen,nwin
                Xdata_win = bsxfun(@times,Xdata_seg3,reshape(win,1,1,winlen,1)); % Xdata_win -> C,T,winlen,nwin
                Xfft{c} = fft(Xdata_win,innfft,3); % Xfft -> C,T,nfft,nwin
            end
            
            filters = [];
            patterns = [];
            alphas = [];
            for fb=1:size(args.freqwnds,1)
                
                [freqs,findx] = getfgrid(Fs,innfft,freqwnds(fb,:));
                I = zeros(innfft,1); I(findx)=1;
                for cc=1:2
                    [C,S,T] = size(X{cc}.data);
                    Xspec{cc} = zeros(C,C);
                    %Xfft_crop -> C,T,numbands,numwin
                    Xfft_crop = Xfft{cc}(:, :,findx, :);
                    F = 2  * real(squeeze(sum(bsxfun(@times,conj(permute(Xfft_crop,[1,5,3,2,4])),permute(Xfft_crop,[5,1,3,2,4])),5))./nwin); % F{c}-> C,C,freqs,T
                    % compute the cross-spectrum  as an average over trials
                    Xspec{cc} = sum(squeeze(mean(F,4)),3);
                end
                
                [V,D] = eig(Xspec{1},Xspec{1}+Xspec{2});
                P = inv(V);
                filters = [filters V(:,[1:nof end-nof+1:end])];
                patterns = [patterns P([1:nof end-nof+1:end],:)'];
                alphas = [alphas repmat(I,1, 2*nof)];
                
            end
            freq_args.nfft = innfft;
            freq_args.win = win;
            freq_args.winlen = winlen;
            freq_args.numoverlap = numoverlap;
            freq_args.nwin = nwin;
            model = struct('filters',{filters},'patterns',{patterns},'alphas',{alphas},'freq_args',{freq_args},'chanlocs',{args.signal.chanlocs});
        end
        
        function features = feature_extract(self,signal,featuremodel)
            freq_args = featuremodel.freq_args;
            [nfft, win, winlen, numoverlap, nwin] = deal(freq_args.nfft, freq_args.win, freq_args.winlen, freq_args.numoverlap, freq_args.nwin);
            X = signal.data;
            [C,S,T] = size(X);            
            numf = size(featuremodel.filters,2);
            features = zeros(T,numf);
            
            Xtemp = permute(X,[3,2,1]); % T,S,C
            Xtemp2 = zeros(T,S,size(featuremodel.filters,2));
            for t=1:T
                Xtemp2(t,:,:) = squeeze(Xtemp(t,:,:)) * featuremodel.filters;
            end
            %Xtemp2 -> T,S,numf
            
            % Using welch method
            signal_idx = bsxfun(@plus,[1:winlen]',[0:nwin-1]*(winlen-numoverlap));
            Xdata_seg = repmat(Xtemp2,1,1,1,nwin);  %Xdata_seg -> T,S,numf,nwin
            Xdata_seg2 = permute(Xdata_seg,[1 3 2  4]); %Xdata_seg2 ->T,numf,S,nwin
            Xdata_seg3 = reshape(Xdata_seg2(:,:,signal_idx),T,numf,[],nwin); % Xdata_seg3 -> T,numf,winlen,nwin
            Xdata_win = bsxfun(@times,Xdata_seg3,reshape(win,1,1,winlen,1)); % Xdata_win -> T,numf,winlen,nwin
            Xfft = fft(Xdata_win,nfft,3); % Xfft -> T,numf,nfft,nwin
            Xdata1 = squeeze(sum(Xfft,4)./nwin); % Xdata1 ->T,numf,nfft
            Xdata = permute(Xdata1,[1,3,2]); % Xdata ->T,nfft,numf
            
            %Xdata = fft(Xtemp2,nfft,2); %Using regular fft
            
            Xdata2 = bsxfun(@times,Xdata,permute(featuremodel.alphas,[3,1,2]));
            Xdata3 = ifft(Xdata2,nfft,2);
            Xdata4 = 2*real(Xdata3(:,1:S,:));
            features = log(squeeze(var(Xdata4,0,2)));
            
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

