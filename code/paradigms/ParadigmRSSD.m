classdef ParadigmRSSD < ParadigmDataflowSimplified
    % Advanced paradigm for time-varying oscillatory processes.
    %
    % This paradigm allows to learn the joint space/time/frequency structure in EEG, under the assumption that
    % an overcomplete (AMICA-type) ICA decomposition can produce a reasonably complete coverage of potentially
    % relevant sources. For each of the resulting (more or less independent) source signals, a full time/freq
    % decomposition is computed in the time epoch of interest, and a regression model is learned which maps
    % from the (very high-dimensional) space of these features onto the desired output variable.
    % The regression is by default logistic, and is "rank-sparse" (using the trace norm regularization), meaning
    % that the time/freq weights learned for each signal component will tend to have low rank (reducing the effective
    % # of parameters) and also a sparse set of component signals will have non-zero weights.
    %
    % In addition, a prior can be imposed on the weights, which is parameterized in space (brain area), time, and 
    % frequency. This makes the overall model a fairly general-purpose approach to time-varying / non-stationary
    % oscillations.
    %
    % Notes:
    %  It takes quite a while to calibrate this type of model. We have recently changed a few things how dipole fits
    %  and amica solutions are represented to increase compatibility with EEGLAB - which has kicked up some dust in the
    %  code; if you get a crash here let us know so we can fix it.
    %
    % References:
    %  [1] Ryota Tomioka and Klaus-Robert Mueller, "A regularized discriminative framework for EEG analysis with application to brain-computer interface",
    %      Neuroimage, 49 (1) pp. 415-432, 2010.
    %  [2] Ryota Tomioka & Masashi Sugiyama, "Dual Augmented Lagrangian Method for Efficient Sparse Reconstruction",
    %      IEEE Signal Procesing Letters, 16 (12) pp. 1067-1070, 2009.
    %  [3] Makeig S, Debener S, Onton J, Delorme A, "Mining event-related brain dynamics"
    %      Trends in Cognitive Science, 8(5):204-210, 2004.
    %
    % Name:
    %   Regularized Spatio-Spectral Dynamics
    %
    %                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
    %                                2011-01-28

    
    methods
   
        function defaults = preprocessing_defaults(self)
            defaults = {'Resampling',128, 'FIRFilter',{[0.5 1],'highpass'}, 'EpochExtraction',[-2 2], 'ICA','beamica', 'DipoleFitting',{'ConfusionRange',7}};
        end
        
        function defaults = machine_learning_defaults(self)
            defaults = {'dal', 2.^(8:-0.25:-1), 'scaling','none'};
        end
        
        function [featuremodel,conditioningmodel,predictivemodel] = calibrate_prediction_function(self,varargin)
            args = arg_define(varargin, ...
                arg_norep('signal'), ...
                arg_sub({'fex','FeatureExtraction'},{},...
                    {arg({'spectral_prior','SpectralPrior'},@(f)1,[],'Spectral prior. Likelihood function of frequency in Hz.','guru',true), ...
                     arg({'temporal_prior','TemporalPrior'},@(t)1,[],'Temporal prior. Likelihood function of time in s.','guru',true), ...
                     arg({'spatial_prior','SpatialPrior'},@(p)1,[],'Spatial prior. Likelihood function of MNI coordinate vector.','guru',true), ...
                     arg({'anatomical_prior','AnatomicalPrior'},{'Left Cerebrum','Right Cerebrum','Left Cerebellum','Right Cerebellum','Left Brainstem','Right Brainstem','Inter-Hemispheric'},{'-- Hemispheres --','Left Cerebrum','Right Cerebrum','Left Cerebellum','Right Cerebellum','Left Brainstem','Right Brainstem','Inter-Hemispheric','-- Lobes --','Anterior Lobe','Frontal Lobe','Frontal-Temporal Space','Limbic Lobe','Medulla','Midbrain','Occipital Lobe','Parietal Lobe','Pons','Posterior Lobe','Sub-lobar','Temporal Lobe','-- Gyri --','Angular Gyrus','Anterior Cingulate','Caudate','Cerebellar Lingual','Cerebellar Tonsil','Cingulate Gyrus','Claustrum','Culmen','Culmen of Vermis','Cuneus','Declive','Declive of Vermis','Extra-Nuclear','Fastigium','Fourth Ventricle','Fusiform Gyrus','Inferior Frontal Gyrus','Inferior Occipital Gyrus','Inferior Parietal Lobule','Inferior Semi-Lunar Lobule','Inferior Temporal Gyrus','Insula','Lateral Ventricle','Lentiform Nucleus','Lingual Gyrus','Medial Frontal Gyrus','Middle Frontal Gyrus','Middle Occipital Gyrus','Middle Temporal Gyrus','Nodule','Orbital Gyrus','Paracentral Lobule','Parahippocampal Gyrus','Postcentral Gyrus','Posterior Cingulate','Precentral Gyrus','Precuneus','Pyramis','Pyramis of Vermis','Rectal Gyrus','Subcallosal Gyrus','Sub-Gyral','Superior Frontal Gyrus','Superior Occipital Gyrus','Superior Parietal Lobule','Superior Temporal Gyrus','Supramarginal Gyrus','Thalamus','Third Ventricle','Transverse Temporal Gyrus','Tuber','Tuber of Vermis','Uncus','Uvula','Uvula of Vermis'}, 'Anatomical prior. Select anatomical structures that are likely to contain processes of interest.'),...
                     ...
                     arg({'freq_oversample','FrequencyOversampling'},2,[],'Frequency oversampling.','guru',true), ...
                     arg({'sample_rate','SamplingRate'},15,[],'Sampling rate of the spectrum. In Hz.'), ...
                     arg({'freq_range','FrequencyRange'},[2.5 50],[],'Frequency range. In Hz, automatically clipped.'), ...
                     arg({'freq_scale','FrequencyScale'},'log',{'linear','log'},'Frequency scale.','guru',true), ...
                     arg({'win_size','WindowLength'},0.5,[],'Overlapped window length. For time-frequency estimation, in seconds.'), ...
                     arg({'wavelet_cycles','WaveletCycles'},[3 0.5],[],'Wavelet Cycles. At lower and upper frequency band edge.','guru',true) ...
                     arg({'specmap','SpectralMap'},'sqrt',{'sqrt','log','linear'},'Transform of power.','guru',true) ...
                    }, 'Parameters for the feature-adaptation function. These parameters control how features are statistically adapted and extracted from the filtered data before they are passed int othe machine learning stage.'), ...
                arg_sub({'cond','Conditioning'},{},@self.feature_adapt_conditioning,'Feature conditioning parameters. Allows to further process features for better usability with classifiers.'), ...
                arg_sub({'ml','MachineLearning'},{'Learner',self.machine_learning_defaults()},@ml_train,'Machine learning stage of the paradigm. Operates on the feature vectors that are produced by the feature-extraction stage.'));
            
            try
                data = self.rssd_load_overcomplete(args.signal);
                % read out component dipole fits
                if ~isempty(data.dipfit)
                    if isfield(data.dipfit,'multimodel')
                        dipfits = [data.dipfit.multimodel{:}];
                    else
                        dipfits = data.dipfit.model;
                    end
                else
                    dipfits = [];
                end
                if ~isempty(args.fex.anatomical_prior) && ~isequal(args.fex.anatomical_prior,false) && ~isempty(dipfits)
                    % if an anatomical prior was given, we can pre-prune the potenial ERSPs
                    ersprange = [];                    
                    for k=1:size(data.icaweights,1)
                        matches{k} = intersect(dipfits(k).structures,args.fex.anatomical_prior);
                        if ~isempty(matches{k})
                            ersprange(end+1) = k; end
                    end
                else
                    % otherwise we need to compute them all
                    ersprange = 1:size(data.icaweights,1);
                end
                
                % compute ERSPs
                disp('Now computing time/frequency decompositions...');
                [T,X,freqs,times] = evalc('self.ersp_compute(data,args,ersprange)'); %#ok<ASGLU>
                retain_ics = [];
                structures = {};
                probabilities = [];
                summed_probabilities = [];
                for k = ersprange
                    % get left/right normalization vectors lhs/rhs
                    lhs = args.fex.spectral_prior(freqs) .* sum(reshape(var(X{k},[],2),size(X{k},1),[]),2).^-0.25;
                    rhs = args.fex.temporal_prior(times) .* sum(reshape(var(X{k},[],1),[],size(X{k},3)),2).^-0.25;
                    % make sure that they are well-behaved
                    lhs(~isfinite(lhs)) = 0; medlhs = median(lhs); lhs(lhs>3*medlhs) = medlhs;
                    rhs(~isfinite(rhs)) = 0; medrhs = median(rhs); rhs(rhs>3*medrhs) = medrhs;
                    % turn into a scaling matrix
                    prior{k} =  diag(lhs) * ones(length(freqs),length(times)) * diag(rhs);
                    % incorporate the spatial prior
                    if ~isempty(dipfits)
                        prior{k} = prior{k} * args.fex.spatial_prior(dipfits(k).posxyz); end
                end
                for k = ersprange
                    % incorporate the anatomical prior
                    if ~isempty(args.fex.anatomical_prior) && ~isequal(args.fex.anatomical_prior,false) && ~isempty(dipfits)
                        [matches,idx] = intersect(dipfits(k).structures,args.fex.anatomical_prior); %#ok<ASGLU>
                        % sum the probabilities for being in each of the accepted structures (can be > 1 as the structures are highly overlapping)
                        prior{k} = sum(dipfits(k).probabilities(idx)) * prior{k};
                        structures{k} = dipfits(k).structures(idx);
                        probabilities{k} = dipfits(k).probabilities(idx);
                        summed_probabilities(k) = sum(dipfits(k).probabilities(idx));
                        % figure; topoplot(data.icawinv(:,k),data.chanlocs(data.icachansind),'electrodes','labels'); title([hlp_tostring(xdipfits(k).structures(idx)) ' - ' hlp_tostring(dipfits(k).probabilities(idx))]);
                    end
                    if ~all(all(prior{k}==0))
                        retain_ics(end+1) = k; end
                end
                featuremodel.prior = prior;
                featuremodel.retain_ics = retain_ics;
                featuremodel.args = args;
                blocksizes = cellfun(@size,featuremodel.prior,'UniformOutput',false);
                % keep track of some inspection information
                global tracking; %#ok<TLEV>
                tracking.inspection.rssd_prior = prior;
                tracking.inspection.rssd_mask = retain_ics;
                tracking.inspection.rssd_freqs = freqs;
                tracking.inspection.rssd_times = times;
                tracking.inspection.rssd_structures = structures;
                tracking.inspection.rssd_probabilities = probabilities;
                tracking.inspection.rssd_summed_probabilities = summed_probabilities;
                % update machine learning parameters
                args.ml.learner.shape = vertcat(blocksizes{retain_ics}); %vertcat(blocksizes{cellfun(@prod,blocksizes)});
                args.ml.learner.scaling = 'none';
                % extract features & target labels
                features = self.feature_extract(args.signal, featuremodel);
                targets = set_gettarget(args.signal);

                % adapt and apply feature conditioning
                conditioningmodel = self.feature_adapt_conditioning('features',features,'targets',targets,args.cond);
                [features,targets] = self.feature_apply_conditioning(features,targets,conditioningmodel);
                
                % run the machine learning stage
                predictivemodel = ml_train('data',{features,targets}, args.ml);
            catch
                disp('Glitch in ParadigmRSSD. halting...');
                keyboard
            end
        end        
        
        function features = feature_extract(self,signal,featuremodel)
            try
                data = self.rssd_load_overcomplete(signal);
                % compute ERSP features
                [T,features] = evalc('self.ersp_compute(data,featuremodel.args,featuremodel.retain_ics)'); %#ok<ASGLU>
                % scale each component-ERSP by the respective prior & sparsify
                for c=featuremodel.retain_ics
                    features{c} = bsxfun(@times,features{c},featuremodel.prior{c}); end
                % generate block-compressed version of the data
                features = reshape([features{featuremodel.retain_ics}],[],data.trials)';
            catch
                disp('Glitch in para_rssd. halting...');
                keyboard;
            end
        end
        
        % compute ERSP features
        function [X,freqs,times] = ersp_compute(self,data,args,icrange)
            args = args.fex;
            % compute IC activations...
            sig = zeros(size(data.icaweights,1),data.pnts,data.trials);
            for t=1:data.trials
                sig(:,:,t) = (data.icaweights*data.icasphere)*data.data(data.icachansind,:,t); end
            for k = icrange
                % compute a time-frequency decomposition
                [X{k},freqs,times] = hlp_microcache('tfc',@timefreq,squeeze(sig(k,:,:)),data.srate, ...
                    'winsize',data.srate*args.win_size, 'tlimits',[data.xmin data.xmax], ...
                    'detrend','off', 'itctype','phasecoher', 'wavelet',args.wavelet_cycles, ...
                    'padratio',args.freq_oversample, 'freqs',max(0,min(data.srate/2,args.freq_range)), ...
                    'freqscale',args.freq_scale, 'nfreqs',[], 'timestretch',{[],[]}, ...
                    'causal','off', 'ffttaper','hanning', 'subitc','off', 'wletmethod','dftfilt3', ...
                    'ntimesout',args.sample_rate*(data.xmax-data.xmin));
                % turn into power
                X{k} = X{k} .*conj(X{k});
                switch args.specmap
                    case 'sqrt'
                        X{k} = sqrt(X{k});
                    case 'log'
                        X{k} = log(X{k});
                    case 'linear'
                    otherwise
                        error('Unknown spectral transform.');
                end
                % inspect: k,figure;for j=1:size(X{k},3) imagesc(X{k}(:,:,j));waitforbuttonpress; end
            end
        end
        
        % load overcomplete Amica decomposition if applicable
        function data = rssd_load_overcomplete(self,data)
            if isfield(data.etc,'amica')
                % there is an amica decomposition: use it
                for m=1:size(data.etc.amica.W,3)
                    tmpW{m} = data.etc.amica.W(:,:,m);
                    tmpA{m} = data.etc.amica.A(:,:,m);
                end
                data.icaweights = cat(1,tmpW{:});
                data.icawinv = cat(2,tmpA{:});
                if isfield(data.etc.amica,'dipfit')
                    data.dipfit = data.etc.amica.dipfit; end
            end
        end
        
        function layout = dialog_layout_defaults(self)
            layout = {'SignalProcessing.Resampling.SamplingRate','SignalProcessing.IIRFilter.Frequencies', '', ...                
                'SignalProcessing.ICA.CleaningLevel.DataSetting', '', 'SignalProcessing.ICA.Variant.AmicaVersion', ...
                'SignalProcessing.ICA.Variant.NumModels','SignalProcessing.ICA.Variant.UseGridEngine', ...
                'SignalProcessing.ICA.Variant.NumProcessors', '', 'SignalProcessing.DipoleFitting.ConfusionRange', '', ...
                'SignalProcessing.EpochExtraction','', ...
                'Prediction.FeatureExtraction.SpectralPrior', 'Prediction.FeatureExtraction.TemporalPrior', ...
                'Prediction.FeatureExtraction.SpatialPrior','Prediction.FeatureExtraction.AnatomicalPrior', '', ...
                'Prediction.FeatureExtraction.FrequencyRange','Prediction.FeatureExtraction.WindowLength', ...
                'Prediction.FeatureExtraction.WaveletCycles','Prediction.FeatureExtraction.SamplingRate', '', ...
                'Prediction.MachineLearning.Learner.Lambdas', 'Prediction.MachineLearning.Learner.LossFunction','Prediction.MachineLearning.Learner.Regularizer'};
        end
        
        
        % --- helper functions for SCORD ---
        
        function result = export_model(self,model,filename)
            % get the raw weights
            w = model.predictivemodel.model.w;
            % and the block structure of the corresponding weight matrix
            shape = model.predictivemodel.model.shape;
            % extract the blocks
            ix = 0;
            W = {};
            for s=1:size(shape,1)
                ival = shape(s,1)*shape(s,2);
                W{s} = reshape(w(ix+1:ix+ival),shape(s,:));
                ix = ix+ival;
            end
            % get the IC references
            icidx = model.featuremodel.retain_ics;
            
            % get the ICA filter state
            icastate = utl_get_argument(utl_find_filter(model.tracking.filter_graph,'flt_ica'));
            % get dipoles, chanlocs and scalp maps
            dipoles = icastate.amica.model;
            chanlocs = icastate.root_chanlocs;
            scalpmaps = reshape(icastate.amica.A,size(icastate.amica.A,1),[]);

            % extract the subset that is actually used
            dipoles = dipoles(icidx);
            maps = scalpmaps(:,icidx);
            
            % construct result
            dipole_weights = W;
            dipole_models = dipoles;
            dipole_maps = maps;
            chanlocs = chanlocs;
            if exist('filename','var')
                save(filename,'dipole_weights','dipole_models','dipole_maps','chanlocs'); end
            result = struct('dipole_weights',{dipole_weights},'dipole_models',{dipole_models},'dipole_maps',{dipole_maps},'chanlocs',{chanlocs});
        end
                
    end
end
