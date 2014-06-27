classdef ParadigmRSSD2 < ParadigmDataflowSimplified
    % Advanced paradigm for time-varying oscillatory processes.
    %
    % This paradigm allows to learn the joint space/time/frequency structure in EEG, under the
    % assumption that an overcomplete ICA decomposition can produce a reasonably complete coverage
    % of potentially relevant sources. For each of the resulting (more or less independent) source
    % signals, a full time/freq decomposition is computed in the time epoch of interest, and a
    % regression model is learned which maps from the (very high-dimensional) space of these
    % features onto the desired output variable. The regression is by default logistic, and is
    % "rank-sparse" (using the trace norm regularization), meaning that the time/freq weights
    % learned for each signal component will tend to have low rank (reducing the effective
    % # of parameters) and also a sparse set of component signals will have non-zero weights.
    %
    % In addition, a prior can be imposed on the weights, which is parameterized in space (brain area), time, and 
    % frequency. This makes the overall model a fairly general-purpose approach to time-varying / non-stationary
    % oscillations.
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
    %   Regularized Spatio-Spectral Dynamics v2
    %
    %                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
    %                                2011-01-28

    
    methods
        function defaults = preprocessing_defaults(self)
            defaults = { ...
                'Resampling',128, ...
                'FIRFilter',{[0.5 1],'highpass'}, ...
                'EpochExtraction',[-2 2], ...
                'ICA','beamica', ...
                ...%'DipoleFitting',{'ConfusionRange',7}, ...
                'Projection','on', ...
                'ERSPTransform',{'WindowStep',1/15,'SpectralMap','sqrt'}, ...
                };
        end
        
        function defaults = machine_learning_defaults(self)
            defaults = {'dal', 2.^(4:-0.25:-3), 'scaling','none'};
            %defaults = {'logreg', 'variant',{'lars','ElasticMixing',0.5}};
        end
                
        function model = feature_adapt(self,varargin)
            args = arg_define(varargin, ...
                    arg_norep('signal'), ...
                    arg({'spectral_prior','SpectralPrior'},@(f)1,[],'Spectral prior. Likelihood function of frequency in Hz.','guru',true), ...
                    arg({'temporal_prior','TemporalPrior'},@(t)1,[],'Temporal prior. Likelihood function of time in s.','guru',true), ...
                    arg({'spatial_prior','SpatialPrior'},@(p)1,[],'Spatial prior. Likelihood function of MNI coordinate vector.','guru',true), ...
                    arg({'anatomical_prior','AnatomicalPrior'},{'Left Cerebrum','Right Cerebrum','Left Cerebellum','Right Cerebellum','Left Brainstem','Right Brainstem','Inter-Hemispheric'},{'-- Hemispheres --','Left Cerebrum','Right Cerebrum','Left Cerebellum','Right Cerebellum','Left Brainstem','Right Brainstem','Inter-Hemispheric','-- Lobes --','Anterior Lobe','Frontal Lobe','Frontal-Temporal Space','Limbic Lobe','Medulla','Midbrain','Occipital Lobe','Parietal Lobe','Pons','Posterior Lobe','Sub-lobar','Temporal Lobe','-- Gyri --','Angular Gyrus','Anterior Cingulate','Caudate','Cerebellar Lingual','Cerebellar Tonsil','Cingulate Gyrus','Claustrum','Culmen','Culmen of Vermis','Cuneus','Declive','Declive of Vermis','Extra-Nuclear','Fastigium','Fourth Ventricle','Fusiform Gyrus','Inferior Frontal Gyrus','Inferior Occipital Gyrus','Inferior Parietal Lobule','Inferior Semi-Lunar Lobule','Inferior Temporal Gyrus','Insula','Lateral Ventricle','Lentiform Nucleus','Lingual Gyrus','Medial Frontal Gyrus','Middle Frontal Gyrus','Middle Occipital Gyrus','Middle Temporal Gyrus','Nodule','Orbital Gyrus','Paracentral Lobule','Parahippocampal Gyrus','Postcentral Gyrus','Posterior Cingulate','Precentral Gyrus','Precuneus','Pyramis','Pyramis of Vermis','Rectal Gyrus','Subcallosal Gyrus','Sub-Gyral','Superior Frontal Gyrus','Superior Occipital Gyrus','Superior Parietal Lobule','Superior Temporal Gyrus','Supramarginal Gyrus','Thalamus','Third Ventricle','Transverse Temporal Gyrus','Tuber','Tuber of Vermis','Uncus','Uvula','Uvula of Vermis'}, 'Anatomical prior. Select anatomical structures that are likely to contain processes of interest.'),...
                    arg({'vectorize_features','VectorizeFeatures'},true,[],'Vectorize feature tensors. This is for classifiers that cannot handle matrix or tensor-shaped features.'),...
                    arg({'normalize_features','NormalizeFeatures'},true,[],'Normalize time/frequency features. If enabled, features will be normalized by a rank-1 normalization matrix (rather than pixelwise).'));
%
            model = rmfield(args,'signal');
            % determine data rescaling factors
            if args.normalize_features
                tdata = permute(args.signal.data,[1 2 4 3]);
                sdata = permute(args.signal.data,[1 4 2 3]);
                for c=size(args.signal.data,1):-1:1
                    % temporal scale vector
                    temp = 1./median(abs(bsxfun(@minus,tdata(c,:,:),median(tdata(c,:,:),3))),3);
                    % spectral scale vector
                    spec = 1./median(abs(bsxfun(@minus,sdata(c,:,:),median(sdata(c,:,:),3))),3);
                    % time/freq scaling matrix
                    model.scaling(c,:,:) = sqrt(temp')*sqrt(spec);
                end
            end
            % determine model shape
            model.shape = size(args.signal.data);
            model.shape = model.shape([2 4 1]);
        end
        
        function [features,shape] = feature_extract(self,signal,featuremodel)
            signal.data = permute(signal.data,[1 2 4 3]);
            % optionally apply normalization
            if featuremodel.normalize_features
                features = bsxfun(@times,signal.data,featuremodel.scaling);
            else
                features = signal.data;
            end
            features = permute(features,[2 3 1 4]);
            % determine feature shape
            siz = size(features);
            shape = siz(1:3);
            % do final vectorization if desired
            if featuremodel.vectorize_features
                features = reshape(features,[],signal.trials)'; end
        end
        
        function layout = dialog_layout_defaults(self)
            layout = {'SignalProcessing.Resampling.SamplingRate','SignalProcessing.FIRFilter.Frequencies', '', ...                
                'SignalProcessing.ICA.DataCleaning.DataSetting', '', 'SignalProcessing.ICA.Variant', ...
                'SignalProcessing.EpochExtraction','', ...
                'SignalProcessing.ERSPTransform','', ...
                'Prediction.FeatureExtraction.SpectralPrior', 'Prediction.FeatureExtraction.TemporalPrior', ...
                'Prediction.FeatureExtraction.SpatialPrior','Prediction.FeatureExtraction.AnatomicalPrior', '', ...
                'Prediction.MachineLearning.Learner.Lambdas', 'Prediction.MachineLearning.Learner.LossFunction','Prediction.MachineLearning.Learner.Regularizer'};
        end
    end
end


% TODO: reimplement priors

%         function [featuremodel,conditioningmodel,predictivemodel] = calibrate_prediction_function(self,varargin)
%             args = arg_define(varargin, ...
%                 arg_norep('signal'), ...
%                 arg_sub({'fex','FeatureExtraction'},{},...
%                     {arg({'spectral_prior','SpectralPrior'},@(f)1,[],'Spectral prior. Likelihood function of frequency in Hz.','guru',true), ...
%                      arg({'temporal_prior','TemporalPrior'},@(t)1,[],'Temporal prior. Likelihood function of time in s.','guru',true), ...
%                      arg({'spatial_prior','SpatialPrior'},@(p)1,[],'Spatial prior. Likelihood function of MNI coordinate vector.','guru',true), ...
%                      arg({'anatomical_prior','AnatomicalPrior'},{'Left Cerebrum','Right Cerebrum','Left Cerebellum','Right Cerebellum','Left Brainstem','Right Brainstem','Inter-Hemispheric'},{'-- Hemispheres --','Left Cerebrum','Right Cerebrum','Left Cerebellum','Right Cerebellum','Left Brainstem','Right Brainstem','Inter-Hemispheric','-- Lobes --','Anterior Lobe','Frontal Lobe','Frontal-Temporal Space','Limbic Lobe','Medulla','Midbrain','Occipital Lobe','Parietal Lobe','Pons','Posterior Lobe','Sub-lobar','Temporal Lobe','-- Gyri --','Angular Gyrus','Anterior Cingulate','Caudate','Cerebellar Lingual','Cerebellar Tonsil','Cingulate Gyrus','Claustrum','Culmen','Culmen of Vermis','Cuneus','Declive','Declive of Vermis','Extra-Nuclear','Fastigium','Fourth Ventricle','Fusiform Gyrus','Inferior Frontal Gyrus','Inferior Occipital Gyrus','Inferior Parietal Lobule','Inferior Semi-Lunar Lobule','Inferior Temporal Gyrus','Insula','Lateral Ventricle','Lentiform Nucleus','Lingual Gyrus','Medial Frontal Gyrus','Middle Frontal Gyrus','Middle Occipital Gyrus','Middle Temporal Gyrus','Nodule','Orbital Gyrus','Paracentral Lobule','Parahippocampal Gyrus','Postcentral Gyrus','Posterior Cingulate','Precentral Gyrus','Precuneus','Pyramis','Pyramis of Vermis','Rectal Gyrus','Subcallosal Gyrus','Sub-Gyral','Superior Frontal Gyrus','Superior Occipital Gyrus','Superior Parietal Lobule','Superior Temporal Gyrus','Supramarginal Gyrus','Thalamus','Third Ventricle','Transverse Temporal Gyrus','Tuber','Tuber of Vermis','Uncus','Uvula','Uvula of Vermis'}, 'Anatomical prior. Select anatomical structures that are likely to contain processes of interest.'),...
%                     }, 'Parameters for the feature-adaptation function. These parameters control how features are statistically adapted and extracted from the filtered data before they are passed int othe machine learning stage.'), ...
%                 arg_sub({'cond','Conditioning'},{},@self.feature_adapt_conditioning,'Feature conditioning parameters. Allows to further process features for better usability with classifiers.'), ...
%                 arg_sub({'ml','MachineLearning'},{'Learner',self.machine_learning_defaults()},@ml_train,'Machine learning stage of the paradigm. Operates on the feature vectors that are produced by the feature-extraction stage.'));
%             
%             try
%                 data = self.rssd_load_overcomplete(args.signal);
%                 % read out component dipole fits
%                 if ~isempty(data.dipfit)
%                     if isfield(data.dipfit,'multimodel')
%                         dipfits = [data.dipfit.multimodel{:}];
%                     else
%                         dipfits = data.dipfit.model;
%                     end
%                 else
%                     dipfits = [];
%                 end
%                 if ~isempty(args.fex.anatomical_prior) && ~isequal(args.fex.anatomical_prior,false) && ~isempty(dipfits)
%                     % if an anatomical prior was given, we can pre-prune the potenial ERSPs
%                     ersprange = [];                    
%                     for k=1:size(data.icaweights,1)
%                         matches{k} = intersect(dipfits(k).structures,args.fex.anatomical_prior);
%                         if ~isempty(matches{k})
%                             ersprange(end+1) = k; end
%                     end
%                 else
%                     % otherwise we need to compute them all
%                     ersprange = 1:size(data.icaweights,1);
%                 end
%                 
%                 % compute ERSPs
%                 disp('Now computing time/frequency decompositions...');
%                 [T,X,freqs,times] = evalc('self.ersp_compute(data,args,ersprange)'); %#ok<ASGLU>
%                 retain_ics = [];
%                 structures = {};
%                 probabilities = [];
%                 summed_probabilities = [];
%                 for k = ersprange
%                     % get left/right normalization vectors lhs/rhs
%                     lhs = args.fex.spectral_prior(freqs) .* sum(reshape(var(X{k},[],2),size(X{k},1),[]),2).^-0.25;
%                     rhs = args.fex.temporal_prior(times) .* sum(reshape(var(X{k},[],1),[],size(X{k},3)),2).^-0.25;
%                     % make sure that they are well-behaved
%                     lhs(~isfinite(lhs)) = 0; medlhs = median(lhs); lhs(lhs>3*medlhs) = medlhs;
%                     rhs(~isfinite(rhs)) = 0; medrhs = median(rhs); rhs(rhs>3*medrhs) = medrhs;
%                     % turn into a scaling matrix
%                     prior{k} =  diag(lhs) * ones(length(freqs),length(times)) * diag(rhs);
%                     % incorporate the spatial prior
%                     if ~isempty(dipfits)
%                         prior{k} = prior{k} * args.fex.spatial_prior(dipfits(k).posxyz); end
%                 end
%                 for k = ersprange
%                     % incorporate the anatomical prior
%                     if ~isempty(args.fex.anatomical_prior) && ~isequal(args.fex.anatomical_prior,false) && ~isempty(dipfits)
%                         [matches,idx] = intersect(dipfits(k).structures,args.fex.anatomical_prior); %#ok<ASGLU>
%                         % sum the probabilities for being in each of the accepted structures (can be > 1 as the structures are highly overlapping)
%                         prior{k} = sum(dipfits(k).probabilities(idx)) * prior{k};
%                         structures{k} = dipfits(k).structures(idx);
%                         probabilities{k} = dipfits(k).probabilities(idx);
%                         summed_probabilities(k) = sum(dipfits(k).probabilities(idx));
%                         % figure; topoplot(data.icawinv(:,k),data.chanlocs(data.icachansind),'electrodes','labels'); title([hlp_tostring(xdipfits(k).structures(idx)) ' - ' hlp_tostring(dipfits(k).probabilities(idx))]);
%                     end
%                     if ~all(all(prior{k}==0))
%                         retain_ics(end+1) = k; end
%                 end
%                 featuremodel.prior = prior;
%                 featuremodel.retain_ics = retain_ics;
%                 featuremodel.args = args;
%                 blocksizes = cellfun(@size,featuremodel.prior,'UniformOutput',false);
%                 
%                 % update machine learning parameters
%                 args.ml.learner.shape = vertcat(blocksizes{retain_ics}); %vertcat(blocksizes{cellfun(@prod,blocksizes)});
%                 
%         end        
%         
%         function features = feature_extract(self,signal,featuremodel)
%             try
%                 data = self.rssd_load_overcomplete(signal);
%                 % compute ERSP features
%                 [T,features] = evalc('self.ersp_compute(data,featuremodel.args,featuremodel.retain_ics)'); %#ok<ASGLU>
%                 % scale each component-ERSP by the respective prior & sparsify
%                 for c=featuremodel.retain_ics
%                     features{c} = bsxfun(@times,features{c},featuremodel.prior{c}); end
%                 % generate block-compressed version of the data
%                 features = reshape([features{featuremodel.retain_ics}],[],data.trials)';
%             catch
%                 disp('Glitch in para_rssd. halting...');
%                 keyboard;
%             end
%         end
