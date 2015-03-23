classdef ParadigmWindowmeans < ParadigmDataflowSimplified
    % Standard paradigm for slow cortical potentials, using per-channel multi-window signal averages.
    %
    % The windowed means paradigm is a general method for capturing slow-changing cortical potentials,
    % most importantly in reaction to events (then called Event-Related Potentials / ERPs). It is
    % comprehensively described in [1]; The default parameters match one of its first applications, in
    % [2].
    %
    % The paradigm is implemented as a sequence of signal (pre-)processing, feature extraction and
    % machine learing stages. Signal processing usually includes spectral filtering (e.g., lowpass
    % filtering) and occasionally spatial filtering, either for dimensionality reduction (e.g., by
    % selecting channels) or for the extraction of sparsity, independence or other feature qualities
    % (e.g., via independent component analysis). The defining property of the paradigm is the feature
    % extraction, in which windowed averages of (pre-processed) signal data, per channel, are computed
    % and used as features for the subsequent machine learning stage. The dimensionality of the feature
    % space is (# of channels) x (# of windows), which can easily be high enough to exceed the
    % capabilities of simpler classifiers or lead to over-fitting. For these reasons, either very robust
    % classifiers need to be used (such as shrinkage LDA or logistic regression) or strong assumptions
    % must be imposed in the machine learning stage (e.g. sparsity or group sparsity), or the number of
    % windows and channels must be carefully controlled / optimized. The paradigm can also be applied to
    % spectral data, by the use of the fourier filter (in one of the non-complex modes), possibly in
    % combination with the data selection filter. A related paradigm is para_dal_lf, and its
    % generalization para_dal, both of which do not require explicitly specified windows, but can
    % operate on raw data by means of their powerful regularization. The paradigm usually requires a
    % fair amount of manual (or automatic) tuning, in which the optimal window boundaries are determined
    % based on task data. Another parameter that is usually adapted to the task is the length and
    % location of the data epoch under consideration.
    %
    % The paradigm is widely applicable to event-related slow-changing brain dynamics, including, for
    % example, the perception of self-induced errors [3], machine-induced errors and/or suprisal [4,5],
    % prediction of movement intent [2], or (c)overt attention. It can also be used to detect brain
    % processes without a preceding event (i.e. asynchronously) when sufficient amounts of data from the
    % 'nothing'/'rest' condition is included in the calibration data.
    %
    % Simple Example: Consider the goal of predicting whether a person perceives an event as being
    % erroneous (and possibly unexpected), or not. A typical calibration data set for this task would
    % cover a sequence of events, some erroneous, some not, and each event is encoded in the data as an
    % EEGLAB event with type 'err' or 'noerr'. According to the literature [6,7], the assumptions is
    % that these types of events should be reflected in the EEG as a slow cortical potential (e.g., the
    % f-ERN) within 250ms to 600ms following the event. An appropriate predictive model could be
    % obtained as follows:
    %
    %   calib = io_loadset('data sets/john/errors.eeg')
    %   myapproach = {'Windowmeans' 'SignalProcessing', {'EpochExtraction',[0 0.8],'SpectralSelection',[0.1 15]}, ...
    %                 'Prediction',{'FeatureExtraction',{'TimeWindows',[0.25 0.3; 0.3 0.35; 0.35 0.4; 0.4 0.45; 0.45 0.55; 0.55 0.6]}}};
    %   [loss,model,stats] = bci_train('Data',calib, 'Approach',myapproach, 'TargetMarkers',{'err','noerr'});
    %
    %
    % Complex Example: Consider the goal of anticipating a self-paced finger movement (for simplicity
    % only of one hand) of a person. A biological basis for this is the Bereitschaftspotential
    % (readiness potential, [8]). This is a difficult problems, since the detection should happen as
    % early as possible (especially before EMG onset), and because the detection should be reasonably
    % robust against false positives in an asynchronous setting. A possible calibration data set would
    % contain sporadic events in which the subject pressed a button ('press'), with periods of no
    % activity of varying length in between. Surrogate events will be placed in the data to mark epoch
    % windows of the two conditions 'no-press' and 'pre-press', using the function set_insert_markers.
    % Epochs will only be extracted for the surrogate events. Several pre-press data epochs will be
    % generated that end between 125ms to 100ms prior to each movement, and several no-press epochs will
    % be generated that lie well between any two movements. An IIR low-pass filter will be used due to
    % its low latency (replacing the paradigm's default FFT-based filter), and several fine-grained
    % windows will be placed at the very end (the "tip") of the epoch. In addition, several longer
    % "baseline" windows of different lengths will be placed in earlier parts of the epoch , to serve as
    % an adaptively chosen baseline (against which the tip of the epoch can be compared). Logistic
    % regression will be used as a classifier.
    %
    %   % load data with 'press' events
    %   calib = io_loadset('data sets/john/buttonpresses.eeg')
    %   % insert 7 'no-press' events safely between any two 'press' events
    %   calib = set_insert_markers(calib,'SegmentSpec',{'press' 3 -0.5 'press'}, 'Event','no-press', 'Count',7);
    %   % insert 7 'pre-press' events shortly before any 'press' event
    %   calib = set_insert_markers(calib,'SegmentSpec',{-0.125 -0.100 'press'}, 'Event','pre-press', 'Count',7);
    %   % define approach
    %   myapproach = {'Windowmeans' 'SignalProcessing', {'EpochExtraction',[-2 0],'SpectralSelection','off','IIRFilter',{[2.5 14],'lowpass'}}, ...
    %                 'Prediction',{'FeatureExtraction',{'TimeWindows',[-1.6 -0.5; -1.2 -0.5; -0.5 0.45; -0.2 -0.175; -0.025 0]}, ...
    %                               'MachineLearning',{'Learner',{'logreg'}}}};
    %   % learn a model
    %   [loss,model] = bci_train('Data',calib, 'Approach', myapproach, 'TargetMarkers',{'no-press','pre-press'})
    %
    %
    % References:
    %  [1] Benjamin Blankertz, Steven Lemm, Matthias Sebastian Treder, Stefan Haufe, and Klaus-Robert Mueller.
    %      "Single-trial analysis and classification of ERP components -- a tutorial."
    %       Neuroimage, 2010
    %  [2] Blankertz, B., Curio, G., Mueller, K.-R. "Classifying single trial EEG: towards brain computer interfacing."
    %      Adv Neural Inf Process Syst 14:157-164.
    %  [3] Benjamin Blankertz, Christin Sch�fer, Guido Dornhege, and Gabriel Curio.
    %      "Single Trial Detection of EEG Error Potentials: A Tool for Increasing BCI transmission rates"
    %  [4] Pierre W. Ferrez and Jose del R. Millan, "Error-Related EEG Potentials Generated during Simulated Brain-Computer Interaction",
    %      IEEE Trans. on Biomedical Engineering, 55(3):923-929, 2008
    %  [5] Zander T.O., Kothe C., Welke S., Roetting M. "Utilizing Secondary Input from Passive Brain-Computer Interfaces for Enhancing Human-Machine Interaction"
    %      In Hofmann A. (Ed.): Lecture Notes in Computer Science, Springer, Berlin Heidelberg, 2009.
    %  [6] Holroyd, C.B., Coles, M.G.. "The neural basis of human error processing: reinforcement learning, dopamine, and the error-related negativity"
    %      Psychological Review, 109, 679-709, 2002
    %  [7] Gehring, W.J., Coles, M.G.H., Meyer, D.E., Donchin, E.
    %      "The error-related negativity: an event-related brain potential accompanying errors."
    %      Psychophysiology 27, 34-41.
    %  [8] Deecke, L.; Groezinger, B.; Kornhuber H.H. "Voluntary finger movement in man: Cerebral potentials and theory."
    %      Biol Cybern 23: 99?119, 1976
    %
    % Name:
    %   Windowed Means
    %
    %                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
    %                               2010-04-29
    
    methods
      
        function defaults = preprocessing_defaults(self)
            defaults = {'SpectralSelection',[0.1 5],'EpochExtraction',[-1.28 0],'Resampling',100};
        end
                
        function model = feature_adapt(self,varargin)
            arg_define(varargin, ...
                arg_norep('signal'), ...
                arg({'wnds','TimeWindows'},[-0.15 -0.10;-0.10 -0.05;-0.05 0],[],'Epoch intervals to take as features. Matrix containing one row for the start and end of each time window over which the signal mean (per every channel) is taken as a feature. Values in seconds.','cat','Feature Extraction'), ...
                arg_nogui({'split_modalities','SplitModalities'},{},[],'Optionally separate out given modalities. Cell array of channel types. This passes a cell array of channel indices called modality_ranges down to the classifier. Can also be set to true to split all channel types present in the signal.','guru',true,'type','expression'), ...
                arg({'vectorize_features','VectorizeFeatures'},true,[],'Vectorize features. If disabled, sophisticated classifiers will recognize the channel/windows structure.'));
            model.wnds = wnds;
            model.vectorize_features = vectorize_features;
            model.chanlocs = signal.chanlocs;
            model.cov = cov(signal.data(:,:)');
            % allow passing per-modality channel ranges into the classifier (supported only by very few classifiers)
            if ~isempty(split_modalities)
                if isequal(split_modalities,true)
                    split_modalities = unique({signal.chanlocs.type}); end
                for m=1:length(split_modalities)
                    model.modality_ranges{m} = find(strcmp({signal.chanlocs.type},split_modalities{m})); end
            end
        end
        
        function features = feature_extract(self,signal,featuremodel)
            features = utl_picktimes(signal.data,(featuremodel.wnds-signal.xmin)*signal.srate);
            if featuremodel.vectorize_features
                features = reshape(features,[],size(signal.data,3))'; end
        end
        
        function visualize_model(self,varargin) %#ok<*INUSD>
            args = arg_define([0 3],varargin, ...
                arg_norep({'parent','Parent'},[],[],'Parent figure.'), ...
                arg_norep({'fmodel','FeatureModel'},[],[],'Feature model. This is the part of the model that describes the feature extraction.'), ...
                arg_norep({'pmodel','PredictiveModel'},[],[],'Predictive model. This is the part of the model that describes the predictive mapping.'), ...
                arg({'patterns','PlotPatterns'},false,[],'Plot patterns instead of filters. Whether to plot spatial patterns (forward projections) rather than spatial filters.'), ...
                arg({'paper','PaperFigure'},false,[],'Use paper-style font sizes. Whether to generate a plot with font sizes etc. adjusted for paper.'));
            arg_toworkspace(args);
            parent = args.parent;

            % no parent: create new figure
            if isempty(parent)
                myparent = figure('Name','Per-window weights'); end
            
            % number of pairs, and index of pattern per subplot
            np = size(fmodel.wnds,1);
            horz = ceil(sqrt(np));
            vert = ceil(np/horz);
            
            % get the weights
            if isfield(pmodel.model,'w')
                weights = pmodel.model.w;
            elseif isfield(pmodel.model,'W')
                weights = pmodel.model.W;
            elseif isfield(pmodel.model,'weights')
                weights = pmodel.model.weights;
            else
                error('Cannot find model weights.');
            end
            
            % check if weights contains a bias value
            if numel(weights)==length(fmodel.chanlocs)*np+1
                weights = weights(1:end-1);
            elseif numel(weights)~=length(fmodel.chanlocs)*np
                error('The model is probably not linear');
            end
            
            % turn into matrix, and optionally convert to forward projections
            weights = reshape(weights,length(fmodel.chanlocs),np);
            if patterns
                weights = fmodel.cov*weights;  end
            
            % display
            for p=1:np
                subplot(horz,vert,p,'Parent',parent);
                topoplot(weights(:,p),fmodel.chanlocs,'maplimits',[-max(abs(weights(:))) max(abs(weights(:)))]);
                t=title(['Window' num2str(p) ' (' num2str(fmodel.wnds(p,1)) 's to ' num2str(fmodel.wnds(p,2)) 's)']);
                if args.paper
                    set(t,'FontUnits','normalized');
                    set(t,'FontSize',0.5);                    
                    set(gca,'FontUnits','normalized');
                    set(gca,'FontSize',0.5);
                end
            end
        end
                
        function layout = dialog_layout_defaults(self)
            layout = {'SignalProcessing.Resampling.SamplingRate', 'SignalProcessing.EpochExtraction', ...
                'SignalProcessing.SpectralSelection.FrequencySpecification', '', ...
                'Prediction.FeatureExtraction.TimeWindows', '', 'Prediction.MachineLearning.Learner'};
        end
        
    end
end
                