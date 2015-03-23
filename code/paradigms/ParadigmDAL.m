classdef ParadigmDAL < ParadigmDataflowSimplified
    % Advanced paradigm for slow cortical potentials and/or oscillatory processes.
    %
    % The DAL paradigm implements the framework described in [1], and can be customized to give
    % state-of-the-art results in almost any current BCI situation. Note that DAL is the name of the
    % optimization method [2], and not an accepted or recognized name for BCI paradigms using it (but is
    % used here for the lack of a better name). Since this method offers a lot of flexibility,
    % simplified versions are also available, restricting it to either slow cortical potentials
    % (DALERP) or to oscillatory processes (DALOSC).
    %
    % The DAL paradigm can use optional signal (pre-)processing stages, but blurs the boundary between
    % the subsequent feature extraction and machine learning stages. Features are defined using regions,
    % in time and frequency, of the data epoch (which must be chosen depending on the task at hand), and
    % an optimization criterion which includes a 'loss' term and a 'regularization' term, each defined
    % w.r.t. these features. The loss term can be chosen depending on whether the task is to predict
    % categories (for classifaction) or real values (for regression). The regularization term allows to
    % impose assumptions about the structure of the data (in the epoch), and, especially, in what
    % respect(s) it is assumed to be "simple" (to control for overfitting). Overfitting control
    % counter-balances the extreme flexibility of the model, which can essentially assign a different
    % weight for every sample and channel of the data. Regularization is very time-consuming (since a
    % parameter search is necessary), so that it may be practical to temporarily disable regularization
    % to get quick (yet suboptimal) results.
    %
    % The measure of simplicity depends on the type of data being supplied. For raw sensor-space data,
    % the assumption is that a) there are only relatively few persistent informative source
    % constellations over time, which create the patterns measured at sensor sites, and b) there are
    % relatively few informative temporal activity patterns across activated sources. This assumption is
    % implemented in the dual-spectral, 'ds', regularizer. For independent component data, the
    % assumption is usually that only few components are informative (which is the group sparsity with
    % row groups, 'glr', regularizer); the same regularizer can also be used to select the most
    % informative sensors. Another possible assumption is that only few channels/components, and only
    % few time-points (or, e.g. spectral components when operating on spectrally-transformed data) are
    % relevant; this is the 'l1' regularizer. Despite these powerful constraints, the method can be made
    % to overfit, especially when many different (possibly irrelevant) regions are configured. Due to
    % its ability to select the few most relevant structures in the data, the method can not just be
    % used to obtain predictive models, but also to investigate neuroscientific questions about the
    % underlying processes - for example, which areas of the brain (or collections of areas) are
    % information w.r.t. some known latent (condition) variable.
    %
    % The region selection includes some advanced parameters (order and norm) which can be customized
    % and/or parameter-searched if desired (in [1], a parameter search is recommended), but reasonable
    % defaults are assigned which cover a large fraction of use cases.
    %
    % Example 1: Consider the well-known BCI task of classifying motor imagery; the goal is to predict
    % whether the user is imagining a movement of his left or right hand. A calibration data set would
    % typically include events with types such as 'left-imag' and 'right-imag', to indicate the time and
    % type of instruction stimuli presented to the subject (see, e.g. [3]). The subject would then
    % imagine either one or the other movement over the course of several seconds. Several EEG features
    % can be used to infer the type of movement, including the readiness potential ([4], which is a slow
    % cortical potential) and event-related synchronization/desynchronization ([5], which are
    % oscillatory processes). Depending on the choice of frequency bands (here, separate regions are
    % used for the readiness potential (0.5-10Hz) and alpha (7-15Hz) and beta (15-30Hz) rhythms),
    % appropriate order and normalization coefficients are automatically chosen by para_dal, and a
    % combined feature matrix as in [1] is constructed, from which the relevant portions are
    % automatically selected. The split in alpha and beta rhythm is here for demonstration purposes; it
    % is not necessarily a superior choice in practice.
    %
    %   calib = io_loadset('data sets/john/gestures.eeg')
    %   myapproach = {'DAL', 'SignalProcessing',{'EpochExtraction',[0.5 3.5]}, ...
    %       'Prediction'{'FeatureExtraction',{'WindowFreqs',[0.5 10; 7 15; 15 30]}}}; 
    %   [loss,model,stats] = bci_train('Data',calib, 'Approach',myapproach, 'TargetMarkers',{'left-imag','right-imag'})
    %
    %
    % Example 2: Consider the working-memory load example from para_multiband_csp. Here, the assumption
    % is that at least the theta (4-6Hz) and alpha bands (7-15Hz) are relevant for predicting the number
    % of items held in working memory by the subject. In this case, we will not pose the problem as one
    % of classification, but rather as one of regression (where the target variable takes on values in
    % {1,2,3}).
    %
    %   calib = io_loadset('data sets/mary/nback.eeg')
    %   myapproach = {'DAL', 'SignalProcessing',{'EpochExtraction',[-0.5 2.5]}, ...
    %       'Prediction'{'FeatureExtraction',{'WindowFreqs',[4 6; 7 15]}, ...
    %                    'MachineLearning',{'Learner',{'dal' 'LossFunction','squared'}}}}; 
    %   [loss,model,stats] = bci_train('Data',calib, 'Approach',myapproach}, 'TargetMarkers',{'n1','n2','n3'});
    %
    %
    % Example 3: Another (hypothetical) scenario would be one in which complex temporal structure in
    % oscillatory and slow-changing processes is expected, for example during the visual processing of
    % faces versus houses (see, e.g., [6]). We assume a data set annotated with events of either the
    % 'face', the 'house', or the 'noise' type, indicating the presentation of a stimulus from the
    % respective category. It has been shown that there are significant differences in relevant ERP
    % components, as well as modulations in the 5-15Hz band and further effects, over the course to
    % 400ms. We first obtain a performance estimate, to find out to what degree these features can be
    % used for category prediction.
    %
    %   calib = io_loadset('data sets/john/facesvshouses.eeg')
    %   myapproach = {'DAL', 'SignalProcessing',{'EpochExtraction',[0 0.4]}, ...
    %       'Prediction'{'FeatureExtraction',{'WindowFreqs',[0.5 5; 5 15; 5 15; 5 15],'WindowTimes',[0 0.4; 0.05 0.15; 0.15 0.25; 0.25 0.4]}, ...
    %                    'MachineLearning',{'Learner',{'dal' 'LossFunction','squared'}}}}; 
    %   [loss,model,stats] = bci_train('Data',calib, 'Approach',myapproach, 'TargetMarkers',{'face','house','noise'})
    %
    % The resulting model can now be visualized using the techniques detailed in [1] to find the
    % relevant source projections and their temporal activations.
    %
    %
    % References:
    %  [1] Ryota Tomioka and Klaus-Robert Mueller, "A regularized discriminative framework for EEG analysis with application to brain-computer interface",
    %      Neuroimage, 49 (1) pp. 415-432, 2010.
    %  [2] Ryota Tomioka & Masashi Sugiyama, "Dual Augmented Lagrangian Method for Efficient Sparse Reconstruction",
    %      IEEE Signal Proccesing Letters, 16 (12) pp. 1067-1070, 2009.
    %  [3] Blankertz, B., Dornhege, G., Krauledat, M., M�ller, K., and Curio,G.
    %      "The non-invasive Berlin Brain-Computer interface: Fast acquisition of effective performance in untrained subjects."
    %      NeuroImage 37, 2 (Aug. 2007), 539?550.
    %  [4] Deecke, L.; Groezinger, B.; Kornhuber H.H. "Voluntary finger movement in man: Cerebral potentials and theory."
    %      Biol Cybern 23: 99?119, 1976
    %  [5] Pfurtscheller, G., and da Silva, L. "Event-related EEG/MEG synchronizaion and desynchronization: basic principles."
    %      Clin Neurophysiol 110 (1999), 1842-1857.
    %  [6] Guillaume A. Rousselet, Jesse S. Husk, Patrick J. Bennett and Allison B. Sekuler, "Single-trial EEG dynamics of object and face visual processing"
    %      NeuroImage 36 (3), 843-862, 2007
    %
    % Name:
    %   Dual-Augmented Lagrangian
    %
    %                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
    %                           2010-06-25

    
    methods
      
        function defaults = preprocessing_defaults(self)
            defaults = {'EpochExtraction',[0.5 3.5], 'Resampling',100};
        end
        
        function defaults = machine_learning_defaults(self)
            defaults = {'dal', 'Lambdas',2.^(4:-0.25:-3), 'NumFolds',5,'FoldMargin',1};
        end
        
        function [featuremodel,conditioningmodel,predictivemodel] = calibrate_prediction_function(self,varargin)
            args = arg_define(varargin, ...
                arg_norep('signal'), ...
                arg_sub({'fex','FeatureExtraction'},{},...
                    {arg_nogui('chan_prior'), ... 
                     arg_nogui('time_prior'), ... 
                     arg({'freqwnds','WindowFreqs'},[0.5 5; 7 30],[0 0.1 200 1000],'Frequency bands of interest. Matrix containing one row for the start and end of each frequency band on which DAL shall operate. Values in Hz.','cat','Feature Extraction'), ...
                     arg({'timewnds','WindowTimes'},[],[],'Time intervals of interest. Matrix containing one row for the start and end of the time window for each region. Values in seconds. If both this and the freqwnds parameter are non-empty, they should have the same number of rows.','cat','Feature Extraction'), ...
                     arg({'norms','WindowNorms'},[],[],'Normalization coefficients. One pair for each window, see [1] (Nx2 matrix, or []).','cat','Feature Extraction'), ...
                     arg({'orders','WindowOrders'},[],[],'Per-window order. This is the order (1 or 2) for each signal window (Nx1 matrix, or []).','cat','Feature Extraction'), ...
                     arg({'winfunc','WindowFunction'},'rect',{'barthann','bartlett','blackman','blackmanharris','bohman','cheb','flattop','gauss','hamming','hann','kaiser','nuttall','parzen','rect','taylor','triang','tukey'},'Type of window function. Typical choices are rect (rectangular), hann, gauss, blackman and kaiser; the same window function is assumed for all second-order windows (first-order windows use the rectangular window).'),...
                     arg({'winparam','WindowParameter','param'},[],[],'Parameter of the window function. This is mandatory for cheb, kaiser and tukey and optional for some others.','shape','scalar')}, 'Parameters for the feature-adaptation function. These parameters control how features are statistically adapted and extracted from the filtered data before they are passed int othe machine learning stage; the same window parameter is assumed for all second-order windows (first-order windows use the rectangular window).'), ...
                arg_sub({'cond','Conditioning'},{},@self.feature_adapt_conditioning,'Feature conditioning parameters. Allows to further process features for better usability with classifiers.'), ...
                arg_sub({'ml','MachineLearning'},{'Learner',self.machine_learning_defaults()},@ml_train,'Machine learning stage of the paradigm. Operates on the feature vectors that are produced by the feature-extraction stage.'));

            if ~isempty(args.fex.freqwnds) && ~isempty(args.fex.timewnds) && size(args.fex.freqwnds,1) ~= size(args.fex.timewnds,1)
                error('If both time and frequency windows are specified, both arrays must have the same number of rows (together they define the windows in time and frequency).'); end

            % shorten some names...
            [data,learner,twnds,fwnds,norms,orders,chan_prior,time_prior] = deal(args.signal,args.ml.learner,args.fex.timewnds,args.fex.freqwnds,args.fex.norms,args.fex.orders,args.fex.chan_prior,args.fex.time_prior);
            
            % find out whether we are using a dual-spectral regularizer or not
            dual_spectral = ~(any(strcmp(learner,'regularizer')) && ~any(strcmp(learner(find(strcmp(learner,'regularizer'),1,'last')+1),{'ds','dual-spectral'})));

            % replicate empty arrays
            if isempty(twnds)
                twnds = zeros(size(fwnds,1),0); end
            if isempty(fwnds)
                fwnds = zeros(size(twnds,1),0); end
            if isempty(norms)
                norms = zeros(size(fwnds,1),0); end
            if isempty(orders)
                orders = zeros(size(fwnds,1),0); end

            % for each region...
            for r = 1:max(size(fwnds,1),size(twnds,1))
                % start building a region descriptor
                reg.freq = arg_report('vals',@flt_spectrum,{'freq',fwnds(r,:)});
                % we assume the region to be first-order if any retained frequency is below 3, otherwise second-order
                reg.order = orders(r,:);
                if isempty(reg.order)
                    reg.order = quickif(any(reg.freq.freq<3),1,2); end
                % we only use the window function & parameter for second-order regions
                if reg.order == 2
                    reg.time = arg_report('vals',@flt_window,{'time',{twnds(r,:),args.fex.winfunc,args.fex.winparam}});
                else
                    reg.time = arg_report('vals',@flt_window,{'time',{twnds(r,:)}});
                end

                % squared data needs a different power than raw data
                reg.norm = norms(r,:);                
                if isempty(reg.norm)
                    reg.norm = quickif(reg.order==1,{-0.25,-0.25},{-0.5,-0.5}); end
                if isnumeric(reg.norm) && length(reg.norm)==2
                    reg.norm = {reg.norm(1),reg.norm(2)}; end
                
                % compute the filtered data ...
                flt = exp_eval_optimized(flt_spectrum('signal',flt_window('signal',data,reg.time),reg.freq));
                
                % and derive the preprocessing matrices from it
                if dual_spectral
                    X = num2cell(flt.data,[1 2]);
                    % in this case we use matrix powers of the inverse pooled covariance matrices for scaling
                    if reg.order == 1
                        % first-order: spatio-temporal scaling
                        reg.preproc = {hlp_diskcache('featuremodels',@cov_shrink,cat(2,X{:})')^reg.norm{1}, hlp_diskcache('featuremodels',@cov_shrink,cat(1,X{:}))^reg.norm{2}};
                        % apply spatio-temporal prior
                        if ~isempty(chan_prior)
                            reg.preproc{1} = reg.preproc{1} .* diag(chan_prior); end
                        if ~isempty(time_prior)
                            reg.preproc{2} = reg.preproc{2} .* diag(time_prior); end
                        reg.shape = [size(flt.data,1) size(flt.data,2)];
                    else
                        % second-order: just spatial scaling
                        tmp = cov(cat(2,X{:})');
                        reg.preproc = {tmp^reg.norm{1}, tmp^reg.norm{2}};
                        % apply spatial prior (temporal one does not apply)
                        if ~isempty(chan_prior)
                            reg.preproc{1} = reg.preproc{1} .* diag(chan_prior.^2);
                            reg.preproc{2} = reg.preproc{2} .* diag(chan_prior.^2);
                        end
                        if ~isempty(time_prior)
                            error('temporal prior does not apply in the case of a second-order detector'); end
                        reg.shape = [size(flt.data,1) size(flt.data,1)];
                    end
                else
                    X = flt.data;
                    % in this case we use per-channel / per-timepoint scaling matrices
                    if reg.order == 1
                        % first-order: spatio-temporal scaling
                        reg.preproc = {diag(mean(squeeze(var(X,[],2)),2))^-0.5, diag(mean(squeeze(var(X,[],1)),2))^-0.5};
                        % apply spatio-temporal prior
                        if ~isempty(chan_prior)
                            reg.preproc{1} = reg.preproc{1} .* diag(chan_prior); end
                        if ~isempty(time_prior)
                            reg.preproc{2} = reg.preproc{2} .* diag(time_prior); end
                        reg.shape = [size(flt.data,1) size(flt.data,2)];
                    else
                        % second-order: just spatial scaling
                        tmp = diag(mean(squeeze(var(X,[],2)),2))^-0.5;
                        reg.preproc = {tmp, tmp};
                        % apply spatial prior (temporal one does not apply)
                        if ~isempty(chan_prior)
                            reg.preproc{1} = reg.preproc{1} .* diag(chan_prior.^2);
                            reg.preproc{2} = reg.preproc{2} .* diag(chan_prior.^2);
                        end
                        if ~isempty(time_prior)
                            error('temporal prior does not apply in the case of a second-order detector'); end
                        reg.shape = [size(flt.data,1) size(flt.data,2)];
                    end
                end
                
                % compute & append the block scale of the region (we want to make sure that the
                % data is properly normalized so that our lambda range does not have to be tuned 
                % to fit the data scale)
                X = flt.data; [lhs,rhs] = reg.preproc{:};
                if reg.order == 1
                    Y = zeros(size(X));
                    for t=1:size(X,3)
                        Y(:,:,t) = lhs*X(:,:,t)*rhs; end
                else
                    Y = zeros(size(X,1),size(X,1),size(X,3));
                    for t=1:size(X,3)
                        Y(:,:,t) = lhs*cov(X(:,:,t)')*rhs; end
                end
                reg.preproc = [reg.preproc 1/sum(sum(sqrt(var(Y,[],3)))/(size(Y,1)*size(Y,2)))];
                
                % finally aggregate the region descriptors
                regions{r} = reg;                
            end
            vectorize = ~strcmp(learner.arg_selection,'dal');
            featuremodel = struct('regions',{regions}, 'vectorize',{vectorize});
            tmp = cellfun(@(x)x.shape,regions,'UniformOutput',false);
            % update machine learning parameters
            args.ml.learner.shape = vertcat(tmp{:});
            args.ml.learner.scaling = 'none';            
            
            % extract features & target labels
            features = self.feature_extract(args.signal, featuremodel);
            targets = set_gettarget(args.signal);

            % adapt and apply feature conditioning
            conditioningmodel = self.feature_adapt_conditioning('features',features,'targets',targets,args.cond);
            [features,targets] = self.feature_apply_conditioning(features,targets,conditioningmodel);
            
            % run the machine learning stage
            predictivemodel = ml_train('data',{features,targets}, args.ml);
        end        
        
        function features = feature_extract(self,signal,featuremodel)
            block = cell(1,length(featuremodel.regions));
            for r = 1:length(featuremodel.regions)
                reg = featuremodel.regions{r};
                % select a region
                flt = exp_eval_optimized(flt_spectrum('signal',flt_window('signal',signal,reg.time),reg.freq));
                % scale the features
                X = flt.data; [lhs,rhs,bs] = reg.preproc{:};
                if reg.order == 1
                    Y = zeros(size(X));
                    for t=1:size(X,3)
                        Y(:,:,t) = bs*(lhs*X(:,:,t)*rhs); end
                else
                    Y = zeros(size(X,1),size(X,1),size(X,3));
                    for t=1:size(X,3)
                        Y(:,:,t) = bs*(lhs*cov(X(:,:,t)')*rhs); end
                end
                % store, for later block-diagonalization
                block{r} = Y;
            end
            % combine all blocks into a block-compressed trial matrix
            features = reshape([block{:}],[],signal.trials)';
            features(~isfinite(features(:))) = 0;
            % and optionally vectorize the result
            if featuremodel.vectorize
                features = double(reshape(features,[],size(features,3))'); end
        end
        
        function layout = dialog_layout_defaults(self)
            layout = {'SignalProcessing.Resampling.SamplingRate', 'SignalProcessing.EpochExtraction', '', ...
                'Prediction.FeatureExtraction.WindowFreqs', 'Prediction.FeatureExtraction.WindowTimes', ...
                'Prediction.FeatureExtraction.WindowFunction', '', ...
                'Prediction.MachineLearning.Learner.Lambdas','Prediction.MachineLearning.Learner.LossFunction',...
                'Prediction.MachineLearning.Learner.Regularizer'};
        end        
    end
end
