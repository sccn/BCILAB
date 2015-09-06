classdef ParadigmDataflowSimplified < ParadigmBaseSimplified
    % Base class for simple dataflow-oriented BCI paradigms.
    %
    % This is the base class that almost all traditional (non-multi-subject and single-stream) BCI
    % paradigms should be derived from.
    %
    % Most such BCI paradigms do the following during online processing: The raw input signal is
    % first sent through a series of filter steps (i.e., signal processing), each here implemented
    % as a flt_***.m plugin function (e.g. spatial and spectral filters). To compute an estimate of
    % (or inference about) the cognitive state of the user at a particular time point in the signal,
    % a "prediction function" is invoked on a segment of that signal. It then typically first
    % extracts a set of features from the signal (yielding a feature vector or matrix) and then
    % passes on the feature vector to a machine learning plugin's prediction function (an
    % ml_predict***.m function) to obtain the desired output. The feature extraction and predictive
    % mapping can also be implemented as a custom combined step without going through a machine
    % learning plugin.
    %
    % The prediction function performs its processing with the help of a (separately computed)
    % "predictive model", which is a MATLAB struct that contains any parameters necessary for online
    % processing. A predictive model is computed in a dedicated calibration (a.k.a. 'training') step
    % based on a calibration recording. Such a recording is typically a 10-60 minute EEGLAB data set
    % which is annotated with target markers that indicate the desired output of the BCI at various
    % time points. During calibration, both the signal processing stages, the feature extraction
    % function, as well as the machine learning part may be adapted to optimally predict the desired
    % outputs for some given input signal data. In a dataflow model, calibration consists of first
    % applying a (custom-defined) sequence of pre-processing (i.e. filter) steps to the raw signal,
    % then optionally adapting any features based on the contents of the processed signal or
    % additional user parameters, and finally applying a machine learning plugin function (an
    % ml_train***.m function) to learn a statistical map from the (labeled) feature vectors (as
    % previously extracted from the data) onto the desired outputs (= labels). The feature
    % adaptation / extraction and machine learning can also be implemented in a signal combined step,
    % without the need to go through a separate plugin function.
    %
    % Derive from this class if your processing consists of a chain of filter steps followed by a
    % (possibly adaptive) feature extraction, followed by a machine learning step, and if in your
    % paradigm, calibration is done based on a single data set (i.e. neither dataset collections nor
    % stream bundles are handled meaningfully by it).
    %
    % * you supply a feature_adapt() function which takes a signal and arbitrary optional arguments
    %   and packages any information needed by the feature-extraction stage into a "feature model"
    %   struct (e.g. user parameters or signal-dependent parameters).
    % * you supply a feature_extract() function which takes a signal and a previously computed feature
    %   model and extracts features for every trial in the signal.
    % * typically, you specify the default signal processing stages to apply before feature extraction
    %   and the default machine learning step to apply after feature extraction (both of which can be
    %   overridden by users of the paradigm), by overriding preprocessing_defaults() and/or
    %   machine_learning_defaults()
    %
    % * optionally, you specify a visualize_model function for the resulting model and a default
    %   specification for a configuration dialog (dialog_layout_defaults, see other paradigms for examples).
    %
    % * if your learning and prediction process (after signal processing) is non-traditional (e.g. you
    %   perform a computation that blends traditional feature extraction and machine learning, or
    %   the notion of data points in your machine learning differs from the marker-locked trials in the
    %   calibration data), you may implement your own calibrate_prediction_function() and
    %   optionally apply_prediction_function() methods. In this case you may ignore the feature_adapt()
    %   and feature_extract() functions.
    %
    % Name:
    %   Basic Dataflow (abstract)
    %
    %                            Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
    %                            2011-08-28
    
    methods
        
        function defaults = preprocessing_defaults(self)
            % Optionally override this function to specify custom defaults for the preprocessing
            % pipeline (flt_pipeline).
            %
            % Out:
            %   Defaults : cell array of name-value pairs (or a struct) controlling the default settings
            %              for flt_pipeline (see this function for detailed information on the syntax).
            %
            % Notes:
            %   Here you can define which of the available signal processing plugins you want to
            %   apply to the raw input signal (e.g. FIR band-pass, epoch extraction, resampling),
            %   and what parameters should be used to call them. This is done by returning a cell
            %   array of the form {filtername, parameters, filtername, parameters, ...} where the
            %   filtername is either the name of a flt_*** function or its human-readable (and
            %   CamelCase) name; any filter declares such a name in a declare_properties(...) line
            %   at the beginning of the function. The CamelCase name is recommended here for
            %   consistency with the GUI display of these filters. The parameters of the filter are
            %   most generally a cell array of the function's inputs, or, if only one input is
            %   passed, may also be just the value itself (e.g. some numeric array). This is also
            %   the same way in which the filter stages to apply can be passed to flt_pipeline.m
            %   (which is the "dispatch function used internally by ParadigmDataflowSimplified,
            %   and which contains some additional documentation about the syntax allowed here).
            %
            %   The order in which filters are listed here is ignored, and instead usually
            %   auto-determined based on some common ordering rules (each filter declares some
            %   required and/or preferred ordering relationships); it may however be overridden partially
            %   or completely by appending a 'FilterOrdering',{'filtername','filtername','filtername', ...}
            %   name-value pair to the defaults cell array that is returned here.
            %
            %   The majority of paradigms perform some continuous signal processing (e.g. IIR
            %   frequency filtering), followed by extraction of epochs around target markers
            %   (implemented by the 'EpochExtraction' filter, see set_makepos.m), possibly followed
            %   by some further epoch-based signal processing (e.g. fourier transform). While it is
            %   possible to perform all processing on continuous data (in this case likely using a
            %   "target channel" instead of target markers), this is currently not a common practice
            %   in BCILAB, so it is more probable that the support for such processing is
            %   insufficiently documented and perhaps partially untested.
            %
            %   While this function defines the default signal processing setup, the user can
            %   partially or completely override the setup when he/she configures the paradigm into
            %   a custom approach (by passing a custom 'SignalProcessing' argument to the paradigm
            %   during calibration).
            
            % by default, here epochs of 1 second length will be extracted relative to the target
            % markers (beginning 0.5 seconds before each target marker, and ending 0.5 seconds after
            % it). See set_targetmarkers for an explanation of what a target marker is.
            defaults = {'EpochExtraction',{'TimeWindow',[-0.5 0.5]}};
        end
        
        
        function featuremodel = feature_adapt(self,varargin)
            % Override this function if you have a custom feature adaption step or need to declare
            % custom feature extraction parameters.
            % FeatureModel = feature_adapt(Signal, Arguments...)
            %
            % This function is called during calibration on the already pre-processed data. The
            % output of this function is a "featuremodel" struct that describes the parameters with
            % which feature extraction shall be done on the data for online processing (as well as
            % for the subsequent feature extraction during calibration, whose outputs feed into
            % machine learning). Any user parameters that shall go into feature extraction should be
            % declared here (using the arg_define() facility), as well as a parameter named 'Signal'
            % which is the pre-processed data. The signal may be used to adapt some of the feature-
            % extraction parameters in data-dependent (and possibly supervised, i.e. also
            % label-dependent) manner. See ParadigmCSP.m for a case with supervised feature adaptation
            % and ParadigmWindowmeans.m for a case without data-driven feature adaptation.
            %
            % The code below demonstrates how to declare a simple feature-extraction argument as
            % well as the (mandatory) signal argument, and how the arguments are usually documented
            % in paradigms. All of these arguments are accessible as sub-arguments of the
            % Prediction.FeatureExtraction argument of the paradigm.
            %
            % In:
            %   Signal : a signal as processed by the pre-processing pipeline (flt_pipeline);
            %            this is in almost all configurations an expoched data set
            %
            %   GroupInto : Feature grouping. This controls how the processed data in each epoch is
            %               vectorized into a feature vector; either the resulting vector has a
            %               block with all channels of the first sample in succession, followed by a
            %               block with all channels of the second sample and so on (channel
            %               grouping), or a block with all samples for the first channel, followed
            %               by a block with all samples of the second channel, and so on (samples
            %               grouping). The tensor setting will forward tensor-shaped features (for
            %               signals with arbitrary extra dimensions), and the vectorized setting
            %               will vectorize features, regardless of signal dimensionality. The matrix
            %               setting will use [channels x time-points] feature matrices and give an
            %               error if the signal has additional dimensions.
            %
            %               GroupInto is a sample user parameter; other feature extraction functions may
            %               define arbitrary other parameters (or none -- see some other paradigms for
            %               examples).
            %
            % Out:
            %   FeatureModel : an adapted model that can subsequently be used for feature extraction
            %                  (here: an empty struct)
            %
            % Notes:
            %   Overridden functions should declare their arguments using arg_define().
            
            arg_define(varargin, ...
                arg_norep({'signal','Signal'}), ...
                arg({'group_into','GroupInto'},'channels',{'channels','samples','matrix','tensor','vectorized'},'Feature grouping. This controls how the processed data in each epoch is vectorized into a feature vector; either the resulting vector has a block with all channels of the first sample in succession, followed by a block with all channels of the second sample and so on (channel grouping), or a block with all samples for the first channel, followed by a block with all samples of the second channel, and so on (samples grouping). The tensor setting will forward tensor-shaped features (for signals with arbitrary extra dimensions), and the vectorized setting will vectorize features, regardless of signal dimensionality.'));
            
            % pack our group_into parameter into the featuremodel for later use in feature_extract().
            featuremodel.group_into = group_into;
        end
        
        
        function features = feature_extract(self,signal,featuremodel)
            % Override this function to implement your feature extraction step.
            % Features = feature_extract(Signal, FeatureModel, Arguments...)
            %
            % This function implements the paradigm's feature extraction: given a pre-processed
            % signal and a (previously determined) featuremodel, return an array of feature vectors
            % (or in rare cases feature matrices) for subsequent use by the machine learning.
            %
            % This function is must return one feature vector for each target value in the signal.
            % Thus, if the signal is epoched, a feature vector must be returned for each epoch
            % (during online use, this function will usually only see a single epoch). If the signal
            % is continuous (i.e., no EpochExtraction was used in the filter setup), one feature
            % vector should be be returned for each sample in the signal (and - importantly - the
            % target channel, if any, must be ignored).
            %
            % In:
            %   Signal : a signal as processed by the preprocessing pipeline (flt_pipeline);
            %            this is in almost all configurations an epoched data set
            %
            %   FeatureModel : the adapted feature model, as specified by the FeatureAdaption step (if
            %                  any)
            %
            % Out:
            %   Features : extracted feature vectors (usually one per trial / target value in the data)
            %              these may be in any form supported by ml_train (and the default training function
            %              in particular), most frequently [#Trials x #Features]
            %
            
            switch featuremodel.group_into
                case 'channels'
                    % pass on feature vectors grouped by channels
                    if ndims(signal.data)>3
                        error('Your signal has extra dimensions: to use this type of signal, you need to set GroupInto either to ''tensor'' or to ''vectorized''.'); end
                    features = squeeze(reshape(signal.data,[],1,size(signal.data,3)))';
                case 'samples'
                    % pass on feature vectorss grouped by samples
                    if ndims(signal.data)>3
                        error('Your signal has extra dimensions: to use this type of signal, you need to set GroupInto either to ''tensor'' or to ''vectorized''.'); end
                    features = squeeze(reshape(permute(signal.data,[2 1 3]),[],1,size(signal.data,3)))';
                case 'matrix'
                    % pass on [CxTxN] feature matrices
                    if ndims(signal.data)>3
                        error('Your signal has extra dimensions: to use this type of signal, you need to set GroupInto either to ''tensor'' or to ''vectorized''.'); end
                    features = signal.data;
                case 'tensor'
                    % pass on [AxBxCx...xN] feature tensors
                    features = permute(signal.data,[1:2 4:ndims(signal.data) 3]);
                case 'vectorized'
                    % pass on [NxF] feature vectors
                    features = permute(signal.data,[3 1:2 4:ndims(signal.data)]);
                    features = features(:,:);
                otherwise
                    error('Unsupported setting for the GroupInto parameter: %s',featuremodel.group_into);
            end
        end
        
        
        function conditioningmodel = feature_adapt_conditioning(self,varargin)
            % Override this function if you need to implement feature conditioning steps.
            % ConditioningModel = feature_adapt_conditioning(Signal, Arguments...)
            %
            % This function serves to simplify feature representations so as to be better suited for
            % use with classifiers. This type of processing is not frequently used and the
            % default implementation currently handles only some basic cases.
            %
            % In:
            %   Features : feature representations, as allowed by ml_train
            %
            %   Targets : target-value representations, as allowed by ml_train
            %
            %   PruneTrivialFeatures : Prune trivial features. This prunes features which are constant
            %                          across the whole training set. Currently restricted to vectorized features. (default: false)
            %
            %   EqualizeClasses : ensure that the number of exemplars per class is equal by reducing trials.
            %                     (default: false)
            %
            % Out:
            %   ConditioningModel : a model that holds the parameters for feature conditioning.
            %
            % Notes:
            %   Overridden functions should declare their arguments using arg_define().
            
            args = arg_define(varargin, ...
                arg_norep({'features','Features'}), ...
                arg_norep({'targets','Targets'}), ...
                arg({'prune_trivial','PruneTrivialFeatures'},false,[],'Prune trivial features. This prunes features which are constant across the whole training set.'), ...
                arg({'equalize_classes','EqualizeClasses'},false,[],'Equalize class ratios. This removes trials of the larger class(es) such that the number of exemplars for all classes is equal.'));
            
            features = args.features;
            targets = args.targets;
            
            if args.prune_trivial
                conditioningmodel.prune_indices = find(sum(bsxfun(@minus,features(1,:),features))==0); end
            if args.equalize_classes
                classes = unique(targets);
                num_exemplars = sum(bsxfun(@eq,classes',targets));
                reduce_to = min(num_exemplars);
                subset = {};
                for k = 1:length(classes)
                    shuffled = shuffle(find(targets==classes(k)));
                    subset{k} = shuffled(1:reduce_to);
                end
                subset = sort(vertcat(subset{:}));
                conditioningmodel.subset = subset;
            end
            
            conditioningmodel.prune_trivial = args.prune_trivial;
            conditioningmodel.equalize_classes = args.equalize_classes;
        end
        
        function [features,targets] = feature_apply_conditioning(self,features,targets,conditioningmodel)
            % This function implements the actual feature conditioning.
            % [Features,Targets] = feature_apply_conditioning(Features,Targets,ConditioningModel)
            %
            % In:
            %   Features : feature representations, as allowed by ml_train
            %
            %   Targets : target-value representations, as allowed by ml_train (can be empty)
            %
            %   ConditioningModel : the model generated by feature_adapt_conditioning
            %
            % Out:
            %   Features : conditioned feature representations, as allowed by ml_train
            %
            %   Targets : conditioned target-value representations, as allowed by ml_train
            
            if conditioningmodel.prune_trivial
                features(:,conditioningmodel.prune_indices) = []; end
            if conditioningmodel.equalize_classes
                features = features(conditioningmodel.subset,:);
                if ~isempty(targets)
                    targets = targets(conditioningmodel.subset); end
            end
        end
        
        function defaults = machine_learning_defaults(self)
            % Optionally override this function to specify custom defaults for the machine learning
            % step (ml_train***.m).
            %
            % Similarly to preprocessing_defaults(), this function specifies the default settings
            % to use for machine learning. Usually, this involves selecting the machine learning plugin
            % to apply via its acronym (the *** in the respective ml_train***.m) and optionally
            % declaring the default user arguments to use for it, typically as name-value pairs.
            % These arguments will be processed by the function ml_train.m, which is the dispatch
            % function used internally by ParadigmDataflowSimplified.
            %
            % Out:
            %   Defaults : cell array of {Learner, Arguments...} for ml_train; where Learner is the shortcut
            %              name of the respective learning function (e.g. 'logreg' for ml_trainlogreg), and
            %              and the remaining elements are arguments that will be passed as user arguments
            %              into the respective machine learning function (usually as name-value pairs,
            %              see ParadigmWPI.m or ParadigmRSSD.m for examples).
            
            % by default, LDA (Linear Discriminant Analysis, ml_trainlda) is used with no arguments;
            % if called with no arguments, ml_trainlda will run "shrinkage LDA" which is a good default
            % linear classifier
            defaults = {'lda'};
        end
        
        
        function layout = dialog_layout_defaults(self)
            % Optionally override this function to specify a custom GUI dialog layout.
            %
            % Each BCI paradigm should ideally have a dialog that exposes its key user-configurable
            % settings (this dialog is brought by the GUI). In BCILAB, the dialog is auto-generated
            % from the arguments of calibrate_simple() (actually most generally calibrate()) and
            % their sub-arguments. These are rich argument declarations made via arg_define
            % (comprising both the names, default values, extended help texts and optionally valid
            % ranges of any parameter of the paradigm).
            %
            % This function lists a selection of arguments (by their declared CamelCase names) as a
            % cell array, in the order of appearance in the dialog (one argument typically
            % translates into a one-line label/widget combo, except if it has sub-arguments),
            % possibly interleaved with '' separators (translating into blank lines in the dialog for
            % optical separation). Sub-arguments of the top-level arguments in calibrate_simple() can
            % be accessed by dot notation (as in 'Prediction.FeatureExtraction.GroupInto).
            %
            % Out:
            %   Layout : config layout; This is a cell array of parameter names, optionally with '' interleaved.
            %            Each parameter name results in one or more lines (and corresponding entry
            %            fields) inserted into the GUI dialog for the respective parameter (multiple
            %            lines if the denoted parameter has sub-parameters). Sub-parameters of a
            %            parameter can be referred to by means of dot notation. '' entries generate
            %            blank lines (i.e. spacing) in the dialog. The variable names are those that
            %            are defined (and exposed) by the calibrate_simple() function.
            
            layout = {'SignalProcessing','','Prediction.FeatureExtraction','','Prediction.MachineLearning.Learner'};
        end
        
        
        function visualize_model(self,parent,featuremodel,predictivemodel) %#ok<*INUSD>
            % Optionally override this function to implement your visualization code
            % visualize(Parent,FeatureModel,PredictiveModel)
            %
            % In:
            %   Parent : parent window / figure
            %
            %   FeatureModel : a feature model as generated by feature_adapt and as understood by
            %                  feature_extract
            %
            %   PredictiveModel : a predictive model as generated by ml_train and as understood by
            %                     ml_predict
            
            text(0.5,0.5,'This paradigm does not yet implement a visualization.','HorizontalAlignment','center');
        end
        
        
        function tf = needs_voting(self)
            % Override this function if your feature extraction only works for 2-class data (e.g.,
            % like CSP). This is the case when the feature adaptation is done in a supervised manner
            % (i.e. using labels) and can only deal with two classes.
            tf = false;
        end
        
        
        function [featuremodel,conditioningmodel,predictivemodel] = calibrate_prediction_function(self,varargin)
            % Perform calibration of the prediction function; this includes everything except for signal
            % processing. This function can optionally be overridden if some custom feature-extraction /
            % machine learning data flow is desired; its user parameters may be arbitrarily redefined then.
            %
            % This function invokes the feature adaptation, feature extraction and machine learning
            % during the calibration phase (i.e. everything that is required to determine the
            % parameters of the BCI paradigm's prediction function).
            %
            % This function is what gives rise to the "Prediction" top-level argument of the paradigm;
            % as you see below, it has two sub-arguments: FeatureExtraction and MachineLearning, which
            % themselves are defined by feature_adapt() and ml_train().
            %
            % In:
            %   Signal : a signal as pre-processed according to the paradigm's pre-processing pipeline
            %
            %   FeatureExtraction : User parameters for the feature-extraction stage. These parameters
            %                       control how features are extracted from the filtered data before
            %                       they are passed int othe machine learning stage.
            %
            %   Conditioning : User parameters for an optional feature-conditioning stage. These parameters
            %                  control how features are remapped to features that are subsequently received
            %                  by the machine learning.
            %
            %   MachineLearning : Machine learning stage of the paradigm. Operates on the feature
            %                     vectors that are produced by the feature-extraction stage.
            %
            % Out:
            %   FeatureModel : a feature-extraction model as understood by apply_prediction_function()
            %                  or (if not otherwise customized) by the feature_extract() function
            %                  * special feature: if this contains a non-empty field named shape, this
            %                                     value will be passed on to the machine learning method
            %
            %   ConditioningModel : a model that is sandwiched between feature extraction and machine learning,
            %                       generated by feature_adapt_conditioning and understood by feature_apply_conditioning
            %
            %   PredictiveModel : a predictive model, as understood by apply_prediction_function() or
            %                     (if not otherwise customized) by the ml_predict() function
            %
            %
            % Notes:
            %   You may override this function if your prediction function blends traditional
            %   feature extraction and machine learning or otherwise makes this separation
            %   impractical (for example if you have an unusual mapping between training instances
            %   for machine learning and target values in the data set). This function should
            %   declare its arguments using arg_define().
            
            args = arg_define(varargin, ...
                arg_norep({'signal','Signal'}), ...
                arg_sub({'fex','FeatureExtraction'},{},@self.feature_adapt,'Parameters for the feature-adaptation function. These parameters control how features are statistically adapted and extracted from the filtered data before they are passed into the machine learning stage.'), ...
                arg_sub({'cond','Conditioning'},{},@self.feature_adapt_conditioning,'Feature conditioning parameters. Allows to further process features for better usability with classifiers.'), ...
                arg_sub({'ml','MachineLearning'},{'Learner',self.machine_learning_defaults()},@ml_train,'Machine learning stage of the paradigm. Operates on the feature vectors that are produced by the feature-extraction stage.'));
            
            % adapt features if necessary
            featuremodel = self.feature_adapt('signal',args.signal, args.fex);
            if isfield(featuremodel,'shape') && ~isempty(featuremodel.shape)
                % check if the learner supports a shape parameter...
                if isfield(args.ml.learner,'shape')
                    args.ml.learner.shape = featuremodel.shape; 
                else
                    warn_once('ParadigmDataflowSimplified:ignoring_shape','The learning function does not appear to support a shape parameter, but the paradigm prefers to supply one; ignoring the shape. This warning will not be shown again during this session.');
                end
            end
            if isfield(featuremodel,'modality_ranges') && ~isempty(featuremodel.modality_ranges)
                % check if the learner supports a modality_ranges parameter...
                if isfield(args.ml.learner,'modality_ranges')
                    args.ml.learner.modality_ranges = featuremodel.modality_ranges; 
                else
                    warn_once('ParadigmDataflowSimplified:ignoring_modality_ranges','The learning function does not appear to support a modality_ranges parameter, but the paradigm prefers to supply one; ignoring the modality_ranges. This warning will not be shown again during this session.');
                end
            end
            
            % try to extract some signal-related properties
            featuremodel.signalinfo.chanlocs = args.signal.chanlocs;
            featuremodel.signalinfo.chaninfo = args.signal.chaninfo;
            
            % extract features
            features = self.feature_extract(args.signal, featuremodel);
            
            % extract target labels
            targets = set_gettarget(args.signal);
            
            % adapt and apply feature conditioning
            conditioningmodel = self.feature_adapt_conditioning('features',features,'targets',targets,args.cond);
            [features,targets] = self.feature_apply_conditioning(features,targets,conditioningmodel);
            
            % run the machine learning stage
            predictivemodel = ml_train('data',{features,targets}, args.ml);
            
        end
        
        
        function outputs = apply_prediction_function(self,signal,featuremodel,conditioningmodel,predictivemodel)
            % Apply the feature extraction and final predictive mapping for every trial in the data
            % set (where a trial corresponds 1:1 to the outputs generated by set_gettarget() for the
            % given signal).
            %
            % This function basically is the BCI paradigm's prediction function - the only difference
            % is that the actual prediction function of ParadigmDataflowSimplified (below) may run
            % this function with different sets of parameters in a voting arrangement.
            %
            % In:
            %   Signal : a signal as pre-processed according to the paradigm's pre-processing pipeline
            %
            %   FeatureModel : a feature-extraction model as previously generated by
            %                  calibrate_prediction_function()
            %
            %   ConditioningModel : a feature-conditioning model as previously generated by
            %                       calibrate_prediction_function()
            %
            %   PredictiveModel : a predictive model, as previously generated by
            %                     calibrate_prediction_function()
            %
            % Out:
            %   Outputs : a prediction/estimate for the most recent time point in the data (or one for
            %             every epoch if the signal is epoched); see ml_predict for the allowed formats
            %
            % Notes:
            %   You may override this function if you have a predictive mapping that is not handled
            %   by any of the machine learning plugins (e.g. if you are using a custom computation
            %   in a custom calibrate_prediction_function()).
            
            % predict given the extracted features and the model
            features = self.feature_extract(signal, featuremodel);
            features = self.feature_apply_conditioning(features,[],conditioningmodel);
            outputs = ml_predict(features, predictivemodel);
        end
        
        
        % --- internal implementation ---
        
        function model = calibrate_simple(self,varargin)
            % Calibrates a BCI model based on the given signal and arguments, for later use by the
            % predict_simple function.
            % Model = calibrate_simple(Signal,Arguments...)
            %
            % This function defines the top-level arguments of the paradigm, SignalProcessing
            % (effectively
            %
            % In:
            %   Signal : a single continuous EEGLAB data set (usually annotated with target markers;
            %            see set_targetmarkers for more info)
            %
            %   SignalProcessing : optionally a cell array of custom signal processing parameters,
            %                      as {'filtername',{arguments...}, 'filtername',{arguments..}, ...}
            %                      where the filtername is the name of a flt_*** function or its
            %                      declared CamelCase name (declared by its respective .m file in its
            %                      declare_properties(...) line). CamelCase names are preferred as they
            %                      match the names displayed in the Review/Edit GUI. The arguments is
            %                      a cell array of arguments to the filter, usually name-value pairs.
            %
            %                      Note that these parameters are allowed to override the defaults
            %                      declared by the paradigm in its preprocessing_defaults() function.
            %                      The entire argument list is basically passed to and interpreted by
            %                      flt_pipeline() for execution.
            %
            %   Prediction : optionally a cell array of arguments to calibrate_prediction_function(),
            %                which determines, among others, how feature extraction and/or machine
            %                learning should be performed. This is a cell array of name-value pairs,
            %                and most standard paradigms (which don't override that function) have
            %                here a sub-argument called "FeatureExtraction" (of arguments to feature_adapt())
            %                and one called "MachineLearning" (of arguments to ml_train).
            %
            % Out:
            %   Model : a model struct with a mandatory field .filter_graph and arbitrary other content
            %           * The .filter_graph filed is a 1x1 cell array that contains the desciption of
            %             filter steps that is to be applied to the data, and is usually passed as either
            %             the .tracking.online_expression field of the processed Signal or the processed
            %             Signal itself. If no signal processing is performed by this paradigm, the raw
            %             signal may be passed.
            %
            %           * May have optional fields .prediction_function, .prediction_window and
            %             .prediction_channels - though these are generally auto-deduced
            %
            % Notes:
            %   If you find that you need to override this function (which should be very rare),
            %   it is a better choice to instead inherit directly from ParadigmBaseSimplified.
            
            args = arg_define(varargin, ...
                arg_norep({'signal','Signal'}), ...
                arg_sub({'flt','SignalProcessing'}, self.preprocessing_defaults(), @flt_pipeline, 'Signal processing stages. These parameters control filter stages that run on the signal level; they can be enabled, disabled and configured for the given paradigm. The prediction operates on the outputs of this stage.'), ...
                arg_sub({'pred','Prediction'}, {}, self.cached_method('calibrate_prediction_function'), 'Prediction stage. These parameters control the calibration and processing of stages that run on the output of the signal processing.'), ...
                arg({'arg_dialogsel','ConfigLayout'},self.dialog_layout_defaults(),[],'Parameters displayed in the config dialog. Cell array of parameter names to display (dot-notation allowed); blanks are translated into empty rows in the dialog. Each string refers to an argument (or structure thereof) of the paradigm''s calibrate function. If a structure is identified, all parameters of that struture are listed, except if it is a switchable structure - in this case, a pulldown menu with switch options is displayed.','type','cellstr','shape','row'));
            %@self.calibrate_prediction_function
            
            % first pre-process the data (symbolically)
            % this means that signal is turned into an unevaluated expression (data structure) like
            % flt_resample(flt_fir(signal,[7,30]), 200)
            signal_expression = flt_pipeline('signal',args.signal, args.flt); %#ok<*NODEF>
            
            % evaluate this in an optimized fashion (this effectively evaluates the filter expression
            % with some key optimizations, such as caching of intermediate results, turned on)
            signal = exp_eval_optimized(signal_expression);
            
            % with signal processing done, we now calibrate the prediction function on it, using
            % calibrate_prediction_function(). If the paradigm needs voting, we here call this function
            % on each pair of classes separately.
            model.args = rmfield(args,'signal');
            model.classes = unique(set_gettarget(signal),'rows');
            numclasses = size(model.classes,1);
            if numclasses > 2 && self.needs_voting()
                for i=1:numclasses
                    for j=i+1:numclasses
                        [model.voting{i,j}.featuremodel,model.voting{i,j}.conditioningmodel,model.voting{i,j}.predictivemodel] = self.calibrate_prediction_function('signal',exp_eval(set_picktrials(signal,'rank',{i,j})), args.pred); end
                end
            else
                [model.featuremodel,model.conditioningmodel,model.predictivemodel] = self.calibrate_prediction_function('signal',signal, args.pred);
            end
            
            model.tracking.filter_graph = signal;
        end
        
        
        function outputs = predict_simple(self,signal,model)
            % Override this function to implement your prediction code
            % Outputs = predict_simple(Sginal,Model)
            %
            % In:
            %   Signal : a signal pre-processed according to the model's filter graph
            %
            %   Model : a predictive model as created by your calibrate_simple() function
            %
            % Out:
            %   Outputs : a prediction/estimate for the most recent time point in the data (or one for
            %             every epoch if the signal is epoched); see ml_predict for the allowed formats
            %
            % Notes:
            %   If you find that you need to override this function (which should be very rare),
            %   it is a better choice to instead inherit directly from ParadigmBaseSimplified.
            
            if ~isfield(model,'voting')
                % predict given the extracted features and the model
                outputs = self.apply_prediction_function(signal,model.featuremodel,model.conditioningmodel,model.predictivemodel);
            else
                % 1-vs-1 voting is necessary, construct the aggregate result
                outputs = [];
                % vote, adding up the probabilities from each vote
                for i=1:length(model.classes)
                    for j=i+1:length(model.classes)
                        outcome = self.apply_prediction_function(signal,model.voting{i,j}.featuremodel,model.voting{i,j}.conditioningmodel,model.voting{i,j}.predictivemodel);
                        if isempty(outputs)
                            outputs = {'disc' , zeros(size(outcome{2},1),length(model.classes)), model.classes}; end
                        outputs{2}(:,[i j]) = outputs{2}(:,[i j]) + outcome{2};
                    end
                end
                
                % renormalize probabilities
                outputs{2} = outputs{2} ./ repmat(sum(outputs{2},2),1,size(outputs{2},2));
            end
        end
        
        
        function visualize(self,varargin)
            % Optionally override this function to implement your visualization code
            % visualize(Model)
            %
            % In:
            %   Model : a model as created by your calibrate() function;
            %           a plot or GUI will be produced to inspect the model
            %
            %   PlotOptions : cell array or struct of name-value pairs.

            args = arg_define(varargin, ...
                arg_norep({'model','Model'},struct(),[],'BCI Model to visualize.'), ...
                quickif(arg_supported(@self.visualize_model), ...
                    arg_sub({'plotopts','PlotOptions','plotoptions'},{}, @self.visualize_model, 'Plotting options.'), ...
                    arg({'plotopts','PlotOptions','plotoptions','options','Options'},{},[],'Plotting options. Cell array of name-value pairs.','type','expression')));
            
            if ~iscell(args.plotopts)
                args.plotopts = {args.plotopts}; end
            
            % visualize the model, either using one figure or multiple in case of voting
            if ~isfield(args.model,'voting')
                p = figure();
                self.visualize_model(p,args.model.featuremodel,args.model.predictivemodel,args.plotopts{:});
            else
                numcl = length(args.model.classes);
                numclx = length(args.model.classes)-1;
                for i=1:numcl
                    for j=i+1:numcl
                        p = figure('NumberTitle','off','MenuBar','none','Toolbar','none','Units','normalized', 'Name',sprintf('%d vs. %d',i,j), ...
                            'Position',[(i-0.9)/numclx (j-1-0.9)/numclx 0.8/numclx 0.8/numclx]);
                        self.visualize_model(p,args.model.voting{i,j}.featuremodel,args.model.voting{i,j}.predictivemodel,args.plotopts{:});
                    end
                end
            end
        end
        
    end
end

% disable a warning about self...
%#ok<*MANU>