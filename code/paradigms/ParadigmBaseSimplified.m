classdef ParadigmBaseSimplified < ParadigmBase
    % Simplified base class for BCI paradigms.
    %
    % Derive from this class if in your paradigm, calibration is done based on a single data set
    % (i.e. neither dataset collections nor stream bundles are handled meaningfully by it).
    % * you supply a calibrate_simple() function that takes a single recording and derives a model 
    %   from it (optionally specifying a filter chain for efficient pre-processing)
    % * you supply a predict_simple() function that takes a single data set and calculates outputs 
    %   from it
    % * (optionally, you may specify a visualization function)
    %
    %
    % If your paradigm can be understood as a dataflow model (i.e. a sequence of signal processing,
    % feature extraction, and machine learning), consider deriving from a simpler dataflow base class
    % instead (ParadigmDataflowSimplified).
    %
    % Notes:
    %   Most BCI paradigms (e.g., almost all from before 2008) do not support the full generality of
    %   possible inputs that are allowed by ParadigmBase (the general definition and implementation
    %   contract of a BCI paradigm, see ParadigmBase.m). Instead, firstly, most traditional BCI
    %   paradigms operate during online processing only on a single (multi-channel) signal instead
    %   of multiple concurrent signals, and secondly, most paradigms are initially calibrated based
    %   on a single calibration recording (usually a 30-60 minute recording of calibration data)
    %   instead of a (possibly tagged) collection of multiple data sets. For any such paradigm, the
    %   additional flexibility is useless and the corresponding organizational structure is an
    %   unnecessary implementation burden.
    %
    %   This base class is a refinement of ParadigmBase which omits the additional flexibility and
    %   provides an interface for simple paradigms. Any such paradigm should implement the
    %   calibrate_simple() method and the predict_simple() method, which are the simplified
    %   equivalents the calibrate() and predict() functions in ParadigmBase (without multiple
    %   concurrent streams and without multiple calibration data sets).
    %
    %   There exists a further specialization of this base class, ParadigmDataflowSimplified, which
    %   provides an additional layer of convenience for BCI paradigms that adhere to a standard 
    %   "dataflow" (multi-stage) structure. If a paradigm can be expressed as a sequence of signal 
    %   processing, feature extraction and machine learning, it is recommended to derive from this
    %   more specialized class instead (see ParadigmDataflowSimplified.m).
    %
    % Name:
    %   Simplified Base Paradigm (abstract)
    %
    %                            Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
    %                            2011-08-28
    
    methods
        
        function model = calibrate_simple(self,varargin)
            % Override this function to implement your model calibration code
            % Model = calibrate_simple(Signal,Arguments...)
            %
            % This function defines how a BCI is calibrated (i.e. adapted) to an individual based on
            % a given calibration recording. The calibration recording is a single EEGLAB data set
            % (though may contain arbitrary signal contents) annotated with target markers that
            % contain information about the desired BCI output at various time points in the
            % recording. Based on these markers and the data a relationship between signal contents
            % and the value of some (assumed) 'cognitive state variable' (to be approximated by the
            % BCI's outputs) is inferred by this function and captured in the predictive model 
            % returned by it (a struct with any parameters necessary for predict_simple).
            %
            % Implementations of this function should declare at least a parameter named "Signal",
            % plus any desired custom user options. These parameters should be declared using the
            % arg_define() facility to allow the framework to query the parameter names and defaults
            % (e.g. for GUI display). The returned model also describes any signal (pre-)processing that
            % shall be performed by BCILAB before the predict_simple() function is called by it,
            % by means of its .filter_graph field.
            %
            % In:
            %   Signal : a single continuous EEGLAB data set (annotated with target markers; see
            %            set_targetmarkers for more info)
            %
            %   Arguments : further optional user arguments
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
            %           The rest of the model struct is arbitrary and should contain whatever
            %           parameters are needed by predict_simple() to perform its operation.
            %
            % Notes:
            %   Overridden functions should declare their arguments using arg_define().
            
            error('This paradigm implements no calibrate_simple() function.');
        end
        
        
        function outputs = predict_simple(self,signal,model)
            % Override this function to implement your prediction code
            % Outputs = predict_simple(Signal,Model)
            %
            % This function defines how this BCI paradigm responds to "prediction" queries, i.e.
            % how, for a given piece of signal, an estimate of (or inference about) the user's
            % cognitive state shall be computed. It receives a segment of a pre-processed signal, as
            % well as a predictive model and shall produce an estimate about the most recent sample
            % in the given signal (although it may utilize any previous samples in its computation).
            % The signal is implicitly pre-processed before the function is called as defined by the
            % .filter_graph field of the predictive model (set by the calibrate_simple) function.
            % 
            % The signal may also be "epoched", i.e. contain multiple segments, in which case the 
            % output of this function shall be an array of predictions (one for each segment).
            %
            % In:
            %   Signal : a signal pre-processed according to the model's filter graph
            %
            %   Model : a predictive model as created by your calibrate_simple() function
            %
            % Out:
            %   Outputs : a prediction/estimate for the most recent time point in the data (or one for
            %             every epoch if the signal is epoched); see ml_predict for the allowed formats
            
            error('This paradigm implements no predict_simple() function.');
        end
        
        
        function visualize(self,model,varargin)
            % Optionally override this function to implement your visualization code
            % visualize(Model)
            %
            % In:
            %   Model: a model as created by your calibrate() function;
            %          a plot or GUI will be produced to inspect the model
            %
            
            error('This paradigm implements no visualize() function.');
        end
        
        
        
        % --- implementation ----
        
        function model = calibrate(self,varargin)
            % Implementation of the calibrate() function (see ParadigmBase)
            % Model = calibrate(Collection,GoalIdentifier,Arguments...)
            %
            % This is essentially a wrapper around calibrate_simple() which peels off any unused
            % organizational ballast (for multiple collections and stream bundles) before passing it
            % on. If multiple data sets are passed, they are implicitly merged.
            %
            % In:
            %   Collection : a collection (cell array) of stream bundles
            %
            %   GoalIdentifier: the goal-identifier
            %
            %   Arguments... : further optional user arguments
            %
            % Out:
            %   Model : The returned model struct
            
            % separate the collection argument from the rest
            collection = arg_extract(varargin,{'collection','Collection'},[],{});
            remove = find(strcmpi(varargin(1:2:end),'collection'))-1;
            retain = setdiff(1:length(varargin),[1 + 2*remove,2 + 2*remove]);
            varargin = varargin(retain);
            % also remove the goal identifier argument
            remove = find(strcmpi(varargin,'goal_identifier'));
            varargin([remove remove+1]) = [];
            % and extract the signal argument
            if ~isempty(collection)
                % check if any of the data sets in the collection has more than one stream
                if any(cellfun(@(x)length(x.streams) > 1,collection))
                    error('This paradigm does not support more than one stream in parallel.'); end
                if length(collection) == 1
                    signal = collection{1}.streams{1};
                else
                    % merge all data sets in the collection
                    tmp = cellfun(@(x)x.streams{1},collection,'UniformOutput',false);
                    signal = set_merge(tmp{:});
                end
            else
                signal = {};
            end
            
            % call calibrate_simple with the signal passed in as additional argument
            model = self.calibrate_simple('signal',signal,varargin{:});
        end
        
        
        function outputs = predict(self,bundle,model)
            % Override this function to implement your prediction code
            % Outputs = predict(Bundle,Model)
            %
            % This is a wrapper around the predict_simple() function which peels off the stream bundle
            % representation and just passes on the contained signal.
            %
            % In:
            %   Bundle : stream bundle preprocessed by the model's filter graph
            %
            %   Model : a predictive model
            %
            % Out:
            %   Outputs : the outputs of calibrate_simple()
            
            if length(bundle.streams) > 1
                error('This paradigm does not support more than one stream in parallel.'); end
            outputs = self.predict_simple(bundle.streams{1},model);
        end
        
    end
end


% (turn off a few editor warnings because some actual implementations are missing in this file)
%#ok<*INUSD,*STOUT,*MANU>
