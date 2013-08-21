classdef ParadigmBase
    % Base class for all BCI paradigms.
    %
    % Derive from this class if are looking for complete flexibility in your BCI structure: 
    %  * you supply a calibrate() function that can handle multiple recordings or multiple parallel 
    %    streams and computes a model from them (optionally also specifying a filter chain for 
    %    efficient pre-processing)
    %  * you supply a predict() function that can handle one or more parallel streams and generates 
    %    outputs from it (using information from the previously computed model)
    %  * (optionally, you may also supply a visualization function)
    %
    % If your paradigm is somewhat more traditional, consider deriving instead from one of the
    % simplified base classes. ParadigmBaseSimplified is an appropriate base class if your paradigm
    % handles only a single recording and stream for calibration (like most BCIs) but is not 
    % constrained in any other way. ParadigmDataflowSimplified is the most appropriate base class for 
    % almost all BCIs: takes a single recording and stream as input for calibration, operates online
    % on a single stream, and processes its data through a chain of signal processing plugins, 
    % followed by custom feature extraction, followed by one of the machine learning plugin.
    %
    % Notes:
    %   A BCI paradigm is a complete implementation of a particular type of BCI. The includes both
    %   the real-time processing performed on one or more input signals to produce outputs (e.g. the
    %   estimated value of a "cognitive state variable"), as well as any offline processing
    %   necessary for real-time operation. The offline processing basically amounts to the
    %   calibration of the parameters that are used during online processing (summarized and packed
    %   into what is called a "predictive model") based on one or more calibration recordings. A
    %   miscellaneous function associated with a BCI paradigm is the visualization of the parameters
    %   of a predictive model (function visualize()).
    %
    % Name:
    %   Base Paradigm (abstract)
    %
    %                            Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
    %                            2011-08-28
    
    methods
        
        function model = calibrate(self,varargin)
            % Override this function to implement your model calibration code.
            % Model = calibrate(Collection,GoalIdentifier,Arguments...)
            %
            % The calibrate() function takes one or more annotated recordings (of EEG or other
            % signals) plus any further auxilliary data and generates from these a predictive model
            % that contains all necessary parameters for the online processing code of this paradigm.
            %
            % In the simplest case the input is a single recording annotated with target markers
            % which encode the desired output of the BCI for various time points in the recording
            % (assuming that a statistical mapping can be learned from features of the recorded data
            % onto the desired outputs. The detailed interpretation of the markers is up to the paradigm.
            %
            % In more complex cases, the input can be a structured collection of recordings
            % taken under various conditions (e.g. different people or tasks) plus of a description
            % of the "target" person/task etc. at which the output of the calibration should be 
            % aimed. Auxilliary inputs can be additional measurements, e.g. of sensor properties or
            % individual physiology.
            %
            % A paradigm's implementation of this method should always declare a parameter named 
            % "Collection" and one named "GoalIdentifier", plus any additional user parameters or
            % hyper-parameters. Generally, any implementation should declare its arguments using the
            % arg_define() facility so that the parameter names and defaults can be queried by the 
            % framework (and presented in a GUI, for example).
            %
            % The output is generally a struct that contains any parameters needed by the predict()
            % function, plus a field (.filter_graph) that describes any signal (pre-) processing to
            % be performed on the input signals. This signal processing description is for the
            % framework and is (during online use or offline simulations) applied automatically to
            % the input signal(s) before the paradigm's predict() function would be invoked.
            %
            % In:
            %   Collection : a collection (cell array) of stream bundles
            %                * each stream bundle in the collection represents a different recording
            %                  (e.g. of a different person)
            %                * a stream bundle contains multiple signals overlapping in time; it is
            %                  represented as a struct with a field .streams (itself a cell array of EEGLAB
            %                  data set structs)
            %                * each bundle in the collection may have arbitrary distinguishing meta-data
            %                  fields besides the .stream field, such as subject, day, session, montage,
            %                  etc.; these allow organizing the learning problem to exploit the
            %                  associations between sets
            %
            %   GoalIdentifier: a struct that may contain identifying meta-data of the data to which 
            %                   the model should eventually be applied (i.e. subject, day, montage, etc.)
            %
            %   Arguments... : arbitrary further optional user arguments, defined in a base class
            %
            % Out:
            %   Model : a model struct with a mandatory field .tracking.filter_graph and arbitrary other
            %           contents
            %           * The filter_graph field is a cell array of filter descriptions: one per stream
            %             expected by the predict() function.
            %
            %             After filters have been applied to a data set, the resulting set contains the
            %             final state of the filters and a symbolic description of the filter chain that
            %             was applied to it in the field .tracking.online_expression. The filter_graph
            %             field is a cell array of one or more of these online epxressions.
            %
            %             For simplicity, it may also be given as a (processed or unprocessed) stream
            %             bundle.
            %
            %           * the model may have and optional fields .tracking.prediction_function,
            %             .tracking.prediction_window, and .tracking.prediction_channels (though these
            %             are normally auto-deduced, see utl_complete_model())
            %
            % Notes:
            %   Overridden functions should declare their arguments using arg_define().

            error('This paradigm implements no calibrate() function.');
        end
        
        
        function outputs = predict(self,bundle,model)
            % Override this function to implement your prediction code.
            % Outputs = predict(Bundle,Model)
            %
            % This function implements the BCI paradigm's response to a "prediction" query, i.e., it
            % receives a (pre-processed) signal segment (plus a previously computed predictive model
            % that contains all its other necessary parameters) and produces an estimate of (or
            % inference about) the value of some 'cognitive state variable' based on the given
            % signal's data. This output is "by contract" about the last sample in the segment (but
            % typically is computed based on at least some or all of the other samples in the
            % segment).
            % 
            % More generally, the input can be an epoched data set (i.e. multiple segments), in
            % which case an array of outputs (one per segment) shall be produced (mainly for better
            % efficiency during offline simulations); every paradigm should support this input form. 
            %
            % If supported by the paradigm (as encoded by calibrate() in the .filter_graph field
            % produced by it) the input may consist of multiple available signals that refer to the
            % same (or overlapping) time ranges; if multiple signals are passed, the .xmin/.xmax
            % properties of the signals determine the relative times of the contained samples (in
            % seconds). Most paradigms will not be that general. If multiple signals are supported,
            % it should be expected that the relative timing of input signals can be "off" depending
            % of varying real-world transmission latencies.
            %
            % Generally, whether one or more signals are passed, the input is formatted as a "stream
            % bundle", which is a struct with a field .streams, which in turn is a cell array of all
            % signals (a.k.a. streams) that are available.
            %
            % In:
            %   Bundle : a stream bundle preprocessed by the model's filter graph (i.e. one stream per
            %            cell entry in the model's .filter_graph)
            %
            %   Model : a predictive model as created by your calibrate() function
            %
            % Out:
            %   Outputs : a prediction/estimate for the most recent time point in the data (or one for
            %             every epoch if the bundle is epoched); see ml_predict for the allowed formats

            error('This paradigm implements no predict() function.');
        end
        
        
        function visualize(self,model,varargin)
            % Optionally override this function to implement your visualization code
            % visualize(Model)
            %
            % In:
            %   Model: a model as created by your calibrate() function;
            %          a plot or GUI will be produced to inspect the model
            
            error('This paradigm implements no visualize() function.');
        end
        
    end
end


% (turn off a few editor warnings because some actual implementations are missing in this file)
%#ok<*INUSD,*STOUT,*MANU>
