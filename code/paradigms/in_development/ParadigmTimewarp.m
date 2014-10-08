classdef ParadigmTimewarp < ParadigmDataflowSimplified
    % Experimental paradigm that implements time-warping
    %
    % Name:
    %   Time-Warped Features
    %
    %                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
    %                               2014-06-19
    
    methods
      
        function defaults = preprocessing_defaults(self)
            defaults = {'IIRFilter',{[0.01 0.1],'highpass'},'EpochExtraction',[-5 0],'Resampling',100};
        end
            
        function defaults = machine_learning_defaults(self)
            % defaults = {'logreg','Variant',{'lars','ElasticMixing',0.5}};           % sparse classification
            % defaults = {'logreg','Variant',{'l2','Lambda',search(2.^(-4:1:4))}};  % l2 regularized logistic regression (to do: optimize range)
            % defaults = {'ridge','Lambda',search(2.^(-5:2:15))};                     % regression
            % defaults = {'rvm','Type','classification','Kernel','linear','KernelScale',1,'MaxIterations',100};                     % classification
            defaults = {'rvm','Type','regression','Kernel','linear','KernelScale',1,'MaxIterations',100};                     % classification
            % ........
        end
        
        function model = feature_adapt(self,varargin)
            args = arg_define(varargin, ...
                arg_norep('signal'), ...
                ... % TODO: put your other feature-extraction arguments here...
                arg({'vectorizeFeatures','VectorizeFeatures'},true,[],'Vectorize features. Needed for most of the basic classifiers.'));
            
            % this model struct is passed into feature_extract
            model.args = args;
        end
        
        function features = feature_extract(self,signal,featuremodel)
            % signal = epoched data, 5s up to marker X
            features = signal.data;  % for now use the raw data, not timewarped
            
            % TODO: timewarp features based on signal.event
            % features = timewarp(signal, featuremodel);
            
            % vectorize features for our classifier
            if featuremodel.args.vectorizeFeatures
                features = permute(features,[3 1:2 4:ndims(features)]);
                features = features(:,:);
            end
        end
    end
end
