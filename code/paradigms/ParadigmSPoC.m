classdef ParadigmSPoC < ParadigmDataflowSimplified
    %
    % The Source Power Comodulation (SPoC) paradigm [1] is used to predict continuous quantities (for
    % example, the workload level) from the power (amplitude) of oscillatory source processes in the
    % brain, in a particular frequency and time range of interest.
    %
    % The advantage of SPoC over approaches that operate either on channels or on sources that come
    % from beamforming or a component-based technique such as ICA is that it learns optimized
    % spatial filters in a supervised manner, similarly to how the Common Spatial Patterns algorithm
    % learns spatial filters for the case of binary class labels (in fact, SPoC can be viewed as a
    % generalization of CSP from binary labels to continuous labels).
    %
    % SPoC can be used with any regressor or classifier, although using ridge regression is the most
    % natural choice.
    %
    % Examples:
    %   After an approach has been defined as in one of the following examples, a predictive model can be obtained
    %   (given a calibration data set and a specification of target markers) using bci_train:
    %   [loss,model,stats] = bci_train('Data',io_loadset('calibration_rec.eeg'),'Approach',myapproach','TargetMarkers',{'mymarker1','mymarker2'});
    %
    %   % define a basic SPoC approach, using the defaults
    %   myapproach = 'SPoC';
    %
    %   % use a different frequency range (here: approx the alpha band)
    %   myapproach = {'SPoC' 'SignalProcessing',{'FIRFilter',[6 8 14 16]}};
    %
    %   % also use an IIR filter instead of the default FIR
    %   myapproach = {'SPoC' 'SignalProcessing',{'FIRFilter','off', 'IIRFilter',[5 7 14 17]}};
    %
    %   % use a time window relative to the target markers
    %   myapproach = {'SPoC' 'SignalProcessing',{'EpochExtraction',[-2 2]}};
    %
    %   % use sparse least-angle regression (LARS)
    %   myapproach = {'SPoC' 'Prediction',{'MachineLearning',{'Learner',{'logreg','Variant','lars','Regression',true}}}};
    %
    %   % use delay-embedding to learn spatio-spectral filters (4th order)
    %   myapproach = {'SPoC' 'SignalProcessing',{'DelayEmbed',4,'EpochExtraction',[-4 4]}};
    %
    %
    % References:
    %   [1] Dähne, S., Meinecke, F. C., Haufe, S., Höhne, J., Tangermann, M., Müller, K. R., & Nikulin, V. V.
    %       "SPoC: A novel framework for relating the amplitude of neuronal oscillations to behaviorally relevant parameters."
    %       NeuroImage 86 (2014), 111-122.
    %
    % Name:
    %   Source Power Comodulation
    %
    %                            Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
    %                            2013-11-17
    
    methods
      
        function defaults = preprocessing_defaults(self)
            % define the default pre-processing parameters of this paradigm
            defaults = {'FIRFilter',{'Frequencies',[6 8 28 32],'Type','minimum-phase'}, 'EpochExtraction',[0.5 3.5], 'Resampling',100};
        end
        
        function defaults = machine_learning_defaults(self)
            defaults = 'ridge';
        end                
                
        function model = feature_adapt(self,varargin)
            % adapt a feature representation using the SPoC algorithm
            arg_define(varargin, ...
                arg_norep('signal'), ...
                arg({'patterns','PatternPairs'},3,uint32([1 1 64 1000]),'Number of SPoC patterns (times two).','cat','Feature Extraction','type','expression','shape','row'), ...
                arg({'cov_lambda','CovarianceRegularization','lambda'}, 0.001, [0 1], 'Covariance regularization. This is a regularization parameter for the covariance estimate.'));
            
            if signal.nbchan < patterns
                error('SPoC requires at least as many channels as you request output patterns. Please reduce the number of pattern pairs.'); end
            
            weighted_cov = zeros(signal.nbchan);
            mean_cov = zeros(signal.nbchan);
            for k=signal.trials:-1:1
                weighted_cov = weighted_cov + self.cov_shrinkage(signal.data(:,:,k)',cov_lambda) * signal.epoch(k).target; 
                mean_cov = mean_cov + self.cov_shrinkage(signal.data(:,:,k)',cov_lambda);
            end
            [V,D] = eig(weighted_cov,mean_cov); %#ok<NASGU>
            
            model.filters = V(:,[1:patterns end-patterns+1:end]);
            P = inv(V);
            model.patterns = P([1:patterns end-patterns+1:end],:);
            model.chanlocs = signal.chanlocs;
        end
        
        function features = feature_extract(self,signal,featuremodel)
            % extract log-variance features according to SPoC
            features = zeros(size(signal.data,3),size(featuremodel.filters,2));
            for t=1:size(signal.data,3)
                features(t,:) = log(var(signal.data(:,:,t)'*featuremodel.filters)); end
        end               
        
        function visualize_model(self,parent,featuremodel,predictivemodel,varargin) %#ok<*INUSD>
            % number of pairs, and index of pattern per subplot
            np = size(featuremodel.patterns,1)/2; idx = [1:np 2*np:-1:np+1];
            % for each SPoC pattern...
            for p=1:np*2
                subplot(2,np,p,'Parent',parent);
                topoplot(featuremodel.patterns(idx(p),:),featuremodel.chanlocs);
                title(['SPoC Pattern ' num2str(idx(p))]);
            end
        end
        
        function layout = dialog_layout_defaults(self)
            % define the default configuration dialog layout 
            layout = {'SignalProcessing.Resampling.SamplingRate', 'SignalProcessing.FIRFilter.Frequencies', ...
                'SignalProcessing.FIRFilter.Type', 'SignalProcessing.EpochExtraction', '', ...
                'Prediction.FeatureExtraction.PatternPairs', '', 'Prediction.MachineLearning.Learner'};
        end
        
        function V = cov_shrinkage(self,X,lambda)
            V = cov(X);
            V = V*(1-lambda) + lambda*eye(length(V))*trace(V)/length(V);
        end        
    end
end

