classdef ParadigmCSP < ParadigmDataflowSimplified
    % Common Spatial Pattern(s) algorithm. 
    %
    % The CSP paradigm is based on the design of the Berlin Brain-Computer Interface (BBCI) [1], more
    % comprehensively described in [2], which is mainly controlled by (sensori-)motor imagery. The
    % features exploited by this paradigm in its original form are Event-Related Synchronization and
    % Desynchronization [3] localized in the (sensori-)motor cortex, but the paradigm is not restricted
    % to these applications. CSP was originally introduced in [5] and first applied to EEG in [6].
    %
    % Due to its simplicity, speed and relative robustness, CSP is the bread-and-butter paradigm for
    % oscillatory processes, and if nothing else, can be used to get a quick estimate of whether the
    % data contains information of interest or not. Like para_bandpower, CSP uses log-variance features
    % over a single non-adapted frequency range (which may have multiple peaks), and neither temporal
    % structure (variations) in the signal is captured, nor are interactions between frequency bands.
    % The major strength of the paradigm is its adaptive spatial filter, which is computed using the CSP
    % algorithm.
    %
    % The paradigm is implemented as a standard sequence of signal (pre-)processing (spatial/spectral
    % filtering), feature extraction, and machine learning. The first preprocessing step is frequency
    % filtering, followed by an adaptively learned spatial filter (which is the defining propery of the
    % paradigm), followed by log-variance feature extraction and finally a (usually simple) machine
    % learning step applied to the log-variance features. The spatial filtering projects the channels of
    % the original signal down to a small set of (usually 4-6) surrogate channels, where the (linear)
    % mapping is optimized such that the variance in these channels is maximally informative w.r.t. to
    % the prediction task. The CSP filters can be obtained from the per-class signal covariance matrices
    % by solving a generalized eigenvalue problem (of the form [V,D]=eig(Cov1,Cov1+Cov2)). CSP can also
    % be applied to independent components to rate their importance or for better artifact robustness. A
    % wide range of classifiers can be used with CSP features, the most commonly used one being LDA.
    % There exists a large corpus of CSP variants and extensions, mostly to give better control over
    % spectral filtering, including multiband CSP (para_multiband_csp), Spectrally Weighted CSP
    % (para_speccsp), Invariant CSP, Common Spatio-Spectral Patterns (CSSP), Common Sparse Spectral
    % Spatial Pattern (CSSSP), Regularized CSP, and several others. A more advanced (but also
    % computationally more costly) paradigm than CSP is the Dual-Augmented Lagrange Paradigm
    % (para_dal/para_dal_hf). The length of the data epoch and the choice of a frequency band
    % (defaulting to motor imagery time scales and frequency ranges) are the parameters that are most
    % commonly tuned to the task, both of which can also be found via a small parameter search.
    %
    % Some application areas include detection of major brain rhythm modulations (e.g. alpha, beta), for
    % example related to relaxation/stress, aspects of workload, sensori-motor imagery, visual
    % processing vs. idling and other idle-rhythm-related questions, or emotion recognition. See also
    % [4].
    %
    % Examples:
    %   After an approach has been defined as in one of the following examples, a predictive model can be obtained
    %   (given a calibration data set and a specification of target markers) using bci_train:
    %   [loss,model,stats] = bci_train('Data',io_loadset('calibration_rec.eeg'),'Approach',myapproach','TargetMarkers',{'mymarker1','mymarker2'});
    %
    %   % define a basic CSP approach, using the defaults (7-30 Hz bandpass, shrinkage LDA classifier, 3 pattern pairs)
    %   myapproach = 'CSP';
    %
    %   % use an FIR filter restricted to the alpha band
    %   myapproach = {'CSP' 'SignalProcessing',{'FIRFilter',[6 8 14 15]}};
    %
    %   % use an IIR filter instead of the default FIR
    %   myapproach = {'CSP' 'SignalProcessing',{'FIRFilter','off', 'IIRFilter',[5 7 14 17]}};
    %
    %   % also restrict the model to a stationary subspace
    %   myapproach = {'CSP' 'SignalProcessing',{'FIRFilter','off', 'IIRFilter',[5 7 14 17], 'StationarySubspace',{'StationaryDim',-0.1,'Operation','keep_stationary'}}};
    %
    %   % use a sharp FFT band-pass filter
    %   myapproach = {'CSP' 'SignalProcessing',{'FIRFilter','off','SpectralSelection',[7 15]}};
    %
    %   % use a simple logistic regression classifier (variational Bayes) instead of the LDA
    %   myapproach = {'CSP' 'Prediction',{'MachineLearning',{'Learner','logreg'}}};
    %
    %   % use a simple logistic regression classifier (sparse variational Bayes)
    %   myapproach = {'CSP' 'Prediction',{'MachineLearning',{'Learner',{'logreg','variant','vb-ard'}}}};
    %
    %   % use the sparse logistic regression classifier but applied to a larger set of patterns
    %   myapproach = {'CSP' 'Prediction',{'MachineLearning',{'Learner',{'logreg','variant','vb-ard'}},'FeatureExtraction',{'PatternPairs',6}}};
    %
    %   % using quadratic discriminant analysis
    %   myapproach = {'CSP' 'Prediction',{'MachineLearning',{'Learner','qda'}}};
    %
    %   % using Gaussian mixture models (variational Bayesian Dirichlet process prior)
    %   myapproach = {'CSP' 'Prediction',{'MachineLearning',{'Learner','gmm'}}};
    %
    %   % using relevance vector machines (here with a fixed kernel scale for speed)
    %   myapproach = {'CSP' 'Prediction',{'MachineLearning',{'Learner',{'rvm','kernel','rbf','gamma',1}}}};
    %
    %   % use l1-regularized logreg (which involves a parameter search over the regularization parameter); takes about 2 minutes
    %   myapproach = {'CSP' 'Prediction',{'MachineLearning',{'Learner',{'logreg','variant','l1','lambda',search(2.^(-6:1.5:10))}},'FeatureExtraction',{'PatternPairs',6}}};
    %
    %   % use support vector machines (using the SVMlight package); note: also requires a parameter search (takes about 2 minutes)
    %   myapproach = {'CSP' 'Prediction',{'MachineLearning',{'Learner',{'svmlight','cost',search(2.^(-5:2.5:15)),'gamma',1}}}};
    %
    %   % optimize the size of the stationary subspace to retain; note: cross-validation takes place on continuous data (takes approx 2.5 minutes)
    %   myapproach = {'CSP' 'SignalProcessing',{'FIRFilter','off', 'IIRFilter',[5 7 14 17], 'StationarySubspace',{'StationaryDim',search(-0.5:0.1:-0.1),'Operation','keep_stationary'}}};
    %
    %   % optimize the location of the frequency band manually (note: 2-dimensional parameter space -- takes approx. 7 minutes)
    %   myapproach = {'CSP' 'SignalProcessing',{'FIRFilter','off', 'SpectralSelection',[search(7:10) search(14:2:30)]}};
    %
    % References:
    %  [1] Blankertz, B., Dornhege, G., Krauledat, M., M�ller, K., and Curio,G.
    %      "The non-invasive Berlin Brain-Computer interface: Fast acquisition of effective performance in untrained subjects."
    %      NeuroImage 37, 2 (Aug. 2007), 539?550.
    %  [2] Benjamin Blankertz, Ryota Tomioka, Steven Lemm, Motoaki Kawanabe, and Klaus-Robert Mueller.
    %      "Optimizing spatial filters for robust EEG single-trial analysis."
    %      IEEE Signal Process Mag, 25(1):41-56, January 2008
    %  [3] Pfurtscheller, G., and da Silva, L. "Event-related EEG/MEG synchronizaion and desynchronization: basic principles."
    %      Clin Neurophysiol 110 (1999), 1842-1857.
    %  [4] Buzsaki, G., "Rhythms of the brain"
    %      Oxford University Press US, 2006
    %  [5] Fukunaga K., "Introduction to Statistical Pattern Recognition"
    %      Academic Press, Computer Science and Scientific Computing Series, 1990
    %  [6] Ramoser, H., Gerking, M., and Pfurtscheller, G. "Optimal spatial filtering of single trial EEG during imagined hand movement."
    %      IEEE Trans. Rehab. Eng 8 (2000), 446, 441.
    %
    % Name:
    %   Common Spatial Patterns
    %
    %                            Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
    %                            2010-04-29
    
    
    methods
      
        function defaults = preprocessing_defaults(self)
            % define the default pre-processing parameters of this paradigm
            defaults = {'FIRFilter',{'Frequencies',[6 8 28 32],'Type','minimum-phase'}, 'EpochExtraction',[0.5 3.5], 'Resampling',100};
        end
                
        function model = feature_adapt(self,varargin)
            % adapt a feature representation using the CSP algorithm
            arg_define(varargin, ...
                arg_norep('signal'), ...
                arg({'patterns','PatternPairs'},3,uint32([1 1 64 10000]),'Number of CSP patterns (times two).','cat','Feature Extraction'),...
                arg({'dologtransform','LogTransform','logtransform'},true,[],'Perform log-transform. This is almost always the right thing to do.'),...
                arg({'shrinkage_cov','ShrinkageCovariance','ShrinkageCov'},false,[],'Shrinkage covariance estimator. Whether to use shrinkage to estimate the covariance matrices.'));
            
            if signal.nbchan < patterns
                error('CSP requires at least as many channels as you request output patterns. Please reduce the number of pattern pairs.'); end
            for k=1:2
                trials{k} = exp_eval(set_picktrials(signal,'rank',k));
                if shrinkage_cov
                    covar{k} = hlp_diskcache('featuremodels',@cov_shrink,reshape(trials{k}.data,size(trials{k}.data,1),[])');
                else
                    covar{k} = cov(reshape(trials{k}.data,size(trials{k}.data,1),[])');
                end
                covar{k}(~isfinite(covar{k})) = 0;
            end
            [V,D] = eig(covar{1},covar{1}+covar{2}); %#ok<NASGU>
            model.filters = V(:,[1:patterns end-patterns+1:end]);
            P = inv(V);
            model.patterns = P([1:patterns end-patterns+1:end],:);
            model.cov = cov(signal.data(:,:)');
            model.chanlocs = signal.chanlocs;
            model.logtransform = dologtransform;
        end
        
        function features = feature_extract(self,signal,featuremodel)
            % extract log-variance features according to CSP
            features = zeros(size(signal.data,3),size(featuremodel.filters,2));
            for t=1:size(signal.data,3)
                features(t,:) = var(signal.data(:,:,t)'*featuremodel.filters); end
            if featuremodel.logtransform
                features = log(features); end
        end
        
        function visualize_model(self,varargin) %#ok<*INUSD>
            args = arg_define([0 3],varargin, ...
                arg_norep({'myparent','Parent'},[],[],'Parent figure.'), ...
                arg_norep({'featuremodel','FeatureModel'},[],[],'Feature model. This is the part of the model that describes the feature extraction.'), ...
                arg_norep({'predictivemodel','PredictiveModel'},[],[],'Predictive model. This is the part of the model that describes the predictive mapping.'), ...
                arg({'patterns','PlotPatterns'},true,[],'Plot patterns instead of filters. Whether to plot spatial patterns (forward projections) rather than spatial filters.'), ...
                arg({'paper','PaperFigure'},false,[],'Use paper-style font sizes. Whether to generate a plot with font sizes etc. adjusted for paper.'), ...
                arg_nogui({'nosedir_override','NoseDirectionOverride'},'',{'','+X','+Y','-X','-Y'},'Override nose direction.'));
            arg_toworkspace(args);
            
            % determine nose direction for EEGLAB graphics
            try
                nosedir = args.fmodel.signal.info.chaninfo.nosedir;
            catch
                disp_once('Nose direction for plotting not store in model; assuming +X');
                nosedir = '+X';
            end
            if ~isempty(nosedir_override)
                nosedir = nosedir_override; end            
            % number of pairs, and index of pattern per subplot
            np = size(featuremodel.patterns,1)/2; 
            idx = [1:np 2*np:-1:np+1];            
            % for each CSP pattern...
            for p=1:np*2
                subplot(2,np,p,'Parent',myparent);
                if args.patterns
                    plotdata = featuremodel.patterns(idx(p),:);
                else
                    plotdata = featuremodel.filters(:,idx(p));
                end
                topoplot(plotdata,featuremodel.chanlocs,'nosedir',nosedir);
                t = title(['CSP Pattern ' num2str(idx(p))]);
                if args.paper
                    set(t,'FontUnits','normalized');
                    set(t,'FontSize',0.1);                    
                end
            end
        end
        
        function layout = dialog_layout_defaults(self)
            % define the default configuration dialog layout 
            layout = {'SignalProcessing.Resampling.SamplingRate', 'SignalProcessing.FIRFilter.Frequencies', ...
                'SignalProcessing.FIRFilter.Type', 'SignalProcessing.EpochExtraction', '', ...
                'Prediction.FeatureExtraction.PatternPairs', '', 'Prediction.MachineLearning.Learner'};
        end
        
        function tf = needs_voting(self)
            % standard CSP requires voting to handle more than 2 classes
            tf = true; 
        end
        
    end
end

