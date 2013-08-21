classdef ParadigmBandpower < ParadigmDataflowSimplified
    % Basic paradigm for oscillatory processes, via per-channel logarithmic
    % bandpower (note: fairly primitive by modern standards)
    %
    % The logarithmic Bandpower estimates ("log-BP") paradigm is based on the design of the original
    % Graz Brain-Computer Interface [1][5], which used lateralized motor imagery for control. The 
    % features exploited by this paradigm in its original form are Event-Related Synchronization and
    % Desynchronization [2] localized in the motor cortex, but the paradigm is not restricted to these
    % applications. Similar measures have also been used in [4] although without machine learning.
    %
    % Generally, log-BP can be used as a simple method to operate on oscillatory processes, either in
    % relation to events, or asynchronously. The paradigm is simple in that it does not capture any
    % complex time variations in the oscillations detected, does not capture interactions between
    % multiple frequency bands, and does not contain data-adaptive signal processing or feature
    % extraction (the machine learning component is adaptive, though).
    %
    % The paradigm is implemented as a standard sequence of signal (pre-)processing (spatial/spectral
    % filtering), feature extraction, and machine learning. The defining property of the paradigm is
    % that it extracts, per trial, the per-channel log-variance log(var(X)) as features of the signal.
    % The resulting feature vectors are then passed along to the learner component. By default, the
    % paradigm uses a non-adaptive spatial filter, the surface Laplacian, and a non-adaptive spectral
    % filter. In place of these, a spatial transformation into Indepentent Components and/or adaptive
    % spectral filtering can be used. Since the feature space created by it is typically relatively
    % low-dimensional, almost any classifier in the toolbox can be used with log-BP, including
    % non-linear ones such as SVMperf and QDA. Aside from these optional modifications, the parameters
    % that almost certainly need to be tuned to the task at hand are the length of the data epoch and
    % the choice of a frequency band (defaulting to motor imagery time scales and frequency ranges),
    % both of which can also be found via a small parameter search.
    %
    % The most closely related (more advanced) paradigms are CSP (Common Spatial Patterns), FBCSP
    % (Filter-Bank CSP), and SpecCSP (Spectrally weighted CSP), which add adaptively learned spatial
    % filters in one (CSP) or more (FBCSP) bands, and adaptively learned spatio-spectral filters
    % (SpecCSP), respectively. An even more advanced paradigm for oscillatory processes is DALOSC
    % (named after the Dual-Augmented Lagrangian optimization used), which learns both spatial filters
    % and their relative weightings in a unified cost function (although omitting the log() transform).
    %
    % Some application areas include detection of major brain rhythm modulations (e.g. alpha, beta), for
    % example related to relaxation, aspects of workload, and motor imagery. See also [3].
    %
    % Simple Example: Consider a calibration data set in which a subject is in varying degrees of
    % relaxedness throughout different time periods. Markers 'relaxed', 'non-relaxed' and 'stressed'
    % have been inserted in 2 second intervals into these data ranges. The goal is to learn a model
    % which can predict these states, using alpha rhythm features.
    %
    %   calib = io_loadset('data sets/demos/relaxation.eeg')
    %   myapproach = {'ParadigmBandpower' 'SignalProcessing',{'FIRFilter',[6 8 14 15],'EpochExtraction',[-2 2]}};
    %   [loss,model,stats] = bci_train('Data',calib, 'Approach',myapproach, 'TargetMarkers',{'relaxed','non-relaxed','stressed'})
    %
    % Advanced Example: On the same data set, the relevant information shall be derived using a better
    % (adaptive) spatial filter (ICA) and a better classifier (logistic regression) than the defaults.
    % The classifier shall be sparse, since we assume that only few of the independent components in the
    % data will carry relevant information. A Bayesian variant of the classifier will be used (using
    % Automatic Relevance Determination) to avoid the need for time-consuming regularization. As an
    % alternative, the l1-regularized variant of the classifier could be used.
    %
    %   calib = io_loadset('data sets/demos/relaxation.eeg')
    %   myapproach = {'ParadigmBandpower', 'SignalProcessing',{'FIRFilter',[6 8 14 15],'EpochExtraction',[-2 2], ...
    %       'ICA','on'}, 'Prediction',{'MachineLearning',{'Learner',{'logreg','variant','vb-ard'}}}}
    %   [loss,model,stats] = bci_train('Data',calib, 'Approach',myapproach, 'TargetMarkers',{'relaxed','non-relaxed','stressed'})
    %
    % Examples:
    %   After an approach has been defined as in one of the following examples, a predictive model can be obtained
    %   (given a calibration data set and a specification of target markers) using bci_train:
    %   [loss,model,stats] = bci_train('Data',io_loadset('calibration_rec.eeg'),'Approach',myapproach','TargetMarkers',{'mymarker1','mymarker2'});
    %
    %   % define a log-bandpower approach, using the defaults (7-30 Hz bandpass, shrinkage LDA classifier)
    %   myapproach = 'Bandpower';
    %
    %   % define an approach using a different frequency filter (here ~10-15 Hz, i.e. alpha band)
    %   myapproach = {'Bandpower' 'SignalProcessing',{'FIRFilter',[9 10 14 15]}};
    %
    %   % define an approach using a different relevant time window around the markers (here 0.5s - 3s post-event)
    %   myapproach = {'Bandpower' 'SignalProcessing',{'EpochExtraction',[0.5 3]}};
    %
    %   % define an approach where the surface Laplacian filter is turned off
    %   myapproach = {'Bandpower' 'SignalProcessing',{'SurfaceLaplacian','off'}};
    %
    %   % turn off the FIRFilter filter, and use instead an IIRFilter to achieve the same effect
    %   myapproach = {'Bandpower' 'SignalProcessing',{'FIRFilter','off', 'IIRFilter',[5 7 25 30]}};
    %
    %   % use a different machine learning function (here: shrinkage QDA)
    %   myapproach = {'Bandpower' 'Prediction',{'MachineLearning',{'Learner','qda'}}};
    %
    %   % use a fast sparse Bayesian logistic regression
    %   myapproach = {'Bandpower' 'Prediction',{'MachineLearning',{'Learner',{'logreg',[],'variant','vb-ard'}}}};
    %
    %   % use a (slow!) l1-regularized logistic regression (using parameter search, takes 2-5 minutes)
    %   myapproach = {'Bandpower' 'Prediction',{'MachineLearning',{'Learner',{'logreg','lambda',search(2.^(-6:1.5:10)),'variant','l1'}}}};
    %
    % References:
    %  [1] G. Pfurtscheller, C. Neuper, "Motor imagery and direct brain-computer communication"
    %      Proceedings of the IEEE, Vol. 89, No. 7. (07 August 2002), pp. 1123-1134.
    %  [2] Pfurtscheller, G., and da Silva, L. "Event-related EEG/MEG synchronizaion and desynchronization: basic principles."
    %      Clin Neurophysiol 110 (1999), 1842-1857.
    %  [3] Buzsaki, G., "Rhythms of the brain"
    %      Oxford University Press US, 2006
    %  ﻿[4] D. J. McFarland, L. M. McCane, S. V. David, and J. R. Wolpaw,  “Spatial filter selection for {EEG-based} communication,” 
    %      Electroencephalography and Clinical Neurophysiology, vol. 103, no. 3, pp. 386-394, Sep. 1997.
    %  ﻿[5] J. Kalcher, D. Flotzinger, C. Neuper, S. Gölly, and G. Pfurtscheller, 
    %      “Graz brain-computer interface {II:} towards communication between humans and computers based on online classification of three different {EEG} patterns,” 
    %      Medical & Biological Engineering & Computing, vol. 34, no. 5, pp. 382-388, Sep. 1996.
    %
    % Name:
    %   log-Bandpower
    %
    %                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
    %                           2010-04-29
    
    methods
      
        function defaults = preprocessing_defaults(self)
            defaults = {'FIRFilter',{'Frequencies',[6 8 28 32],'Type','minimum-phase'}, 'SurfaceLaplacian','on', 'EpochExtraction',[0.5 3.5], 'Resampling',100};
        end
                
        function model = feature_adapt(self,varargin)
            arg_define(varargin, ...
                arg_norep({'signal','Signal'}));
            model.chanlocs = signal.chanlocs;
        end
        
        function features = feature_extract(self,signal,featuremodel)
            features = reshape(log(var(signal.data,0,2)),[size(signal.data,1),size(signal.data,3)])';
        end
        
        function visualize_model(self,parent,fmodel,pmodel) %#ok<*INUSD>
            % no parent: create new figure
            if isempty(parent)
                parent = figure('Name','Per-window weights'); end
            if isfield(pmodel.model,'w')
                weights = pmodel.model.w;
            elseif isfield(pmodel.model,'W')
                weights = pmodel.model.W;
            elseif isfield(pmodel.model,'weights')
                weights = pmodel.model.weights;
            else
                title('Cannot find model weights.');
                weights = [];
            end
            % display
            if ~isempty(weights)
                topoplot(weights,fmodel.chanlocs,'maplimits',[-max(abs(weights)) max(abs(weights))]);
            end
        end
        
        function layout = dialog_layout_defaults(self)
            layout = {'SignalProcessing.Resampling.SamplingRate', 'SignalProcessing.SurfaceLaplacian', ...
                'SignalProcessing.FIRFilter.Frequencies', 'SignalProcessing.FIRFilter.Type', ...
                'SignalProcessing.EpochExtraction', '', 'Prediction.MachineLearning.Learner'};
        end
        
    end
end

