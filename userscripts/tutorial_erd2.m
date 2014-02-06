% --- tutorial script for spectral analysis ---

% This tutorial continues the analysis of motor cortex neural idle oscillations, this time
% with more customized scripting examples.
%
% The used data set contains a sequence of trials in which a subject was instructed to imagine
% moving either the left hand, or the right hand, or a foot. The data contains no markers, but 
% a trigger channel whose level indicates the type of condition. Also, there are multiple sessions 
% that should be processed jointly.
%
% Data courtesy of Yijun Wang, SCCN.
%
% Citation: 
%  Y. Wang, B. Hong, X. Gao, and S. Gao, "Implementation of a brain-computer interface based on three states of motor imagery." 
%  Proceedings of 29th International IEEE EMBS Conference, Lyon, France, pp. 5059-5062, 23-26 Aug, (2007).
%
%#ok<*SUSENS,*NOPTS,*ASGLU,*NASGU,*SNASGU,*SAGROW,*SNAFU> % turn off a few editor warnings...

%% --- load training and test data (invoke this first) ---

% load and curate raw data for each session (BioSemi format)
for s = 1:4
    filename = ['data:/tutorial/imag_movements2/session' num2str(s) '.bdf'];
    % load & retain the first 32 channels
    session{s} = exp_eval_optimized(io_loadset(filename,'channels',1:32));
    % override the channel locations (they are incorrect in the original data)
    session{s}.chanlocs = set_infer_chanlocs('data:/tutorial/imag_movements2/mi32.loc');
end

% concatenate the first 3 sessions as training data (leave the 4th session for testing)
traindata = set_merge(session{1:3});

%% --- train 3-class model using the Spec-CSP method ---

% define approach (here Spec-CSP is restricted to the alpha band)
myapproach = {'SpecCSP' 'Prediction',{'FeatureExtraction',{'SpectralPrior',[7 15]}}};

% learn model, across 3 classes (0=foot, 16=left, 30=right)
[trainloss,lastmodel,laststats] = bci_train('Data',traindata,'Approach',myapproach,'TargetMarkers',{'0','16','30'});
disp(['training mis-classification rate: ' num2str(trainloss*100,3) '%']);

% visualize results
bci_visualize(lastmodel);



%% --- test/evaluate model on the last session ---

% apply the previously learned model (and also extract the original labels associated with the test
% data for comparisons)
[prediction,loss,teststats,targets] = bci_predict(lastmodel,session{end});

% display the results
disp(['test mis-classification rate: ' num2str(loss*100,3) '%']);
disp(['  predicted classes: ',num2str(round(prediction{2}*prediction{3})')]);  % class probabilities * class values
disp(['  true classes     : ',num2str(round(targets)')]);



%% --- annotate the data set with continuous BCI predictions ---

% append channels that carry the BCI predictions
annotated = bci_annotate(lastmodel,session{end})

% and visualize their time courses together with the rest of the EEG...
% (here, we standardize all channels so that the scales of BCI signals and raw data match up)
pop_eegplot(exp_eval(flt_standardize(annotated)));



%% --- do a pseudo-online simulation at flexible time points ---

% (here, we query the BCI 3.5s after each stimulus marker)
[predictions,latencies] = onl_simulate(session{end},lastmodel,'markers',{'0','16','30'},'offset',3.5)

accuracy = mean(argmax(predictions') == targets');
disp(['online mis-classification rate: ' num2str((1-accuracy)*100,3) '%']);


%% --- do a real-time simulation ---
% ( click into the figure to stop the update (and make sure that your click was registered) )

% play it back in real time
run_readdataset('Dataset',session{end});

% process it in real time using lastmodel, and visualize outputs
run_writevisualization('Model',lastmodel, 'VisFunction','bar(y)');

% make sure that the online processing gets terminated...
disp('Click into the figure to stop online processing.'); 
waitforbuttonpress; onl_clear; close(gcf);


%% --- optionally train an alternative model with parameter search (SLOW!) ---
% note: this takes quite some time: up to 10 minutes
% The search is over possible values for the number of pattern pairs, using Spec-CSP
% (this yields a different number of pattern pairs than the a priori value)

% define approach (here with a parameter search over pattern pairs and a custom spectral prior)
myapproach = {'SpecCSP' 'Prediction',{'FeatureExtraction',{'SpectralPrior',[7 15], 'PatternPairs',search(1,2,3)}}};

% learn model, across 3 classes (0=foot, 16=left, 30=right)
[trainloss,lastmodel,laststats] = bci_train('Data',traindata,'Approach',myapproach,'TargetMarkers',{'0','16','30'});
    
% visualize results
bci_visualize(lastmodel);


%% --- try some more variants of the CSP paradigm (using other signal processing and machine learning methods) ---
% (also test them offline, pseudo-online and in real time)


approaches = [];

% use an FIR filter restricted to the alpha band
approaches.alphafir = {'CSP' 'SignalProcessing',{'FIRFilter',[6 8 14 15]}};
% use an IIR filter instead of the default FIR
approaches.alphaiir = {'CSP' 'SignalProcessing',{'FIRFilter','off', 'IIRFilter',[5 7 14 17]}};
% also restrict the model to a stationary subspace
approaches.stationary = {'CSP' 'SignalProcessing',{'FIRFilter','off', 'IIRFilter',[5 7 14 17], 'StationarySubspace',{'StationaryDim',-0.1,'Operation','keep_stationary'}}};
% use a sharp FFT band-pass filter
approaches.fftspectrum = {'CSP' 'SignalProcessing',{'FIRFilter','off','SpectralSelection',[7 15]}};
% use the basic version (full spectrum from now on)
approaches.basic = 'CSP';
% use a simple logistic regression classifier (variational Bayes) instead of the LDA
approaches.logreg = {'CSP' 'Prediction',{'MachineLearning',{'Learner','logreg'}}};
% use a simple logistic regression classifier (sparse variational Bayes)
approaches.sparselogreg = {'CSP' 'Prediction',{'MachineLearning',{'Learner',{'logreg','variant','vb-ard'}}}};
% use the sparse logistic regression classifier but applied to a larger set of patterns
approaches.bigsparselogreg = {'CSP' 'Prediction',{'MachineLearning',{'Learner',{'logreg','variant','vb-ard'}},'FeatureExtraction',{'PatternPairs',6}}};
% using quadratic discriminant analysis
approaches.qda = {'CSP' 'Prediction',{'MachineLearning',{'Learner','qda'}}};
% using Gaussian mixture models (variational Bayesian Dirichlet process prior)
approaches.vbgmm = {'CSP' 'Prediction',{'MachineLearning',{'Learner','gmm'}}};
% using relevance vector machines (here with a fixed kernel scale for speed)
approaches.rvm = {'CSP' 'Prediction',{'MachineLearning',{'Learner',{'rvm','kernel','rbf','gamma',1}}}};
% use l1-regularized logreg (which involves a parameter search over the regularization parameter); takes 2-3 minutes
approaches.l1logreg = {'CSP' 'Prediction',{'MachineLearning',{'Learner',{'logreg','variant','l1','lambda',search(2.^(-6:1.5:10))}},'FeatureExtraction',{'PatternPairs',6}}};
% use support vector machines (using the SVMlight package); note: also requires a parameter search (takes about 2 minutes)
% some people are getting a segfault for this method on Windows; probably not BCILAB's fault...
if ~ispc
    approaches.svmlight = {'CSP' 'Prediction',{'MachineLearning',{'Learner',{'svmlight','cost',search(2.^(-5:2.5:15)),'gamma',1}}}}; end
% optimize the size of the stationary subspace to retain; note: cross-validation takes place on continuous data (takes approx 2.5 minutes)
approaches.optspace = {'CSP' 'SignalProcessing',{'FIRFilter','off', 'IIRFilter',[5 7 14 17], 'StationarySubspace',{'StationaryDim',search(-0.5:0.1:-0.1),'Operation','keep_stationary'}}};
% optimize the location of the frequency band manually (note: 2-dimensional parameter space -- takes approx. 7 minutes)
approaches.optflt = {'CSP' 'SignalProcessing',{'FIRFilter','off', 'SpectralSelection',[search(7:10) search(14:2:30)]}};
% use a hierarchical kernel learning classifier (note: this takes ~55 minutes and ~13GB ram!)
% (this is 3x as fast offline & online for only 2 classes)
approaches.hkl = {'CSP' 'Prediction',{'MachineLearning',{'Learner','hkl'}}};


% for each of the above approaches...
for app = fieldnames(approaches)'
    fprintf(['\n==== now testing "' app{1} '" ====\n']);
    fprintf([utl_printapproach(approaches.(app{1})) '\n\n']);
    % train & cross-validate
    [trainloss,lastmodel,laststats] = bci_train('Data',traindata,'Approach',approaches.(app{1}),'TargetMarkers',{'0','16','30'})    
    disp(['training mis-classification rate: ' num2str(trainloss*100,3) '%']);
    % test offline
    [prediction,loss,teststats,targets] = bci_predict(lastmodel,session{end});
    disp(['offline test mis-classification rate: ' num2str(loss*100,3) '%']);
    % test pseudo-online
    [predictions,latencies] = onl_simulate(session{end},lastmodel,'markers',{'0','16','30'},'offset',3.5);
    disp(['pseudo-online test mis-classification rate: ' num2str(mean(argmax(predictions') ~= targets')*100,3) '%']);
    % visualize in real time
    run_readdataset('Dataset',session{end}); run_writevisualization('Model',lastmodel, 'VisFunction','bar(y)'); 
    waitforbuttonpress; onl_clear; close(gcf);
end


%% --- try a few other paradigms that are applicable to oscillatory processes ---

approaches = [];

% using the Bandpower paradigm, extracting simple per-channel logarithmic band-power estimates
approaches.logbp = {'Bandpower', 'SignalProcessing',{'FIRFilter',[6 8 28 32], 'SurfaceLaplacian','off'}};
% Bandpower with the surface Laplacian as basic spatial filter
approaches.surflap = {'Bandpower', 'SignalProcessing',{'FIRFilter',[6 8 28 32], 'SurfaceLaplacian','on'}};
% Bandpower with stationary subspaces as spatial filters
approaches.stationary = {'Bandpower', 'SignalProcessing',{'FIRFilter',[6 8 28 32], 'SurfaceLaplacian','off', 'StationarySubspace',{'Operation','separate'}}};

% using the SpectralMeans paradigm for the same frequency window, on raw channels
% Note: if you are running out of memory under Win32, you can skip the specm approaches and continue
% with specpca (you may also run env_clear_memcaches in the command line to free up space).
pproaches.specm_chans = {'Spectralmeans', 'Prediction',{'FeatureExtraction',{'FreqWindows',[7 30]}}};
% Spectralmeans using the surface Laplacian as spatial filters
approaches.specm_surflap = {'Spectralmeans', 'SignalProcessing',{'SurfaceLaplacian','on'},'Prediction',{'FeatureExtraction',{'FreqWindows',[7 30]}}};
% Spectralmeans using stationary subspaces as spatial filters
approaches.specm_stat = {'Spectralmeans', 'SignalProcessing',{'StationarySubspace',{'Operation','separate'}}, 'Prediction',{'FeatureExtraction',{'FreqWindows',[7 30]}}};
% Spectralmeans, but now using more frequency windows per component (7-15, 15-25, 7-30)
approaches.specm_stat_big = {'Spectralmeans', 'SignalProcessing',{'StationarySubspace',{'Operation','separate'}}, 'Prediction',{'FeatureExtraction',{'FreqWindows',[7 15; 15 25; 7 30]}}};
% Spectralmeans, now also using sparse logistic regression as classifier
approaches.specm_stat_big_sparse_logreg = {'Spectralmeans', 'SignalProcessing',{'StationarySubspace',{'Operation','separate'}},  ...
    'Prediction',{'FeatureExtraction',{'FreqWindows',[7 15;15 25;7 30]},'MachineLearning',{'Learner',{'logreg','variant','vb-ard'}}}};

% using the basic DataflowSimplified paradigm, and operating on raw PCA features of the spectra of stationary components (using sparse logistic regression)
approaches.specpca = {'DataflowSimplified', 'SignalProcessing',{'IIRFilter',{[0.1 2],'highpass'}, 'EpochExtraction',[0.5 3.5], 'SpectralTransform',{'multitaper',true,false,80}, 'Resampling',100, 'StationarySubspace',{'Operation','separate'},'EpochPCA',10},  ...
    'Prediction',{'MachineLearning',{'Learner',{'logreg','variant','vb-ard'}}}}; 

% using the vanilla FBCSP paradigm for three frequency windows (filter-bank CSP)
approaches.fbcsp = {'FBCSP' 'Prediction',{'FeatureExtraction',{'FreqWindows',[7 15; 15 25; 7 30]}}};
% using FBCSP with logistic regression
approaches.fbcsp_logreg = {'FBCSP' 'Prediction',{'FeatureExtraction',{'FreqWindows',[7 15; 15 25; 7 30]},'MachineLearning',{'Learner','logreg'}}};
% using FBCSP with a somewhat elaborate time windowing approach
approaches.fbcsp_hann = {'FBCSP' 'SignalProcessing',{'EpochExtraction',[0 4]},'Prediction',{'FeatureExtraction',{'FreqWindows',[7 15; 15 25; 7 30],'WindowFunction','hann'}}};

% using RCSP for Tikhonov-regularized CSP (takes approx. 2 minutes)
approaches.trcsp = {'RCSP' 'Prediction',{'FeatureExtraction',{'beta',0,'gamma',0}}};
% using RCSP for diagonal-loading CSP (takes approx. 2 minutes)
approaches.DLCSPcv = {'RCSP' 'Prediction',{'FeatureExtraction',{'beta',0,'alpha',0}}};
% using RCSP for diagonal-loading CSP (with analytically-derived shrinkage, fast)
approaches.DLCSPauto = {'RCSP' 'Prediction',{'FeatureExtraction',{'beta',0,'alpha',0,'gamma','auto'}}};
% automatically choose a good classifier to use with these features (takes approx. 5-11 minutes)
approaches.DLCSPauto_classifier = {'RCSP' 'Prediction',{'FeatureExtraction',{'beta',0,'alpha',0,'gamma','auto'},'MachineLearning',{'Learner',search('lda','qda','gmm','logreg')}}};

% using the DALOSC paradigm (oscillatory rank-regularized logistic regression on second-order signal features), restricted to a somewhat smaller parameter search range than default
% (as it's otherwise too slow in this multi-class constellation; still takes approx. 5-10 minutes)
approaches.dalosc = {'DALOSC' 'Prediction',{'MachineLearning',{'Learner',{'dal','lambdas',2.^(10:-0.33:-3)}}}};

% using the general DAL paradigm (which supports multiple windows) in its vanilla configuration, but with a restricted search interval
approaches.dalbasic = {'DAL' 'Prediction',{'MachineLearning',{'Learner',{'dal','lambdas',2.^(10:-0.33:-3)}}}};

% Spectralmeans, using LARS (sparse logistic regression) with elastic net regularizer as classifier on coherence, phase and power-spectral density features between stationary components (for selected frequency bands); takes ca. 2 minutes
% (requires ~5GB RAM
approaches.specm_cohen = {'Spectralmeans', 'SignalProcessing',{'IIRFilter',{[0.1 2],'highpass'}, 'StationarySubspace',{'StationaryDim',16,'Operation','keep_stationary'},'SpectralTransform','off','CoherenceTransform','on'},  ...
    'Prediction',{'FeatureExtraction',{'FreqWindows',[4 7; 8 15; 15 25; 7 30]},'MachineLearning',{'Learner',{'logreg','variant',{'lars','ElasticMixing',0.5}}}}};


% for each of the above approaches...
for app = fieldnames(approaches)'
    fprintf(['\n==== now testing "' app{1} '" ====\n']);
    try fprintf([utl_printapproach(approaches.(app{1})) '\n\n']); catch,end
    % train & cross-validate
    [trainloss,lastmodel,laststats] = bci_train('Data',traindata,'Approach',approaches.(app{1}),'TargetMarkers',{'0','16','30'})    
    disp(['training mis-classification rate: ' num2str(trainloss*100,3) '%']);
    % test offline
    [prediction,loss,teststats,targets] = bci_predict(lastmodel,session{end});
    disp(['offline test mis-classification rate: ' num2str(loss*100,3) '%']);
    % test pseudo-online
    [predictions,latencies] = onl_simulate(session{end},lastmodel,'markers',{'0','16','30'},'offset',3.5);
    disp(['pseudo-online test mis-classification rate: ' num2str(mean(argmax(predictions') ~= targets')*100,3) '%']);
    % visualize in real time
    run_readdataset('Dataset',session{end}); run_writevisualization('Model',lastmodel, 'VisFunction','bar(y);ylim([0 1])'); 
    waitforbuttonpress; onl_clear; close(gcf);
end


%% --- try session-wise cross-validation ---

[trainloss,lastmodel,laststats] = bci_train('Data',session,'Approach','Bandpower','TargetMarkers',{'0','16','30'})    
% note: testing on the training data in the following, just to check correct model format
[prediction,loss,teststats,targets] = bci_predict(lastmodel,session{end});
[predictions,latencies] = onl_simulate(session{end},lastmodel,'markers',{'0','16','30'},'offset',3.5);


%% --- test parallel online processing and training ---
% note: this is relatively new code, might display some issues

% first train a model
myapproach1 = {'FBCSP' 'Prediction',{'FeatureExtraction',{'FreqWindows',[7 15; 15 25; 7 30]}}};
[trainloss,lastmodel,laststats] = bci_train('Data',traindata,'Approach',myapproach1,'TargetMarkers',{'0','16','30'});
% now run it online in the background
run_readdataset('MatlabStream','background','Dataset',session{end}); run_writevisualization('Model',lastmodel, 'VisFunction','@(y,f)bar3(get(f,''CurrentAxes''),y)','SourceStream','background','PredictorName','back_pred');

% now train, test, simulate and run another model while the first one is running
myapproach2 = {'Spectralmeans', 'SignalProcessing',{'StationarySubspace',{'Operation','separate'}}, 'Prediction',{'FeatureExtraction',{'FreqWindows',[7 30]}}};
[trainloss,lastmodel,laststats] = bci_train('Data',traindata,'Approach',myapproach2,'TargetMarkers',{'0','16','30'})
disp(['training mis-classification rate: ' num2str(trainloss*100,3) '%']);
% test offline
[prediction,loss,teststats,targets] = bci_predict(lastmodel,session{end});
disp(['offline test mis-classification rate: ' num2str(loss*100,3) '%']);
% test pseudo-online
[predictions,latencies] = onl_simulate(session{end},lastmodel,'markers',{'0','16','30'},'offset',3.5);
disp(['pseudo-online test mis-classification rate: ' num2str(mean(argmax(predictions') ~= targets')*100,3) '%']);
% visualize in real time
run_readdataset('Dataset',session{end}); run_writevisualization('Model',lastmodel, 'VisFunction','@(y,f)bar3(get(f,''CurrentAxes''),y)'); waitforbuttonpress; onl_clear;

