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
    filename = ['bcilab:/userdata/tutorial/imag_movements2/session' num2str(s) '.bdf'];
    % load & retain the first 32 channels
    session{s} = exp_eval_optimized(io_loadset(filename,'channels',1:32));
    % override the channel locations (they are incorrect in the original data)
    session{s}.chanlocs = set_infer_chanlocs('bcilab:/userdata/tutorial/imag_movements2/mi32.loc');
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
run_writevisualization('Model',lastmodel);

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
% For each defined approach we test it offline, pseudo-online and in real time.
% Note that most of these variants are for illustrative purposes only and will actually reduce the
% performance of the default CSP method (as the default settings are pretty good).

 approaches = [];
 
% --- here we apply a few alternative temporal (and in particular spectral) filters ---

% use an FIR filter restricted to the alpha band
approaches.alphafir = {'CSP' 'SignalProcessing',{'FIRFilter',[6 8 14 15]}};
% use an IIR filter instead of the default FIR
approaches.alphaiir = {'CSP' 'SignalProcessing',{'FIRFilter','off', 'IIRFilter',[5 7 14 17]}};
% use an IIR filter with free-form response (peaks both in alpha and beta band)
approaches.alphabetaiir = {'CSP' 'SignalProcessing',{'FIRFilter','off', 'IIRFilter',{[5 7 14 17 20 25 30;0 1 1 0 0.6 0.6 0],'freeform','yule-walker'}}};
% use a sharp FFT-based band-pass filter, applied to each epoch
approaches.fftspectrum = {'CSP' 'SignalProcessing',{'FIRFilter','off','SpectralSelection',[7 15]}};
% use the default frequency selection
approaches.defaults = {'CSP'};
% apply a window in the time domain (covers whole epoch)
approaches.windowed = {'CSP' 'SignalProcessing',{'WindowSelection','hann'}};
% apply to data resampled to 70 Hz
approaches.resampled = {'CSP' 'SignalProcessing',{'Resampling',70}};
% apply a delay-embedding, that is, append channels that contain lagged (shifted) versions
approaches.delayembed = {'CSP' 'SignalProcessing',{'DelayEmbedding',3}};
% try again, but now search for the optimal lag, and exclude intermediate lags (same as CSSP, but here with limited search range); takes ca. 2m
approaches.delayembed_search = {'CSP' 'SignalProcessing',{'DelayEmbedding',{'NumLags',search(1:5),'IncludeIntermediates',false}}};
% apply the delay-embedding to resampled data (note: order of filters in SignalProcessing does not matter, is automatically reordered)
approaches.delayembed_res = {'CSP' 'SignalProcessing',{'DelayEmbedding',3,'Resampling',70}};
% standardize the channels in a moving window, using defaults (note: this is more useful when going across sessions or subjects)
approaches.standardize = {'CSP' 'SignalProcessing',{'Standardization','on'}};
% sphere/whiten the data in a 20s moving window
approaches.sphere = {'CSP' 'SignalProcessing',{'Standardization',{'Sphere',true,'WindowLength',20}}};

% --- here we apply some spatial filters ---
% (note: since CSP is already an optimal spatial filter there is little gain by inserting extra stages, except for additional assumptions/constraints)

% apply a random projection matrix (note: this should give the same result as basic CSP since CSP
% finds an optimal linear transformation regardless of how the data has been linearly transformed
% beforehand, as long as no information is lost)
approaches.project = {'CSP' 'SignalProcessing',{'Projection',randn(32)}};
% re-reference to common average reference (can be slightly worse than original as the signal's rank is reduced by 1)
approaches.reref = {'CSP' 'SignalProcessing',{'Rereferencing','on'}};
% re-reference to common median reference (robust reference, not strictly a rank reduction)
approaches.reref_med = {'CSP' 'SignalProcessing',{'Rereferencing',{'ReferenceType','median'}}};
% surface laplacian
approaches.laplace = {'CSP' 'SignalProcessing',{'SurfaceLaplacian','on'}};
% select a subset of channels (note: the channels in the data are actually not in the 10-20 system, so we have to look up the closest channels in the data)
approaches.subset = {'CSP' 'SignalProcessing',{'ChannelSelection',{'Channels',{'C3','Cz','C4'},'FindClosest',true}}};
% perform both a surface laplacian, and channel selection, and enforce a desired order on the two filters (here we want to select channels *after* the surface laplacian has been applied)
approaches.laplace_subset = {'CSP' 'SignalProcessing',{'FilterOrdering',{'SurfaceLaplacian','ChannelSelection'}, 'SurfaceLaplacian','on', 'ChannelSelection',{'Channels',{'C3','Cz','C4'},'FindClosest',true}}};
% restrict the model to a stationary subspace
approaches.stationary = {'CSP' 'SignalProcessing',{'StationarySubspace',{'StationaryDim',-0.1,'Operation','keep_stationary'}}};

% --- here we apply some data cleaning filters ---

% apply all-inclusive data cleaning with current default settings (note: if you're classifying partially based on artifacts, 
% or if the thresholds are too tight then your performance may suffer)
approaches.defaultclean = {'CSP' 'SignalProcessing',{'DataCleaning','on'}}; % to use a particular BCILAB version's defaults: {'CSP' 'SignalProcessing',{'DataCleaning','1.1-beta'}};
% apply data cleaning and override some parameters (here: disable removal of noisy time windows); note the double cell array since we're too lazy to explicitly assign to the (first) argument of the DataCleaning named 'DataSetting'
approaches.defaultclean_nownd = {'CSP' 'SignalProcessing',{'DataCleaning',{{'1.1-beta','BadWindowRemoval','off'}}}};
% start with the 'highpass' setting and enable just bad window removal, setting the thresholds to -5 and +7 standard deviations
approaches.justwindows57 = {'CSP' 'SignalProcessing',{'DataCleaning',{{'highpass','BadWindowRemoval',[-5 7]}}}};
% start with the 'highpass' setting and enable just bad subspace removal, setting the threshold to +7 standard deviations
approaches.justsubspace = {'CSP' 'SignalProcessing',{'DataCleaning',{{'highpass','BadSubspaceRemoval',7}}}};

% --- here we apply a few of the different classifiers ---

% use the basic version (full spectrum from now on), using an LDA (Linear Discriminant Analysis) classifier
approaches.lda = 'CSP';
% use LDA but account for unbalanced classes in training data (assume that the same class balance will hold on test data)
% note: on these data the number of trials for both classes are approx. the same, so there's little difference
approaches.lda_unbalanced = {'CSP' 'Prediction',{'MachineLearning',{'Learner',{'lda','WeightedBias',true}}}};
% use a simple logistic regression classifier (variational Bayes) instead of the LDA
approaches.logreg = {'CSP' 'Prediction',{'MachineLearning',{'Learner','logreg'}}};
% logreg again, but this time using the 1-vs-rest voting scheme to handle 3 classes using binary classifiers (instead of the default 1-vs-1)
approaches.logreg_1vR = {'CSP' 'Prediction',{'MachineLearning',{'Learner',{'logreg','VotingScheme','1vR'}}}};
% use a simple logistic regression classifier (sparse variational Bayes)
approaches.sparselogreg = {'CSP' 'Prediction',{'MachineLearning',{'Learner',{'logreg','variant','vb-ard'}}}};
% use the sparse logistic regression classifier but applied to a larger set of patterns
approaches.bigsparselogreg = {'CSP' 'Prediction',{'MachineLearning',{'Learner',{'logreg','Variant','vb-ard'}},'FeatureExtraction',{'PatternPairs',6}}};
% use another sparse logistic regression classifier (not variational Bayes but using the very fast LARS algorithm)
approaches.bigsparselogreg_lars = {'CSP' 'Prediction',{'MachineLearning',{'Learner',{'logreg','Variant','lars'}},'FeatureExtraction',{'PatternPairs',6}}};
% implement sparse logistic regression using the Bayesian generalized linear models inference and estimation toolbox (glm-ie) and find the overall sparsity level using evidence maximization
approaches.bigsparselogreg_glm = {'CSP' 'Prediction',{'MachineLearning',{'Learner',{'glm','Type','classification','Priors',{'Term1',{'Laplace','Scales',2.^(-4:1:4)}}}},'FeatureExtraction',{'PatternPairs',6}}};
% using quadratic discriminant analysis
approaches.qda = {'CSP' 'Prediction',{'MachineLearning',{'Learner','qda'}}};
% using Gaussian mixture models (variational Bayesian Dirichlet process prior)
approaches.vbgmm = {'CSP' 'Prediction',{'MachineLearning',{'Learner','gmm'}}};
% using a different algorithm for Gaussian mixture models (fixed number of clusters per class, using the EM algorithm)
approaches.emgmm = {'CSP' 'Prediction',{'MachineLearning',{'Learner',{'gmm','Variant','Greedy Expectation-Maximization','NumClusters',2}}}};
% using relevance vector machines (here with a slightly coarsened range for the kernel scale to increase processing speed); 
% note that the RVM by default searches the optimal gamma parameter internally using evidence maximization (aka Empirical Bayes), so no need for a search()
approaches.rvm = {'CSP' 'Prediction',{'MachineLearning',{'Learner',{'rvm','kernel','rbf','gamma',2.^(-16:0.5:10)}}}};

% --- the following methods demonstrate parameter search (these are usually slower than the above) ---

% use l1-regularized logreg (which involves a parameter search over the regularization parameter); takes 2-3 minutes
approaches.l1logreg = {'CSP' 'Prediction',{'MachineLearning',{'Learner',{'logreg','variant','l1','lambda',search(2.^(-6:1.5:10))}},'FeatureExtraction',{'PatternPairs',6}}};
% use a basic linear SVM
approaches.svmlin = {'CSP' 'Prediction',{'MachineLearning',{'Learner',{'svm','kernel','linear','gamma',1}}}};
% use support vector machines (using the SVMlight package); note: also requires a parameter search (takes about 2 minutes)
% some people are getting a segfault for this method on Windows; probably not BCILAB's fault...
if ~ispc
    approaches.svmlight = {'CSP' 'Prediction',{'MachineLearning',{'Learner',{'svmlight','cost',search(2.^(-5:2.5:15)),'gamma',1}}}}; end
% optimize the size of the stationary subspace to retain; note: cross-validation takes place on continuous data (takes approx 2.5 minutes)
approaches.optspace = {'CSP' 'SignalProcessing',{'FIRFilter','off', 'IIRFilter',[5 7 14 17], 'StationarySubspace',{'StationaryDim',search(-0.5:0.1:-0.1),'Operation','keep_stationary'}}};
% optimize the location of the frequency band manually (note: 2-dimensional parameter space -- takes approx. 7 minutes)
approaches.optflt = {'CSP' 'SignalProcessing',{'FIRFilter','off', 'SpectralSelection',[search(7:10) search(13:2:30)]}};
% use a hierarchical kernel learning classifier (note: this takes ~25 minutes and ~13GB ram!)
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
    run_readdataset('Dataset',session{end}); run_writevisualization('Model',lastmodel); 
    waitforbuttonpress; clear laststream; close(gcf);
end


%% --- try a few other paradigms that are applicable to oscillatory processes ---

approaches = [];

% --- the Bandpower paradigm takes the log-variance of each channel (or component) in the whole epoch ---

% using the Bandpower paradigm, extracting simple per-channel logarithmic band-power estimates
approaches.logbp = {'Bandpower', 'SignalProcessing',{'FIRFilter',[6 8 28 32], 'SurfaceLaplacian','off'}};
% Bandpower with the surface Laplacian as basic spatial filter
approaches.surflap = {'Bandpower', 'SignalProcessing',{'FIRFilter',[6 8 28 32], 'SurfaceLaplacian','on'}};
% Bandpower with stationary subspaces as spatial filters
approaches.stationary = {'Bandpower', 'SignalProcessing',{'FIRFilter',[6 8 28 32], 'SurfaceLaplacian','off', 'StationarySubspace',{'Operation','separate'}}};
% Bandpower of decorrelated ("sphered") channels
approaches.sphered = {'Bandpower', 'SignalProcessing',{'FIRFilter',[6 8 28 32], 'SurfaceLaplacian','off', 'ICA',{'Variant','sphere','DataCleaning','off'}}};
% Bandpower of random spatial projections of the signal, using logistic regression
approaches.randproj = {'Bandpower', 'SignalProcessing',{'FIRFilter',[6 8 28 32], 'Projection',randn(256,32)}, 'Prediction',{'MachineLearning',{'Learner','logreg'}}};

% --- SpectralMeans allows for specifying one or more spectral windows in which to take the power ---

% using the SpectralMeans paradigm for the same frequency window, on raw channels
% Note: if you are running out of memory under Win32, you can skip the specm approaches and continue
% with specpca (you may also run env_clear_memcaches in the command line to free up space).
approaches.specm_chans = {'Spectralmeans', 'Prediction',{'FeatureExtraction',{'FreqWindows',[7 30]}}};
% same but with log-transform enabled (more likely to be linearly separable)
approaches.specm_chans_log = {'Spectralmeans', 'Prediction',{'FeatureExtraction',{'FreqWindows',[7 30],'LogTransform',true}}};
% on channels, but with multiple frequencies (too many features for LDA)
approaches.specm_chans_multi = {'Spectralmeans', 'Prediction',{'FeatureExtraction',{'FreqWindows',[7 15; 15 25; 7 30]}}};
% multiple frequencies using a sparse classifier (no problem with too many features)
approaches.specm_chans_multi_sparse = {'Spectralmeans', 'Prediction',{'FeatureExtraction',{'FreqWindows',[7 15; 15 25; 7 30]},'MachineLearning',{'Learner',{'logreg','variant','vb-ard'}}}};
% Spectralmeans using the surface Laplacian as spatial filters
approaches.specm_surflap = {'Spectralmeans', 'SignalProcessing',{'SurfaceLaplacian','on'},'Prediction',{'FeatureExtraction',{'FreqWindows',[7 30]}}};

% --- the Filter-Bank CSP paradigm applies CSP to multiple selectable frequency bands ---

% using the vanilla FBCSP paradigm for three frequency windows (filter-bank CSP)
approaches.fbcsp = {'FBCSP' 'Prediction',{'FeatureExtraction',{'FreqWindows',[7 15; 15 25; 7 30]}}};
% using FBCSP with a somewhat fancier time windowing approach
approaches.fbcsp_hann = {'FBCSP' 'SignalProcessing',{'EpochExtraction',[0 4]},'Prediction',{'FeatureExtraction',{'FreqWindows',[7 15; 15 25; 7 30],'WindowFunction','hann'}}};
% using FBCSP with logistic regression
approaches.fbcsp_logreg = {'FBCSP' 'Prediction',{'FeatureExtraction',{'FreqWindows',[7 15; 15 25; 7 30]},'MachineLearning',{'Learner','logreg'}}};
% ... using robust covariance estimation
approaches.fbcsp_robust = {'FBCSP' 'Prediction',{'FeatureExtraction',{'FreqWindows',[7 15; 15 25; 7 30],'RobustCovariance',true},'MachineLearning',{'Learner','logreg'}}};
% ... using shrinkage covariance estimation (takes ca. 1m)
approaches.fbcsp_shrink = {'FBCSP' 'Prediction',{'FeatureExtraction',{'FreqWindows',[7 15; 15 25; 7 30],'ShrinkageCovariance',true},'MachineLearning',{'Learner','logreg'}}};

% --- the Regularized CSP paradigm is a regularized version of CSP (works with less data) ---

% using RCSP for Tikhonov-regularized CSP (takes approx. 2 minutes)
approaches.trcsp = {'RCSP' 'Prediction',{'FeatureExtraction',{'beta',0,'gamma',0}}};
% using RCSP for diagonal-loading CSP (takes approx. 2 minutes)
approaches.DLCSPcv = {'RCSP' 'Prediction',{'FeatureExtraction',{'beta',0,'alpha',0}}};
% using RCSP for diagonal-loading CSP (with analytically-derived shrinkage, fast)
approaches.DLCSPauto = {'RCSP' 'Prediction',{'FeatureExtraction',{'beta',0,'alpha',0,'gamma','auto'}}};
% automatically choose a good classifier to use with these features (takes approx. 5-11 minutes)
approaches.DLCSPauto_classifier = {'RCSP' 'Prediction',{'FeatureExtraction',{'beta',0,'alpha',0,'gamma','auto'},'MachineLearning',{'Learner',search('lda','qda','gmm','logreg')}}};

% --- the DAL method is a rank-regularized logistic regression applied to the covariance matrix of each epoch ---
% ... where it can implicitly perform the action of a spatial filter (Dual-Augmented Lagrangian method for Oscillatory Processes)

% using the DALOSC paradigm (oscillatory rank-regularized logistic regression on second-order signal features), restricted to a somewhat smaller parameter search range than default
% (as it's otherwise too slow in this multi-class constellation; still takes approx. 5-10 minutes)
approaches.dalosc = {'DALOSC' 'Prediction',{'MachineLearning',{'Learner',{'dal','lambdas',2.^(10:-0.33:-3)}}}};
% using the generic DAL paradigm is an extension to multiple frequency bands, and in the low-frequency case operates on event-related potentials (like DALERP, see also tutorial_erp1)
approaches.dalgeneric = {'DAL' 'Prediction',{'FeatureExtraction',{'WindowFreqs',[0.5 5; 7 15; 15 25; 7 30]},'MachineLearning',{'Learner',{'dal','lambdas',2.^(10:-0.33:-3)}}}};

% --- the DataflowSimplified paradigm passes the pre-processed epochs directly to the machine learning stage ---
% (normally this is the foundation on which other more specific paradigms are implemented)

% using the basic DataflowSimplified paradigm, and operating on raw PCA features of the spectra of stationary components (using sparse logistic regression)
approaches.specpca = {'DataflowSimplified', 'SignalProcessing',{'IIRFilter',{[0.1 2],'highpass'}, 'EpochExtraction',[0.5 3.5], 'SpectralTransform',{'multitaper',true,false,80}, 'Resampling',100, 'StationarySubspace',{'Operation','separate'},'EpochPCA',10},  ...
    'Prediction',{'MachineLearning',{'Learner',{'logreg','variant','vb-ard'}}}}; 
% operate on the analytic amplitudes of three characteristic frequencies over time using simple l2-regularized logistic regression
approaches.ampphase = {'DataflowSimplified', 'SignalProcessing',{'Resampling',70,'AnalyticPhasor',{{'hilbert',{[3 4 6 7],[7 8 12 13],[18 20 22 25]}}}, 'EpochExtraction',[0.5 3.5]},  ...
    'Prediction',{'MachineLearning',{'Learner',{'logreg','variant','l2'}}}};

% --- miscellaneous approaches ---

% Spectralmeans, using LARS (sparse logistic regression) with elastic net regularizer as classifier on coherence, phase and power-spectral density features between stationary components (for selected frequency bands); takes ca. 2 minutes
% (requires ~5GB RAM
approaches.specm_cohen = {'Spectralmeans', 'SignalProcessing',{'IIRFilter',{[0.1 2],'highpass'}, 'StationarySubspace',{'StationaryDim',16,'Operation','keep_stationary'},'SpectralTransform','off','CoherenceTransform','on'},  ...
    'Prediction',{'FeatureExtraction',{'FreqWindows',[4 7; 8 15; 15 25; 7 30]},'MachineLearning',{'Learner',{'logreg','variant',{'lars','ElasticMixing',0.5}}}}};

% --- examples of very high-dimensional feature spaces (here: time/frequency representations for each channel, using various classifiers) ---
% WARNING: these are quite slow (see time estimates)

% operate on the event-related time/frequency representation (ERSP) of each channel using fast l2-regularized logistic regression
approaches.ersp_l2 = {'DataflowSimplified', 'SignalProcessing',{'Resampling',100, 'ERSPTransform',{'SpectralMap','log'}, 'EpochExtraction',[0.5 3.5]}, ...
    'Prediction',{'FeatureExtraction','vectorized','MachineLearning',{'Learner',{'logreg','variant','l2'}}}};
% using fast l1-regularized logistic regression (no search ver reg. parameter)
approaches.ersp_l1 = {'DataflowSimplified', 'SignalProcessing',{'Resampling',100, 'ERSPTransform',{'SpectralMap','log'}, 'EpochExtraction',[0.5 3.5]}, ...
    'Prediction',{'FeatureExtraction','vectorized','MachineLearning',{'Learner',{'logreg','variant','l1'}}}};
% using elastic net regularized logistic regression (with fast parameter search); takes ca. 10m
approaches.ersp_lars = {'DataflowSimplified', 'SignalProcessing',{'Resampling',100, 'ERSPTransform',{'SpectralMap','log'}, 'EpochExtraction',[0.5 3.5]}, ...
    'Prediction',{'FeatureExtraction','vectorized','MachineLearning',{'Learner',{'logreg','variant',{'lars','ElasticMixing',0.5}}}}};
% using logistic regression with gaussian prior and empirically estimated hyper-parameter; takes ca. 30m
approaches.ersp_glm_gauss = {'DataflowSimplified', 'SignalProcessing',{'Resampling',100, 'ERSPTransform',{'SpectralMap','log'}, 'EpochExtraction',[0.5 3.5]}, ...
    'Prediction',{'FeatureExtraction','vectorized','MachineLearning',{'Learner',{'glm','Type','classification','Lambdas',2.^(-4:1:4),'Priors',{'Term1','Gaussian'}}}}};
% using logistic regression with elastic net prior and empirically estimated hyper-parameters
approaches.ersp_glm_enet = {'DataflowSimplified', 'SignalProcessing',{'Resampling',100, 'ERSPTransform',{'SpectralMap','log'}, 'EpochExtraction',[0.5 3.5]}, ...
    'Prediction',{'FeatureExtraction','vectorized','MachineLearning',{'Learner',{'glm','Type','classification','Lambdas',2.^(-4:1:4),'Priors',{'Term1','Gaussian','Term2',{'Laplace','Scales',2.^(-4:1:4)}}}}}};
% using logistic regression with combined markov random field and sparse priors and empirically estimated hyper-parameters
approaches.ersp_glm_mrf_lap = {'DataflowSimplified', 'SignalProcessing',{'Resampling',100, 'ERSPTransform',{'SpectralMap','log'}, 'EpochExtraction',[0.5 3.5]}, ...
    'Prediction',{'FeatureExtraction','vectorized','MachineLearning',{'Learner',{'glm','Type','classification','Lambdas',2.^(-4:1:4),'Priors',{'Term1',{'Gaussian','LinearOperator','@(x)[vec(diff(x,2));vec(diff(x,3))]'},'Term2',{'Laplace','Scales',2.^(-4:1:4)}}}}}};
% using logistic regression with combined markov random field and sparse priors and empirically estimated hyper-parameters, with separate smoothness across time and frequency
approaches.ersp_glm_mrf_lap = {'DataflowSimplified', 'SignalProcessing',{'Resampling',100, 'ERSPTransform',{'SpectralMap','log'}, 'EpochExtraction',[0.5 3.5]}, ...
    'Prediction',{'FeatureExtraction','vectorized','MachineLearning',{'Learner',{'glm','Type','classification','Priors',{ ...
    'Term1','Gaussian', ...
    'Term2',{'Gaussian','AppendToX',false,'LinearOperator','@(x)vec(diff(x,2))','Scales',2.^(-4:1:4)}, ...
    'Term3',{'Gaussian','AppendToX',false,'LinearOperator','@(x)vec(diff(x,3))','Scales',2.^(-4:1:4)}, ...
    'Term4',{'Laplace','Scales',2.^(-4:1:4)}}}}}};


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
    run_readdataset('Dataset',session{end}); run_writevisualization('Model',lastmodel); 
    try waitforbuttonpress; catch,end; clear laststream; close(gcf);
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

