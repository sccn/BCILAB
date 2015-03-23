%% --- tutorial script for ERP analysis ---

% This tutorial demonstrates the use of BCILAB to learn basic predictive models using 
% ERP-like phenomena, here operating on the error negativity (in a Flanker task).
%
% In this data set, there are two different markers for the correct response ('S101' and 'S102',
% the first for left-hand responses, and the second for right-hand responses), as well as and 
% two different markers for incorrect responses ('S201','S202', again left/right hand).
%
% A second data set, used for testing here, was recorded from an identical twin.
%
% The data is courtesy of Grainne McLoughlin, King's College.
%
% Reference:
%   Study "Neurophysiology of Attention and Activity in Twins" (Institute of Psychiatry, 
%   King's College London). Principal Investigator: Gráinne McLoughlin. 
%
%#ok<*ASGLU,*NASGU,*SNASGU> % turn off a few editor warnings...



%% --- using the WindowMeans method ---
% define markers; here, two groups of markers are being defined; the first group represents class 1
% (correct responses), and the second group represents class 2 (incorrect responses).
mrks = {{'S101','S102'},{'S201','S202'}};

% define ERP windows of interest; here, 7 consecutive windows of 50ms length each are being
% specified, starting from 250ms after the subject response
wnds = [0.25 0.3;0.3 0.35;0.35 0.4; 0.4 0.45;0.45 0.5;0.5 0.55;0.55 0.6];

% define load training data (BrainVision format)
traindata = io_loadset('bcilab:/userdata/tutorial/flanker_task/12-08-002_ERN.vhdr');

% define approach
myapproach = {'Windowmeans' 'SignalProcessing', {'Resampling','off','EpochExtraction',[-0.2 0.8],'SpectralSelection',[0.1 15]}, 'Prediction',{'FeatureExtraction',{'TimeWindows',wnds}}};

%learn model 
[trainloss,lastmodel,laststats] = bci_train('Data',traindata,'Approach',myapproach,'TargetMarkers',mrks);
disp(['training mis-classification rate: ' num2str(trainloss*100,3) '%']);

% visualize results
bci_visualize(lastmodel)


%% --- test on the twin's data set ---

% define test data
testdata = io_loadset('bcilab:/userdata/tutorial/flanker_task/12-08-001_ERN.vhdr');
[prediction,loss,teststats,targets] = bci_predict(lastmodel,testdata);

% result visualization
disp(['test mis-classification rate: ' num2str(loss*100,3) '%']);
disp(['  predicted classes: ',num2str(round(prediction{2}*prediction{3})')]);  % class probabilities * class values
disp(['  true classes     : ',num2str(round(targets)')]);


%% --- do a pseudo-online simulation ---

% simulate online processing with 10 updates per second, and epoch extraction relative to the same target markers
prediction2 = onl_simulate('Data',testdata,'Model',lastmodel,'UpdateRate',10,'TargetMarkers',{'S101','S102','S201','S202'});
fprintf('mean difference in predicted classes: %f\n',mean(abs(argmax(prediction{2}')-argmax(prediction2'))));

% also simulate sliding-window online processing, without locking to markers (updating 5x per second for faster analysis)
prediction3 = onl_simulate('Data',testdata,'Model',lastmodel,'UpdateRate',5);

%% --- do a real-time test ---
% (note that you would query the online stream only after distinct events)
% ( click into the figure to stop the update (and make sure that your click was registered) )

% play it back in real time
run_readdataset('Dataset',testdata);

% process it in real time using lastmodel, and visualize outputs
run_writevisualization('Model',lastmodel, 'VisFunction','bar(y);ylim([0 1])');

% make sure that the online processing gets terminated...
disp('Click into the figure to stop online processing.'); 
waitforbuttonpress; onl_clear; close(gcf);



%% --- train an alternative model, using a sparse classifier ---

% define markers; here, two groups of markers are being defined; the first group represents class 1
% (correct responses), and the second group represents class 2 (incorrect responses).
mrks = {{'S101','S102'},{'S201','S202'}};

% define ERP windows of interest; here, 7 consecutive windows of 50ms length each are being
% specified, starting from 250ms after the subject response
wnds = [0.25 0.3;0.3 0.35;0.35 0.4; 0.4 0.45;0.45 0.5;0.5 0.55;0.55 0.6];

% define load training data (BrainVision format)
traindata = io_loadset('bcilab:/userdata/tutorial/flanker_task/12-08-002_ERN.vhdr');

% define approach
myapproach = {'Windowmeans' ...
    'SignalProcessing', {'EpochExtraction',[0 0.8], 'SpectralSelection',[0.1 15]}, ...
    'Prediction',{'FeatureExtraction',{'TimeWindows',wnds}, 'MachineLearning',{'Learner', {'logreg', [],'Variant','vb-ard'}}}};

%learn model 
[trainloss,lastmodel,laststats] = bci_train('Data',traindata,'Approach',myapproach,'TargetMarkers',mrks);

% visualize results
bci_visualize(lastmodel)


%% --- try a few alternative ERP approaches ---

wnds = [0.25 0.3;0.3 0.35;0.35 0.4; 0.4 0.45;0.45 0.5;0.5 0.55;0.55 0.6];
traindata = io_loadset('bcilab:/userdata/tutorial/flanker_task/12-08-002_ERN.vhdr');
testdata = io_loadset('bcilab:/userdata/tutorial/flanker_task/12-08-001_ERN.vhdr');

approaches = [];

% the basic approach from earlier
approaches.wmeans = {'Windowmeans' 'SignalProcessing',{'EpochExtraction',[0 0.8],'SpectralSelection',[0.1 15]}, ...
    'Prediction',{'FeatureExtraction',{'TimeWindows',wnds}}};
% using hierarchical discriminant component analysis
approaches.hdca = {'Windowmeans' 'SignalProcessing',{'EpochExtraction',[0 0.8],'SpectralSelection',[0.1 15]}, ...
   'Prediction',{'FeatureExtraction',{'TimeWindows',wnds,'VectorizeFeatures',false}, 'MachineLearning',{'Learner','hdca'}}};
% now with logistic regression instead of LDA
approaches.wmeanslr = {'Windowmeans' 'SignalProcessing',{'EpochExtraction',[0 0.8],'SpectralSelection',[0.1 15]}, ...
   'Prediction',{'FeatureExtraction',{'TimeWindows',wnds}, 'MachineLearning',{'Learner','logreg'}}};
% now with linear Support Vector Machines (takes ca. 2 minutes)
approaches.wmeanssvm = {'Windowmeans' 'SignalProcessing',{'EpochExtraction',[0 0.8],'SpectralSelection',[0.1 15]}, ...
    'Prediction',{'FeatureExtraction',{'TimeWindows',wnds}, 'MachineLearning',{'Learner',{'svm','kernel','linear','gamma',1}}}};

% now with sparse logistic regression (using few weights)
approaches.wmeanslrard = {'Windowmeans' 'SignalProcessing',{'EpochExtraction',[0 0.8],'SpectralSelection',[0.1 15]}, ...
    'Prediction',{'FeatureExtraction',{'TimeWindows',wnds}, 'MachineLearning',{'Learner',{'logreg' 'variant','vb-ard'}}}};
% now with group-sparse logistic regression (using few time windows, takes ca. 3 minutes)
approaches.wmeanslrgl = {'Windowmeans' 'SignalProcessing',{'EpochExtraction',[0 0.8],'SpectralSelection',[0.1 15]}, ...
    'Prediction',{'FeatureExtraction',{'TimeWindows',wnds}, 'MachineLearning',{'Learner',{'dal','scaling','std','regularizer','grouplasso-columns','shape',[64 NaN]}}}};

% using the DALERP paradigm (which operates on the raw signal using logistic regression, but imposes a low-rank assumption), takes ca. 5 minutes
approaches.dalerp = {'DALERP','SignalProcessing',{'Resampling',60,'IIRFilter','off','spectrum',[0.1 15],'EpochExtraction',[0 0.8]}, ...
    'Prediction',{'MachineLearning',{'Learner',{'dal','solver','cg'}}}};

% using the basic dataflow paradigm, but running on a multi-level wavelet transform of the signal (using SVMs)
approaches.waveletssvm = {'DataflowSimplified' 'SignalProcessing',{'EpochExtraction',[0 0.8],'WaveletTransform','on'}, ...
   'Prediction',{'MachineLearning',{'Learner',{'svm','kernel','linear','gamma',1}}}};


% some computationally rather expensive approaches below...

% using LARS (sparse least-angle (logistic) regression) on a wavelet decomposition of independent components of the data
approaches.larsica_wavelets = {'DataflowSimplified' 'SignalProcessing',{'Resampling',100,'EpochExtraction',[0 0.8],'WaveletTransform','on','ICA',{'Variant','beamica','TransformData',true}}, ...
    'Prediction',{'MachineLearning',{'Learner',{'logreg','variant','lars'}}}};
% again wavelets, but this time group-sparse logistic regression (takes ~50 minutes!!!)
approaches.waveletsglr = {'DataflowSimplified' 'SignalProcessing',{'Resampling',100,'EpochExtraction',[0 0.8],'WaveletTransform','on'}, ...
    'Prediction',{'MachineLearning',{'Learner',{'dal',2.^(10:-0.125:1),'scaling','std','regularizer','glc','shape',[64 NaN]}}}};
% logistic regression with group-sparse Laplace prior over sphering components as well as sparsity in the wavelet domain (hardcoded hyper-parameters)
approaches.compwavelets = {'DataflowSimplified' 'SignalProcessing',{'Resampling',128,'EpochExtraction',[-0.2 0.79],'ICA',{'Variant','robust_sphere','DataCleaning','off'}}, ...
    'Prediction',{'FeatureExtraction',{'GroupInto','matrix'},'MachineLearning',{'Learner',{'glm', ...
        'ptype','classification', ...
        'lambdas',1, ...
        'priors',{ ...
            'term1',{'Laplace','LinearOperator','@(x)x''','Scales',80}, ...
            'term2',{'Laplace','LinearOperator','@(x)vec(mydwt(x,daubcqf(4,''min''),2))''','Scales',15}, ...
    }}}}};
% again wavelets, using group-sparse logistic regression and an elastic-net like extra penalty (using the proximal classifier)
% cross-validdation takes ca. 8hrs; demo only for syntax
approaches.waveletsprox = {'DataflowSimplified' 'SignalProcessing',{'Resampling',128,'EpochExtraction',[0 0.8],'WaveletTransform','on'}, ...
    'Prediction',{'FeatureExtraction',{'GroupInto','matrix'}, ...
        'MachineLearning',{'Learner',{'proximal' 'Regularizers',{'Term1','l1/l2','Term2','l2'}}}}};


% for each of the above approaches...
for app = fieldnames(approaches)'
    fprintf(['\n==== now testing "' app{1} '" ====\n']);
    fprintf([utl_printapproach(approaches.(app{1})) '\n\n']);
    % train & cross-validate
    [trainloss,lastmodel,laststats] = bci_train('Data',traindata,'Approach',approaches.(app{1}),'TargetMarkers',{{'S101','S102'},{'S201','S202'}})
    disp(['training mis-classification rate: ' num2str(trainloss*100,3) '%']);
    % test offline
    [prediction,loss,teststats,targets] = bci_predict(lastmodel,testdata);
    disp(['offline test mis-classification rate: ' num2str(loss*100,3) '%']);
    % test pseudo-online
    [predictions,latencies] = onl_simulate(testdata,lastmodel,'markers',{'S101','S102','S201','S202'},'offset',0.8);
    disp(['pseudo-online test mis-classification rate: ' num2str(mean(argmax(predictions') ~= targets')*100,3) '%']);
    % visualize in real time
    run_readdataset('Dataset',testdata); run_writevisualization('Model',lastmodel, 'VisFunction','bar3(y)'); 
    waitforbuttonpress; onl_clear; close(gcf);
end
