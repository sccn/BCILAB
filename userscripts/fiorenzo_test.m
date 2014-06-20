%% --- offline processing ---

% clear cache if anything goes wrong!
env_clear_memcaches

% define load training data (BrainVision format)
traindata = io_loadset('data:/tutorial/flanker_task/12-08-002_ERN.vhdr');
% continuous!!!
% markers included!!!



% define markers; here, two groups of markers are being defined; the first group represents class 1
% (correct responses), and the second group represents class 2 (incorrect responses).
mrks = {{'S101','S102'}, {'S201','S202'}};

% define a simple approach (we'll change)
%myapproach = {'DALERP','SignalProcessing',{'Resampling',60,'IIRFilter','off','spectrum',[0.1 15],'EpochExtraction',[0 0.8]}, ...
%    'Prediction',{'MachineLearning',{'Learner',{'dal','solver','cg'}}}};
myapproach = {'Timewarp','SignalProcessing',{'Resampling',60}};

% learn model & evaluate
[trainloss,lastmodel,laststats] = bci_train('Data',traindata,'Approach',myapproach,'TargetMarkers',mrks,'EvaluationMetric','auc')

% visualize classier details
lastmodel.predictivemodel.model
vertcat(laststats.per_fold.pred)
figure;plot(vertcat(laststats.per_fold.targ));hold;plot(vertcat(laststats.per_fold.pred),'r')

%% --- online processing ---v
% see tutorial_lsl.m

