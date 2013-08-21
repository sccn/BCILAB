%% === CONSTANTS ===
timewnd = [0.5 2.5]; 
mrks = {{'769','1'},{'770','2'}};
trainsets = {'data://bcicomp3/IIIa/*_train.set'};
evalsets = {'data://bcicomp3/IIIa/*_test.set'};

% === DEFINE APPROACHES ===
% (note: use unique names for all your approaches)

approaches = [];
approaches.BP_wideband_LDA_3ch = {'Bandpower' 'SignalProcessing',{'ChannelSelection',{{'C3','Cz','C4'}}, 'EpochExtraction',timewnd}};
approaches.CSP_wideband_shrinklda = {'CSP' 'SignalProcessing',{'EpochExtraction',timewnd}};
approaches.CSP_wideband_qda = {'CSP' 'SignalProcessing',{'EpochExtraction',timewnd},'Prediction',{'MachineLearning',{'Learner','qda'}}};

% === RUN BATCH ANALYSIS ===

results = bci_batchtrain('StudyTag','BCICompA','Data',trainsets,'PredictSets',evalsets,'Approaches',approaches,'TargetMarkers',mrks,'ReuseExisting',true, ...
    'LoadArguments',{'type','EEG'}, 'TrainArguments',{'EvaluationScheme',{'chron',10,1}},'StoragePattern','batchresults/%caller-%study/%approach-%set.mat');

