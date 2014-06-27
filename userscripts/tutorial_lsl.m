%% === From another MATLAB: stream EEG and Markers via the lab streaming layer ===

% NOTE: instead of playing back data you would normally turn on an LSL application
% that streams data from some piece of hardware
testdata = io_loadset('data:/tutorial/flanker_task/12-08-001_ERN.vhdr');
play_eegset_lsl(testdata);


%% === calibrate a BCI model ===

% define markers; here, two groups of markers are being defined; the first group represents class 1
% (correct responses), and the second group represents class 2 (incorrect responses).
mrks = {{'S101','S102'},{'S201','S202'}};

% define ERP windows of interest; here, 7 consecutive windows of 50ms length each are being
% specified, starting from 250ms after the subject response
wnds = [0.25 0.3;0.3 0.35;0.35 0.4; 0.4 0.45;0.45 0.5;0.5 0.55;0.55 0.6];

% load your calibration data here
% NOTE: typically this is a calibration recording that you have previously recorded
traindata = io_loadset('data:/tutorial/flanker_task/12-08-002_ERN.vhdr');

% define an approach
myapproach = {'Windowmeans' 'SignalProcessing', {'Resampling','off','EpochExtraction',[-0.2 0.8],'SpectralSelection',[0.1 15]}, 'Prediction',{'MachineLearning',{'Learner',{'lda',0.1,'regularization','shrinkage'}},'FeatureExtraction',{'TimeWindows',wnds}}};

% learn model 
tic;[trainloss,lastmodel,laststats] = bci_train('Data',traindata,'Approach',myapproach,'TargetMarkers',mrks);toc
disp(['training mis-classification rate: ' num2str(trainloss*100,3) '%']);


%% === visualize continuous BCI outputs from EEG stream on LSL ===

% read from an any EEG device on the network, do not require a marker stream
run_readlsl('DataStreamQuery','type=''EEG''', 'MarkerQuery','');

% process it in real time using lastmodel, and visualize outputs
run_writevisualization('Model',lastmodel, 'VisFunction','bar(y);ylim([0 1])');

% make sure that the online processing gets terminated...
disp('Click into the figure to stop online processing.'); 
waitforbuttonpress; onl_clear; close(gcf);

%% === stream continuous BCI outputs over LSL ===

% read from an any EEG device on the network, do not require a marker stream
run_readlsl('DataStreamQuery','type=''EEG''', 'MarkerQuery','');

% process it in real time and emit as stream named "BCI-Continuous"
% this stream can be read anywhere on the network (using LSL for MATLAB, Python, C, C++, Java or C#)
run_writelsl('Model',lastmodel,'LabStreamName','BCI-Continuous');

%% === From a third MATLAB: receive and display BCI outputs ===

% alternatively use 'BCI-At-Markers' to display the marker-locked BCI stream
bci_stream_name = 'BCI-Continuous';  

f=figure;
lib = lsl_loadlib();
disp('Resolving a BCI stream...');
result = {};
while isempty(result)
    result = lsl_resolve_byprop(lib,'name',bci_stream_name,1,1); end
inlet = lsl_inlet(result{1});
disp('Now receiving data...');
while true
    % get data from the inlet (timeout: 1 second)
    [data,timestamp] = inlet.pull_sample(0);
    % and display it
    if timestamp
        fprintf('%.2f\n',data); 
        bar([data-1,1-(data-1)]); ylim([0 1]); drawnow;
    else
        pause(0.01);
    end
end


%% === calculate BCI outputs *only* at markers in a marker stream ===

% read from an any EEG device on the network, and also acquire a marker stream
run_readlsl('DataStreamQuery','type=''EEG''', 'MarkerQuery','type=''Markers''');

% this time predict whenever there is a marker of the given type in the marker stream
run_writevisualization('PredictorName','pred1', 'Model',lastmodel, 'PredictAt',{'S101','S102','S201','S202'}, 'VisFunction','bar(y);ylim([0 1])');
% also send the results out via LSL (note: we're here letting BCILAB do the calculations twice, not the most efficient way)
run_writelsl('PredictorName','pred2', 'Model',lastmodel, 'PredictAt',{'S101','S102','S201','S202'}, 'LabStreamName','BCI-At-Markers');

% make sure that the online processing gets terminated...
disp('Click into the figure to stop online processing.'); 
waitforbuttonpress; onl_clear; close(gcf);

