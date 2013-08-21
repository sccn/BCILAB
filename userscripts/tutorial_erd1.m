%% --- tutorial script for spectral analysis ---

% This tutorial demonstrates the use of BCILAB to learn basic predictive models using 
% spectral properties of the data, here operating on motor cortex idle oscillations.
%
% The used data set contains a sequence of trials in which a subject was instructed to imagine
% moving either the left hand or the right hand. Markers in the data set indicate the timing and 
% type of these instructions (types are: 'StimulusCode_2' and 'StimulusCode_3').
%
% Notes: StimulusCode_1 indicates onset of a "resting" / "relaxing" condition.
%
% Data courtesy of Romain Grandchamp, CERCO, Toulouse, France. Reference:
%   Romain Grandchamp, Arnaud Delorme, "NeuroTRIP: A Framework for Bridging between Open Source Software. Application to Training a Brain Machine Interface," 
%   Fifth International Conference on Signal Image Technology and Internet Based Systems, pp.451-457, 2009
%
%#ok<*ASGLU,*NASGU,*SNASGU> % turn off a few editor warnings...

%% --- using the Common Spatial Pattern method ---

% load the data set (BCI2000 format)
traindata = io_loadset('data:/tutorial/imag_movements1/calib/DanielS001R01.dat');

% define the approach (here: Common Spatial Patterns without any customization)
myapproach = 'CSP';

% learn a predictive model
[trainloss,lastmodel,laststats] = bci_train('Data',traindata,'Approach',myapproach,'TargetMarkers',{'StimulusCode_2','StimulusCode_3'}); 
disp(['training mis-classification rate: ' num2str(trainloss*100,3) '%']);

% visualize results
bci_visualize(lastmodel)

%% --- using the Common Spatial Pattern method with some custom options ---

% load the data set (BCI2000 format)
traindata = io_loadset('data:/tutorial/imag_movements1/calib/DanielS001R01.dat');

% define the approach 
% Note: The settings found in the GUI "Review/Edit Approach" Panel can be translated literally
%       into cell array representations as below. Each paradigm has a few top-level parameter groups
%       (for CSP: SignalProcessing, FeatureExtraction, etc), which in turn have sub-parameters
%       (e.g., SignalProcessing has EpochExtraction, SpectralSelection, Resampling, etc.), and so
%       on. Some parameters are numbers, strings, cell arrays, etc. You only need to specify those 
%	    parameters where you actually want to deviate from the paradigm's defaults.
%
% For illustratory purposes, we use a different window relative to the target markers (0.5s to 3s after),
% and a somewhat customized FIR frequency filter with a pass-band between ~7.5Hz and ~27Hz.
myapproach = {'CSP' 'SignalProcessing',{'EpochExtraction',[0.5 3],'FIRFilter',[7 8 26 28]}};

% learn a predictive model
[trainloss,lastmodel,laststats] = bci_train('Data',traindata,'Approach',myapproach,'TargetMarkers',{'StimulusCode_2','StimulusCode_3'}); 
disp(['training mis-classification rate: ' num2str(trainloss*100,3) '%']);

% visualize results
bci_visualize(lastmodel)


%% --- applying the CSP model to some data (here: same data) ---

% apply the previously learned model to a data set (querying it for each target marker in the data)
[prediction,loss,teststats,targets] = bci_predict(lastmodel,traindata);

% display the results
disp(['test mis-classification rate: ' num2str(loss*100,3) '%']);
disp(['  predicted classes: ',num2str(round(prediction{2}*prediction{3})')]);  % class probabilities * class values
disp(['  true classes     : ',num2str(round(targets)')]);


%% --- the same using pseudo-online processing ---

% This function is quite flexible, but here we only ask it to query the BCI seconds after each of the given markers
[predictions,latencies] = onl_simulate(traindata,lastmodel,'markers',{'StimulusCode_2','StimulusCode_3'},'offset',3);

% now check the prediction accuracy by hand (we knew the correct target label at each of the markers)
accuracy = mean(argmax(predictions') == targets');
disp(['pseudo-online mis-classification rate: ' num2str((1-accuracy)*100,3) '%']);


%% --- using the Spectrally weighted Common Spatial Pattern (Spec-CSP) method ---
% this method automatically learns the spectral filters (within a pre-defined range), but it 
% may run into a local optimum or over-fit spuriously correlated bands

% load the data set
traindata = io_loadset('data:/tutorial/imag_movements1/calib/DanielS001R01.dat');

% define the approach
myapproach = {'SpecCSP' 'SignalProcessing',{'EpochExtraction',[0.5 3]}};

% learn a predictive model
[trainloss,lastmodel,laststats] = bci_train('Data',traindata,'Approach',myapproach,'TargetMarkers',{'StimulusCode_2','StimulusCode_3'}); %#ok<>

% visualize results
bci_visualize(lastmodel)


%% --- test the learned model in simulated real-time processing ---
% ( click into the figure to stop the update (and make sure that your click was registered) )

% load feedback session
testdata = io_loadset('data:/tutorial/imag_movements1/feedback/DanielS001R01.dat');

% play it back in real time
run_readdataset('Dataset',testdata);

% process data in real time using lastmodel, and visualize outputs
run_writevisualization('Model',lastmodel, 'VisFunction','bar(y)');

% make sure that the online processing gets terminated...
disp('Click into the figure to stop online processing.'); 
waitforbuttonpress; onl_clear; close(gcf);


%% --- train an alternative model with parameter search ---
% (over possible values for the number of pattern pairs, using CSP; note: this takes quite some time!)
% (the number of pattern pairs found optimal should be 3 in this case)

% load the data set (BCI2000 format)
traindata = io_loadset('data:/tutorial/imag_movements1/calib/DanielS001R01.dat');

% define approach
myapproach = {'CSP' 'Prediction',{'FeatureExtraction',{'PatternPairs',search(1,2,3)}}};

% learn the model
[trainloss,lastmodel,laststats] = bci_train('Data',traindata,'Approach',myapproach,'TargetMarkers',{'StimulusCode_2','StimulusCode_3'});
disp(['training mis-classification rate: ' num2str(trainloss*100,3) '%']);

% visualize results
bci_visualize(lastmodel);
