%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SCRIPTING EXAMPLE FOR THE SOURCE INFORMATION FLOW TOOLBOX (SIFT)    %%%
%%% SIFT Version: 0.7-alpha                                             %%%
%%%                        
%%% THIS DOES NOT WORK WITH VERSION 0.1-BETA OR LATER
%%%
%%% This example demonstrates how to use SIFT from the command-line or  %%%
%%% in a script. This example applies to SIFT 0.7-alpha.                %%%
%%% For additional information on the below steps, please consult the   %%%
%%% SIFT manual located at http://sccn.ucsd.edu/wiki/SIFT               %%%
%%% Author: Tim Mullen (C) 2011, SCCN, INC, UCSD                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% STEP 1: Load Data

% We will begin by loading up the 'RespWrong.set' and 'RespCorr.set' datasets located in the 
% /Data/ folder within the Sample Data package
% (you can download this package from the SIFT website or at 
% ftp://sccn.ucsd.edu/pub/tim/SIFT/SIFT_SampleData.zip)

EEG = pop_loadset;

%% OPTIONAL STEP: Analyzing Channels

% If you want to analyze channel data rather than components you can use
% this "hack". Namely, we will replace the ICA soution with a "fake" 
% solution that always copies the channel data into the ICA component
% (icaact) field. This ensures that SIFT will always work on the channel
% data rather than component data. NOTE: this functionality will be
% replaced by a more flexible interface in the Beta version.

h=msgbox('Next the EEGLAB options GUI will pop up. Uncheck the option labeled "If set, scale ICA components to RMS in microvolt" and click OK');
waitfor(h);
pop_editoptions;


% Now we will replace the ICA soution with a "fake" solution that always
% copies the channel data into the ICA component (icaact) field
% WARNING, THIS WILL *DELETE* YOUR ICA SOLUTION, uncomment the following
% lines to create a backup of the current ICA solution

% EEG.etc.icaweights_beforeIdentity = EEG.icaweights;
% EEG.etc.icasphere_beforeIdentity = EEG.icasphere;
% EEG.etc.icawinv_beforeIdentity = EEG.icawinv;

[EEG.icaweights EEG.icasphere EEG.icawinv] = deal(eye(EEG.nbchan));
EEG.icaact = [];

% copy the data into ICA 'activations' and verify dataset
EEG = eeg_checkset(EEG,'ica');


%% STEP 2: Define key Processing Parameters

Components = [8 11 13 19 20 23 38 39];   % these are the components/channels to which we'll fit our multivariate model
WindowLengthSec = 0.4;                   % sliding window length in seconds
WindowStepSizeSec = 0.03;                % sliding window step size in seconds
NewSamplingRate = [];                    % new sampling rate (if downsampling)
EpochTimeRange = [-1 1.25];              % this is the time range (in seconds) to analyze (relative to event at t=0)

%% STEP 3: Pre-process the data

ComponentsToKeep = strtrim(cellstr(num2str(Components')));  % convert list of components to cell array of strings

[EEG prepcfg] = pre_prepData('ALLEEG',EEG,'VerbosityLevel',2,'NewSamplingRate',NewSamplingRate,'EpochTimeRange',EpochTimeRange,'NormalizeData',{'Method',{'ensemble'}},'SelectComponents',{'ComponentsToKeep',ComponentsToKeep});

%% STEP 4: Identify the optimal model order

% Here we compute various model order selection criteria for varying model
% orders (e.g. 1 to 30) and visualize the results

% compute model order selection criteria...
IC = pop_est_selModelOrder(EEG,0,'icselector',{'aic','sbc','hq','ris'},'algorithm','vieira-morf','morder',[1 30],'winlen',WindowLengthSec,'winstep',WindowStepSizeSec,'prctWinToSample',100,'verb',1,'plot',0);

% ... and plot the results
handles = vis_plotOrderCriteria(IC);

% If you want to save this figure you can uncomment the following lines:
%
% for i=1:length(handles)
%     saveas(handles(i),sprintf('orderResults%d.fig',i));
% end
% close(handles);

% Finally, we can automatically select the model order which minimizes one
% of the criteria (or you can set this manually based on above figure)
ModelOrder = ceil(mean(IC{1}.hq.popt));

% As an alternative to using the minimum of the selection criteria over 
% model order, you can find the "elbow" in the plot of model order versus
% selection criterion value. This is useful in cases where the selection
% criterion does not have a clear minimum. For example, the lines below
% plot and select the elbow location (averaged across windows) for the AIC 
% criterion
%
% vis_plotOrderCriteria(IC,{},{},'elbow');
% ModelOrder = ceil(mean(IC{1}.aic.pelbow));



%% STEP 5: Fit the VAR model

% Here we can check that our selected parameters make sense
fprintf('===================================================\n');
fprintf('MVAR PARAMETER SUMMARY FOR CONDITION: %s\n',EEG.condition);
fprintf('===================================================\n');
est_dispMVARParamCheck(EEG,struct('morder',ModelOrder','winlen',WindowLengthSec,'winstep',WindowStepSizeSec,'verb',1));

% Once we have identified our optimal model order, we can fit our VAR model.

% Fit a model using a sliding-window approach with the vieira-morf
% lattice filter algorithm
[EEG modfitcfg] = pop_est_fitMVAR(EEG,0,'algorithm','vieira-morf','morder',ModelOrder,'winlen',WindowLengthSec,'winstep',WindowStepSizeSec,'verb',1);


% Note that EEG.CAT.MODEL now contains the model structure with
% coefficients (in MODEL.AR), prediction errors (MODEL.PE) and other
% self-evident information

% Alternately, we can fit the VAR parameters using a Kalman filter (see
% doc est_fitMVARKalman for more info on arguments)
%
% EEG.CAT.MODEL = est_fitMVARKalman(EEG,0,'updatecoeff',0.0005,'updatemode',2,'morder',ModelOrder,'verb',2,'downsampleFactor',50);


%% STEP 6: Validate the fitted model

% Here we assess the quality of the fit of our model w.r.t. the data. This
% step can be slow.

% We can obtain statistics for residual whiteness, percent consistency, and
% model stability ...
[whitestats PC stability] = pop_est_validateMVAR(EEG,0,'checkWhiteness',true,'whitenessCriteria',{'Ljung-Box','ACF','Box-Pierce','Li-McLeod'},...
                                                'checkConsistency',true,'checkStability',true,'alpha',0.05,'prctWinToSample',20,'verb',2,'plot',0);

% ... and then plot the results
handles = vis_plotModelValidation(whitestats,PC,stability);

% If you want to save this figure you can uncomment the following lines:
%
% for i=1:length(handles)
%     saveas(handles(i),sprintf('validationResults%d.fig',i));
% end
% close(handles);


% To automatically determine whether our model accurately fits the data you
% can write a few lines as follows (replace 'acf' with desired statistic):
%
% if ~all(whitestats{1}.acf.w)
%     msgbox('Residuals are not completely white!');
% end


%% STEP 7: Compute Connectivity

% Next we will compute various dynamical quantities, including connectivity,
% from the fitted VAR model. We can compute these for a range of
% frequencies (here 3-45 Hz). See 'doc est_mvarConnectivity' for a complete
% list of available connectivity and spectral estimators.

[EEG conncfg] = pop_est_mvarConnectivity(EEG,'verb',true,'freqs',(1 : 50),'connmethods',{'dDTF08','S','nPDC','RPDC'},'absvalsq',true,'spectraldecibels',true);


%% OPTIONAL STEP 8: Compute Statistics (This step is slow)

NumPermutations = 200;

% reload each of the datasets in the exact same order they appear in ALLEEG
% (or use original, un-preprocessed copy)
for cnd=1:length(EEG)
    EEGfresh(cnd) = pop_loadset;
end

% first we obtain the bootstrap distributions for each condition
for cnd=1:length(EEG)
    PConn_boot(cnd) = stat_surrogate('ALLEEG',EEGfresh(cnd),'configs',struct('prepData',prepcfg(1),'fitMVAR',modfitcfg,'mvarConnectivity',conncfg),'Mode',{'Bootstrap', 'NumPermutations', NumPermutations},'AutoSave',{'FileNamePrefix','SIFT_bootstrap','AutoSaveFrequency',10},'VerbosityLevel',2);
end

% replace connectivity object with estimate of bootstrap mean
for cnd=1:length(EEG)
    EEG(cnd).CAT.Conn = stat_getDistribMean(PConn_boot(cnd));
end

%% NOTE: we can also obtain the phase-randomized null distributions for each condition
for cnd=1:length(EEG)
    PConn_phase(cnd) = stat_surrogate('ALLEEG',EEGfresh(cnd),'configs',struct('prepData',prepcfg(1),'fitMVAR',modfitcfg,'mvarConnectivity',conncfg),'Mode',{'PhaseRand', 'NumPermutations', NumPermutations},'AutoSave',{'FileNamePrefix','SIFT_bootstrap','AutoSaveFrequency',10},'VerbosityLevel',2);
end

%% next we compute p-values and confidence intervals
% (CHOOSE ONE OF THE FOLLOWING)

%% 1) Between-condition test:
%     For conditions A and B, the null hypothesis is either
%     A(i,j)<=B(i,j), for a one-sided test, or
%     A(i,j)=B(i,j), for a two-sided test
%     A p-value for rejection of the null hypothesis can be
%     obtained by taking the difference of the distributions
%     computing the probability
%     that a sample from the difference distribution is non-zero
StatsHab = stat_surrogateStats('BootstrapConnectivity',PConn_boot,'StatisticalTest',{'Hab'},'MultipleComparisonCorrection','none','ConfidenceIntervals',true,'Alpha',0.05,'VerbosityLevel',1);

%% 2) Devation from baseline test
%     For conditions A, the null hypothesis is
%     C(i,j)=baseline_mean(C). This is a two-sided test.
%     A p-value for rejection of the null hypothesis can be
%     obtained by obtaining the distribution of the difference from
%     baseline mean and computing the probability
%     that a sample from this distribution is non-zero
for cnd=1:length(EEG)
    StatsHbase(cnd) = stat_surrogateStats('BootstrapConnectivity',PConn_boot(cnd),'StatisticalTest',{'Hbase' 'Baseline', [-1 -0.25]},'MultipleComparisonCorrection','none','ConfidenceIntervals',true,'Alpha',0.05,'VerbosityLevel',1);
end

%% 3) Presence of absolute connectivity
%     We are testing with respect to a phase-randomized null
%     distribution. A p-value for rejection of the null hypothesis
%     can be obtained by computing the probability that the
%     observed connectivity is a random sample from the null distribution
for cnd=1:length(EEG)
    StatsHnull(cnd) = stat_surrogateStats('BootstrapConnectivity',EEG(cnd).CAT.Conn,'NullDistribution',PConn_phase(cnd),'StatisticalTest',{'Hnull'},'MultipleComparisonCorrection','none','ConfidenceIntervals',true,'Alpha',0.05,'VerbosityLevel',1);
end

%% OPTIONAL STEP 8b: Compute analytic statistics
% This computes analytic alpha-significance thresholds, p-values, and confidence
% intervals for select connectivity estimators (RPDC, nPDC).
% These are asymptotic estimators and may not be accurate for small sample
% sizes. However, they are very fast and usually a reasonable estimate.
StatsAnalytic = stat_analyticStats('ALLEEG',EEG,'Estimator',{'RPDC','nPDC'},'Alpha', 0.01,'verb',true);


%% STEP 9: Visualize the Connectivity estimates in a Time-Frequency Grid

% this plots the difference between conditions
[figureHandles tfgridcfg] = vis_TimeFreqGrid('ALLEEG',EEG,'Conn',[EEG(1).CAT.Conn EEG(2).CAT.Conn],'MatrixLayout',{'Partial','UpperTriangle', 'dDTF08', 'LowerTriangle','dDTF08','Diagonal','S'},'ColorLimits',99.9,'FrequencyScale','log','Baseline',[],'Smooth2D',false);

%% Or if stats are present we can threshold between-cond difference
[figureHandles tfgridcfg] = vis_TimeFreqGrid('ALLEEG',EEG,'Conn',[EEG(1).CAT.Conn EEG(2).CAT.Conn],'Stats',StatsHab,'MatrixLayout',{'Partial','UpperTriangle', 'dDTF08', 'LowerTriangle','dDTF08','Diagonal','S'},'ColorLimits',99.9,'FrequencyScale','log','Baseline',[],'Smooth2D',false,'Thresholding',{'Statistics','ThresholdingMethod','pval','PlotConfidenceIntervals',true},'BackgroundColor',[0 0 0],'TextColor',[1 1 1],'LineColor',[1 1 1]);

%% this plots the deviation from baseline (with stats)
for cnd=1:length(EEG)
    [figureHandles tfgridcfg] = vis_TimeFreqGrid('ALLEEG',EEG(cnd),'Conn',EEG(cnd).CAT.Conn,'Stats',StatsHbase(cnd),'MatrixLayout',{'Partial','UpperTriangle', 'dDTF08', 'LowerTriangle','dDTF08','Diagonal','S'},'ColorLimits',99.9,'FrequencyScale','log','Baseline',[-1.75 -0.5],'Smooth2D',false,'Thresholding',{'Statistics','ThresholdingMethod','pval','PlotConfidenceIntervals',true},'BackgroundColor',[0 0 0],'TextColor',[1 1 1],'LineColor',[1 1 1]);
end

%% this plots a single condition with phase-randomization thresholding
for cnd=1:length(EEG)
    [figureHandles tfgridcfg] = vis_TimeFreqGrid('ALLEEG',EEG(cnd),'Conn',EEG(cnd).CAT.Conn,'Stats',StatsHnull(cnd),'MatrixLayout',{'Partial','UpperTriangle', 'dDTF08', 'LowerTriangle','dDTF08','Diagonal','S'},'ColorLimits',99.9,'FrequencyScale','log','Baseline',[],'Smooth2D',false,'Thresholding',{'Statistics','ThresholdingMethod','pval','PlotConfidenceIntervals',true},'BackgroundColor',[0 0 0],'TextColor',[1 1 1],'LineColor',[1 1 1]);
end
% You can also partially populate the GUI via a call to the pop_ function:
%
%[figureHandles tfgridcfg] = pop_vis_TimeFreqGrid(EEG,'MatrixLayout',{'Partial','UpperTriangle', 'dDTF08', 'LowerTriangle','dDTF08','Diagonal','S'},'ColorLimits',99.9,'Baseline',[-1.75 -0.5],'Smooth2D',true);


%% STEP 10: Visualize the Connectivity estimates in a 3D Brain-Movie

cfg=pop_vis_causalBrainMovie3D(EEG(1),'ConnectivityMethod','dDTF08','FreqCollapseMethod','integrate','FrequenciesToCollapse',(3:8),'NodeColorMapping','CausalFlow','FooterPanelDisplaySpec',{'ICA_ERPenvelope',{'icaenvelopevars', '8'},{'backprojectedchans' 'B2'}},'BrainMovieOptions',{'RenderCorticalSurface',{'VolumeMeshFile' 'standard_BEM_vol.mat', 'Transparency' 0.7}});    %


