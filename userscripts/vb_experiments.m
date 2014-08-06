%% initialize 

% define markers and informative time windows
mrks = {{'S101','S102'},{'S201','S202'}};
wnds = [0.25 0.3;0.3 0.35;0.35 0.4; 0.4 0.45;0.45 0.5;0.5 0.55;0.55 0.6];
chans64 = {'Fp1', 'Fp2', 'F7', 'F3', 'Fz', 'F4', 'F8', 'FC5', 'FC1', 'FC2', 'FC6', 'T7', 'C3', 'Cz', 'C4', 'T8', 'TP9', 'CP5', 'CP1', 'CP2', 'CP6', 'TP10', 'P7', 'P3', 'Pz', 'P4', 'P8', 'PO9', 'O1', 'Oz', 'O2', 'PO10', 'Fpz', 'AF7', 'AF3', 'AF4', 'AF8', 'F5', 'F1', 'F2', 'F6', 'FT9', 'FT7', 'FC3', 'FC4', 'FT8', 'FT10', 'C5', 'C1', 'C2', 'C6', 'TP7', 'CP3', 'CPz', 'CP4', 'TP8', 'PO7', 'PO3', 'POz', 'PO4', 'PO8', 'RE', 'LE', 'VEOG'};
chans20 = {'F3', 'Fz', 'FC5', 'FC1', 'FC2', 'FC6', 'Cz', 'C4', 'CP5', 'CP1', 'P8', 'AF4', 'F1', 'F2', 'F6', 'FC4', 'C5', 'C1', 'C2', 'CP3'};
chans13 = {'Fz', 'FC2', 'FC6', 'Cz', 'C4', 'CP5', 'F1', 'F2', 'F6', 'FC4', 'C1', 'C2', 'CP3'};
chans = chans13;

% load data
traindata = io_loadset('bcilab:/userdata/tutorial/flanker_task/12-08-002_ERN.vhdr');

% preprocess

flt = exp_eval_optimized(set_picktimes(set_makepos(set_targetmarkers(flt_iir(flt_fir(flt_resample(flt_selchans(traindata,chans),60),[0.5 1],'highpass'),[15 20],'lowpass'),mrks),[-0.2 0.8]),wnds));

%% get predictors and target variable
[C,S,T] = size(flt.data); F = C*S;  % #channels, #samples, #trials, #features
Xt = flt.data;                      % tensor-shaped data (#channels x #time-points x #trials)
Xm = reshape(Xt,[],size(Xt,3))';     % matrix-shaped data (#trials x #features)
y = [flt.epoch.target]'*2-3;         % target variable (1 x #trials), \in {-1,+1}

% fixed hyper-parameters: s^2_y = 1 (observation noise) and s^2_beta (prior uncertainty)
s2y = 1; s2beta = 1000;


%% --- using vb_mvn ---
% fixed Gaussian prior for beta
% log p(beta|y) = mvnormpdfln(y,Xm*beta(2:end)+beta(1),s2y*eye(T)) + mvnormpdfln(beta,zeros(F+1,1),s2b*eye(F+1)));
vb_mvn({'beta',zeros(F+1,1)},{@(beta)logmvnpdf(y,Xm*beta(2:end)+beta(1),s2y*eye(T)),@(beta)logmvnpdf(beta,zeros(F+1,1),s2beta*eye(F+1))});

vb_mvn('Variables', {'beta',zeros(F+1,1)}, ...
       'Factors', {vb_factor(@logmvnpdf,y,@(beta)beta(1)+Xm*beta(2:end),s2y*eye(T)), ...
                   vb_factor(@logmvnpdf,@(beta)beta,zeros(F+1,1),s2beta*eye(F+1))});

               
%{@(beta)logmvnpdf(y,beta(1) + Xm*beta(2:end),s2y*eye(T)),@(beta)logmvnpdf(beta,zeros(F+1,1),s2beta*eye(F+1))});

%% --- model 1 ---
disp('solving...');
% standard regression

% s2y = 1, s2b
% beta ~ N(0,s2b);
% y ~ N(Xm*beta(2:end)+beta(1),s2y);
% p(beta|y) ~~ p(y|beta) * p(beta) / Z
%            = N(y|Xm*beta(2:end) + beta(1),s2y) * N(beta|0,s2b)
% log p(beta|y) = mvnormpdfln(y,Xm*beta(2:end)+beta(1),s2y*eye(F)) + mvnormpdfln(beta,zeros(F+1,1),s2b*eye(F+1)));

beta = zeros(F+1,1);


model_logpdf = @(beta) mvnormpdfln(y,Xm*beta(2:end)+beta(1),s2y*eye(T)) + mvnormpdfln(beta,zeros(F+1,1),s2b*eye(F+1));
[mu,prec] = GaussVarApproxBasic_work(zeros(F+1,1),eye(F+1),model_logpdf,100);
disp('done.');

model_logpdf = @(beta) mvnormpdfln(y,Xm*beta(2:end)+beta(1),s2y*eye(T)) + mvnormpdfln(beta,zeros(F+1,1),s2b*eye(F+1));
[mu,prec] = GaussVarApproxBasic_work(zeros(F+1,1),eye(F+1),model_logpdf,100);


%% test some stats
d = 100; n=1000;
s = RandStream('mt19937ar','Seed',0);
X = randn(s,d,n); M = randn(s,d,1); U=randn(d); S=U*U';

mvnormpdfln(X,M,[],S) - logmvnpdf(X,M,S)

%% plot posterior mean model
figure;topoplot_grid(reshape(mu(2:end),[length(chans) length(wnds)]),set_infer_chanlocs(chans));


%% junk
         % log((2*pi)^-k/2*det(sigma)^-1/2))  -1/2 * (y-(Xm*beta+alpha)) * diag(1/s2y) * (y-(Xm*beta+alpha))' + -1/2*

%% --- traditional analysis ---

% define markers; here, two groups of markers are being defined; the first group represents class 1
% (correct responses), and the second group represents class 2 (incorrect responses).

% define ERP windows of interest; here, 7 consecutive windows of 50ms length each are being
% specified, starting from 250ms after the subject response

% define load training data (BrainVision format)
traindata = exp_eval(io_loadset('bcilab:/userdata/tutorial/flanker_task/12-08-002_ERN.vhdr'));

% define approach
myapproach = {'Windowmeans' 'SignalProcessing', {'Resampling','off','EpochExtraction',[-0.2 0.8],'SpectralSelection',[0.1 15]}, 'Prediction',{'FeatureExtraction',{'TimeWindows',wnds}}};

%learn model 
[trainloss,lastmodel,laststats] = bci_train('Data',traindata,'Approach',myapproach,'TargetMarkers',mrks);
disp(['training mis-classification rate: ' num2str(trainloss*100,3) '%']);

% visualize results
bci_visualize(lastmodel)


% get most informative channels
W=reshape(lastmodel.predictivemodel.model.w,[64 7]);
[dum,idx]=sort(abs(W),'descend');
numrows = 5;
bestchans = {traindata.chanlocs(unique(vec(idx(1:numrows,:)))).labels};
hlp_tostring(bestchans)