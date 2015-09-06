% Z_mean{i}: [num_components_for_subj_i x num_clusters]
%            contains posterior prob. of cluster membership
% cmp_keep:  [1 x num_subjects] list of preserved components for each
%           subject
% B_bar_mean{i,j}:  B-spline basis coefficients from cluster i to cluster j
%                   (note this is the transpose of my normal format)
% phi_t:        [num_b_splines x time] time-varying b-spline basis function
% Theta_diag:   [num_b_splines x num_PCs] FPCA weights for
%               self-connectivity
% Theta_offdiag: [num_b_splines x num_PCs] FPCA weights for
%               cross-connectivity
% NOTE: To obtain the 80x1 time-varying group PDC from cluster i to cluster
%       j -- F(i,j):
%       F(j,j) = phi_t'Theta_diag'B_bar_mean{j,j}
%       F(i,j) = phi_t'Theta_offdiag'B_bar_mean{i,j}  where j ~= i
%
% S_bar_mean:   [num_clusters x 3] MNI centroids for each cluster
%
%


%%
clear all; clc;

load real_data
load smoothing_outputs
load mcmc_outputs

%% Get posterior means spatial locations & connectivity coefs & indicators
S_bar_mean=zeros(M,3);
for i=1:size(S_BAR_array,1)
    S_bar_mean=S_bar_mean+S_BAR_array{i}/size(S_BAR_array,1);
end
B_bar_mean=cell(M,M);
for k1=1:M
    for k2=1:M
        B_bar_mean{k1,k2}=zeros(1,Q);
        for iter=1:size(B_BAR_array,1)
            B_bar_mean{k1,k2}=B_bar_mean{k1,k2}+B_BAR_array{iter}{k1,k2}/size(B_BAR_array,1);
        end
    end
end
Z_mean=cell(N,1);
for i=1:N
    Z_mean_i=Z_array{1}{i};
    for iter=2:size(Z_array,1)
        Z_mean_i=Z_mean_i+Z_array{iter}{i};
    end
    Z_mean_i=Z_mean_i/size(Z_array,1);
    Z_mean{i}=Z_mean_i;
end
%% Map posterior connectivities of group components
C_BAR_array=cell(size(B_BAR_array));
for iter=1:size(B_BAR_array,1)
    C_BAR_array{iter}=cell(M);
    for k1=1:M
        for k2=1:M
            if k1==k2
                C_BAR_array{iter}{k1,k2}=phi_t'*Theta_diag*B_BAR_array{iter}{k1,k2}';
            end
            if not(k1==k2)
                C_BAR_array{iter}{k1,k2}=phi_t'*Theta_offdiag*B_BAR_array{iter}{k1,k2}';
            end
        end
    end
end


labels = {'BA19/39 (+/-36,-70,23)' , 'BA31 (12,-41,-34)' , 'BA19 (-32,-70,22)' , 'BA32 (-16,13,34)' , 'BA32/ACC (1,30,29)' , 'BA24/MCC (1,-2,35)' , 'BA32 (20,7,32)' , 'BA6/preCG (52,-1,14)' , 'BA2/pCG (-33,-28,43)' , 'BA20 (-39,-9,-19)' , 'BA7/Cuneus (3,-73,32)' , 'Parietal WM (27,-35,37)'};

%% SELECT SIGNIFICANT COMPONENTS
num_subjects = length(subj_keep);
num_clusters = size(S_bar_mean,1);
T=80;
BASELINE = [-0.7 -0.5];
TAIL = 'both';
erWinCenterTimes = linspace(-0.75, 1.7188, T);

p_cmp = 0.5;  % probably threshold for cluster membership
p_sub = 0.5; % proportion of subjects that must be assigned (with probability > p_cmp) to a cluster in order to keep the cluster
p_M=zeros(num_clusters,num_subjects);  % p_M(k,i) probability that subject i has membership in cluster k
for k=1:num_clusters
    for i=1:num_subjects
        p_M(k,i)=max(Z_mean{i}(:,k));
    end
end
keep=zeros(num_clusters,1);
for k=1:num_clusters
    pr_k=sum(p_M(k,:)>p_cmp)/num_subjects;
    if pr_k>p_sub
        keep(k)=1;
    end
end
S_bar_mean=S_bar_mean(find(keep==1),:);
num_clusters = size(S_bar_mean,1);


%% PUT MNI COORDS IN DIPFIT STRUCTURE

%% insert symmetric dipole on opposite y-hemisphere
% indices of dual equivalent dipoles
dualEquivDipoles = 0;

for k=1:num_clusters
    centroids(k).momxyz = eps*ones(2,3);
    centroids(k).rv = 0;
    centroids(k).clustid = k;
    
    if ismember_bc(k,dualEquivDipoles)
        centroids(k).posxyz = [S_bar_mean(k,:); S_bar_mean(k,:).*[1 -1 1]];
        centroids(k).select = [1 2];
    else
        centroids(k).posxyz = [S_bar_mean(k,:); [0 0 0]];
        centroids(k).select = 1;
    end
end


%% do a dipplot of the cluster centroids

dipplot(centroids,'coordformat','Spherical','plot','on','projlines','on');

%% Create connectivity structure

M_keep=sum(keep);
ind_keep=find(keep==1);
alpha = 0.01;
if strcmpi(TAIL,'both')
    alpha=alpha/2;
end
% t_bl=10;    % time points 1:t_bl form baseline
q_l=alpha; % quantile for lower CI bound
q_u=1-alpha; % quantile for upper CI bound

baseidx = getindex(erWinCenterTimes,BASELINE);

C_bar_bl_mean=cell(size(C_BAR_array,1),1);
for iter=1:size(C_BAR_array,1)
    C_bar_bl_mean{iter}=zeros(M_keep);
    for k1=1:M_keep
        for k2=1:M_keep
            C_bar_bl_mean{iter}(k1,k2)=mean(C_BAR_array{iter}{ind_keep(k1),ind_keep(k2)}(baseidx(1):baseidx(2)));
        end
    end
end

C_bar_lower_keep=cell(M_keep);
C_bar_mean_keep=cell(M_keep);
C_bar_upper_keep=cell(M_keep);
C_bar_sgn_keep=cell(M_keep); % indicates which time points are significantly different from baseline mean
C_bar_keep=cell(M_keep);
for k1=1:M_keep
    for k2=1:M_keep
        C_bar_lower_keep{k1,k2}=zeros(T,1);
        C_bar_mean_keep{k1,k2}=zeros(T,1);
        C_bar_upper_keep{k1,k2}=zeros(T,1);
        C_bar_sgn_keep{k1,k2}=zeros(T,1);
        C_bar_keep{k1,k2}=zeros(T,1);
        for t=1:T
            tmp=zeros(size(C_BAR_array,1),1);
            for iter=1:size(C_BAR_array,1)
                tmp(iter)=C_BAR_array{iter}{ind_keep(k1),ind_keep(k2)}(t);
            end
            C_bar_lower_keep{k1,k2}(t)=prctile(tmp,100*q_l);
            C_bar_mean_keep{k1,k2}(t)=mean(tmp);
            C_bar_upper_keep{k1,k2}(t)=prctile(tmp,100*q_u);
            if strcmpi(TAIL,'both')
                if or((C_bar_lower_keep{k1,k2}(t)-C_bar_bl_mean{iter}(k1,k2))>0,(C_bar_upper_keep{k1,k2}(t)-C_bar_bl_mean{iter}(k1,k2))<0)
                    C_bar_sgn_keep{k1,k2}(t)=1; % significant
                end
            else
                if (C_bar_lower_keep{k1,k2}(t)-C_bar_bl_mean{iter}(k1,k2))>0
                    C_bar_sgn_keep{k1,k2}(t)=1; % significant
                end
            end
            C_bar_keep{k1,k2}(t)=exp(C_bar_mean_keep{k1,k2}(t))-1;
        end
    end
end


for i=1:M_keep
    for j=1:M_keep
        GroupConn.dDTF08(j,i,1,:) = C_bar_keep{i,j};
    end
end

GroupConn.erWinCenterTimes = erWinCenterTimes;
GroupConn.freqs = 5;  % 3-7 Hz


%% Create SIFT EEG structure
EEG_orig = pop_loadset('C:/aaa_wes/brains/multisubject_eeg/For_Wes/Data/eb79RespWrong_Connectivity.set');

SamplingRate = EEG_orig.srate;

EEG = eeg_emptyset;
EEG.data = zeros(EEG_orig.nbchan,T);
EEG.icaweights = eye(num_clusters,EEG_orig.nbchan);
EEG.icasphere =  EEG.icaweights';
EEG.icachansind = 1:num_clusters;
EEG.icaact = [];
EEG.srate = SamplingRate;
EEG.pnts = T;
EEG.trials = 1;
EEG.chanlocs = EEG_orig.chanlocs;
EEG.times = ((0:(EEG.pnts-1))/SamplingRate)*1000;   % ms
EEG.xmin = EEG.times(1);
EEG.xmax = EEG.times(end)/1000;  % sec
EEG.nbchan = EEG_orig.nbchan;
EEG.dipfit = EEG_orig.dipfit;
EEG.dipfit.model = centroids;
EEG.dipfit.chansel = 1:num_clusters;
EEG = eeg_checkset(EEG);
% pop_eegplot(EEG);


EEG.CAT.curComps = 1:num_clusters;
EEG.CAT.curComponentNames = labels(find(keep)); %strtrim(cellstr(num2str(EEG.CAT.curComps')))';
EEG.CAT.Conn = GroupConn;
EEG.CAT.MODEL = struct([]);
EEG.CAT.nbchan = num_clusters;
EEG.CAT.srcdata = EEG.data;


%% Make statistics
for i=1:size(EEG.CAT.Conn.dDTF08,1)
    for j=1:size(EEG.CAT.Conn.dDTF08,2)
        Stats.dDTF08.ci(1,j,i,1,:) = exp(C_bar_lower_keep{i,j})-1;  % lower bound
        Stats.dDTF08.ci(2,j,i,1,:) = exp(C_bar_upper_keep{i,j})-1;  % upper bound
        Stats.dDTF08.pval(j,i,1,:) = ~C_bar_sgn_keep{i,j}; %squeeze(Stats.dDTF08.ci(1,j,i,1,:))<=0 & squeeze(Stats.dDTF08.ci(2,j,i,1,:))>=0;  % insignificant if zero is contained within lower and upper bound
        Stats.alpha = 0.05;
    end
end

%% 
pop_dipplot(EEG,1:num_clusters,'dipolelength',0,'spheres','on','projlines','on','coordformat','Spherical','num','on','dipnames',EEG.CAT.curComponentNames);

%% Visualize
[figureHandles tfgridcfg] = vis_TimeFreqGrid('ALLEEG',EEG,'Conn',EEG.CAT.Conn,'Stats',Stats,'MatrixLayout',{'Partial','UpperTriangle', 'dDTF08', 'LowerTriangle','dDTF08','Diagonal','none'},'TimesToPlot',[-0.7 1],'ColorLimits',[0 0.05],'Baseline',[],'Smooth2D',false,'Thresholding',{'Statistics','plotci',true,'ThresholdingMethod','pval','AlphaSignificance',0.05},'BackgroundColor',[0 0 0],'TextColor',[1 1 1],'LineColor',[0 0 0],'AxesFontSize',14);

%% threshold according to stats
load brainmovieconfig
EEGthresh = EEG;
EEGthresh.CAT.Conn.dDTF08(Stats.dDTF08.pval>=alpha) = 0;

cfg=pop_vis_causalBrainMovie3D(EEGthresh,cfg);
