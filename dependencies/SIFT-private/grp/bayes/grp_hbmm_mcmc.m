% hierarchical bayesian mixture model
function [MCMC_State] = grp_hbmm_mcmc(varargin)

% inputs:
%   B: {subj}(M_i x M_i x Q) matrix of (possibly time-varying)  (obtain from EEG(i).CAT.Conn)
%      connectivities or basis coefficients
%   S: {subj}(M_i x 3) matrix of dipole locations (obtain from EEG(i).dipfit)
%   M_i: [N x 1] vector of number of sources for each subject
%   M:  number of clusters required
% automatically determined:
%   N:  number of subjects (length(B))
%   Q:  number of basis coefficients (size(B{1},3))
%   k:  clustering mode (1 = kmeans, 2=gmm)

arg_define([0 3], varargin, ...
    arg_norep({'B','Connectivity'},mandatory,[],'Connectivity matrices. B{i} is an (M_i x M_i x Q) matrix of (possibly time-varying) connectivity values or basis coefficients for the ith subject. M_i is the number of components for the ith subject.'), ...
    arg_norep({'S','DipoleLocations'},mandatory,[],'Dipole locations. S{i} is an (M_i x 3) matrix of [X Y Z] dipole locations for the ith subject. M_i is the number of components for the ith subject.'), ...
    arg_norep({'MCMC_InitState'},struct([]),[],sprintf('Initial state of MCMC sampler. This can be computed from grp_hbmm_initMCMC or can be the output of a previous call to grp_hmbb_mcmc (i.e. the last state of the MCMC iterator). This structure must contain the following fields:\nZ\t: Mi x M matrices of group indicators\nS_BAR\t: cluster centroid locations\nSIGMA_S\t: cluster centroid variances\nB_BAR\t: group level connectivities\nSIGMA_B\t: variances of connectivities\nN_k\t: number of components that belong to each cluster\nN_k1k2\t: pairwise counts of group level clusters'),'type','struct'), ...
    arg_sub({'hyperparams','Hyperparameters'},[],...
    {...
        arg({'c','ConnMeanPriorVar'},10000,[eps Inf],'Between-cluster mean connectivity prior variance. Variance of gaussian prior for between-cluster mean connectivity. Larger --> more uncertainty'), ...
        arg({'eta','DipoleLocVarPriorShape'},[],[],'Within-cluster dipole location prior variance D.O.F. Degrees of freedom for inverse-wishart prior distribution for within-cluster dipole location covariance matrices. Leave empty to compute from initial clustering.'), ...
        arg({'SS','DipoleLocVarPriorScale'},[],[],'Within-cluster dipole location prior variance scale. Scale matrix of inverse-wishart prior distribution for within-cluster dipole location covariance matrices. Leave empty to compute from initial clustering.'), ...
        arg({'p1','ConnVarPriorShape','a'},1,[],'Between-cluster mean connectivity variance shape.'), ...
        arg({'p2','ConnVarPriorScale','b'},1,[],'Between-cluster mean connectivity variance scale.'), ...
    },'MCMC hyperparameters'), ...
    arg({'nMCMCiters','NumMCMCIters','niter','niters'},1000, [1 Inf], 'Number of MCMC iterations for spline fitting'), ...
    arg({'burnInFraction','BurnInFractionForMCMC'},0.5,[0 0.99],'Fraction of initial MCMC samples to discard (burn in period). The number of MCMC samples is taken to be the smaller of MCMCitersOffDiag and MCMCitersDiag.'), ...
    arg({'thinFactor','ThinningFactor'},1,[1 Inf],'Thinning factor for MCMC. We keep every kth MCMC sample, where k=ThinningFactor. This is useful when we have limited available memory to store MCMC results since successive MCMC estimates are more likely to be correlated.'), ...
    arg({'normlog','NormalizeAndLogTransform'},true,[],'Transform data before smoothinFactorg. Normalize across last dim (time), add 1 and take logarithm. Inverse transform is applied after smoothinFactorg'), ...
    arg({'verb','VerbosityLevel'},2,{int32(0) int32(1) int32(2)},'Verbosity level. 0 = no output, 1 = text, 2 = graphical'), ...
    arg({'appendLastState','AppendLastState'},false,[],'Append new state to initial MCMC state') ... 
    );

%     arg({'basisCoeffVarPrior'},1000,[eps Inf],'Variance of basis coefficient gaussian prior. Larger --> more wiggling allowed'), ...
%     arg({'noiseVarPriorShape'},0.01,[eps Inf],'Shape (D.O.F) of noise variance prior. This is the "alpha" parameter of the inverse gamma prior distribution. Increasing noiseVarPriorShape --> decreased variance of noise variance distribution.'), ...
%     arg({'noiseVarPriorScale'},0.01,[eps Inf],'Scale parameter of noise variance prior. This is the "theta" (1/beta) parameter of inverse gamma prior distribution. Increasing noiseVarPriorScale --> right-shift of distribution --> (increase in expected noise variance). In general MEAN(noiseVariance) = noiseVarPriorScale/noiseVarPriorShape and MODE(noiseVariance) = noiseVarPriorScale/(noiseVarPriorShape-1) for noiseVarPriorShape>=1.'), ...
%     arg({'initNoiseVariance'},0.1,[eps Inf],'Initial noise variance'), ...
    

% set up initial state of MCMC
if isempty(MCMC_InitState)
    error('You must provide an initial state for the MCMC iterator (MCMC_InitState)');
end

% initialize vars based on last state of MCMC_InitState
% this will initialize the following variables:
% 'Z','S_BAR','SIGMA_S','B_BAR','SIGMA_B','N_k','N_k1k2'
varnames = setdiff_bc(fieldnames(MCMC_InitState),'initstate');
if MCMC_InitState.initstate
    % initialize to initial state of MCMC iterator
    for i = 1:length(varnames)
        eval(sprintf('%s=MCMC_InitState.(''%s'');',varnames{i},varnames{i})); 
    end
else
    % initialize to last state of MCMC iterator
    for i = 1:length(varnames)
        eval(sprintf('%s=MCMC_InitState.(''%s''){end};',varnames{i},varnames{i})); 
    end
end
if ~appendLastState
    clear MCMC_InitState; 
end

% define some vars
M   = length(N_k);                      % number of clusters
N   = length(B);                        % number of subjects
M_i = cellfun(@(B_i) size(B_i,1),B);    % number of sources for each subject
Q   = size(B{1},3);                     % connectivity time-series dimension

% initialize prior probabilities of cluster membership
if ~exist('MU','var')
    MU  = ones(M,1)/M; end
% initialize hyperparams
if isempty(hyperparams.eta)
    hyperparams.eta=median(N_k); end
if isempty(hyperparams.SS)
    hyperparams.SS = hyperparams.eta*mean(SIGMA_S,3); end

% compute number of burn-in samples
numBurnInSamples = floor(burnInFraction*nMCMCiters);
niterToKeep      = round((nMCMCiters-numBurnInSamples)/thinFactor);
if verb,
    fprintf(['I will discard %d burn-in samples.\n' ...
             'I will thin the distribution by a factor of %d samples\n', ...
             'The distribution of the estimator will have %d samples\n'], ...
            numBurnInSamples,thinFactor,niterToKeep);
end

%% initialize arrays of posterior draws
Z_array         = cell(niterToKeep,1);
S_BAR_array     = cell(niterToKeep,1);
SIGMA_S_array   = cell(niterToKeep,1);
B_BAR_array     = cell(niterToKeep,1);
SIGMA_B_array   = cell(niterToKeep,1);
MU_array        = cell(niterToKeep,1);
N_k_array       = cell(niterToKeep,1);
N_k1k2_array    = cell(niterToKeep,1);
iter_array      = 0;


%% run MCMC 
for iter=1:nMCMCiters

    % Draw S_BAR (group centroid locations)
    SSinv = double(inverse(hyperparams.SS));
    for k=1:M
        mu_s_bar_k=zeros(3,1);
        Sigma_s_bar_k=double(inverse(N_k(k)*double(inverse(SIGMA_S(:,:,k)))+SSinv));
        for i=1:N
            if max(Z{i}(:,k))==1
                j_k=find(Z{i}(:,k)==1);
                for j=1:length(j_k)
                    mu_s_bar_k=mu_s_bar_k+S{i}(j_k(j),:)';
                end
            end
        end
        mu_s_bar_k=Sigma_s_bar_k/SIGMA_S(:,:,k)*mu_s_bar_k;
        S_BAR(k,:)=mvnrnd(mu_s_bar_k,Sigma_s_bar_k);
    end

    
    % Draw SIGMA_S (group centroid covariance matrices)
    for k=1:M
        eta_k=hyperparams.eta+.5*N_k(k);
        SS_k=hyperparams.SS;
        for i=1:N
            if max(Z{i}(:,k))==1
                j_k=find(Z{i}(:,k)==1);
                for j=1:length(j_k)
                    SS_k=SS_k+.5*(S{i}(j_k(j),:)-S_BAR(k,:))'*(S{i}(j_k(j),:)-S_BAR(k,:));
                end
            end
        end
        SS_k=.5*(SS_k+SS_k');
        inv_SS_k=double(inverse(SS_k));
        inv_SS_k=.5*(inv_SS_k+inv_SS_k');
        SIGMA_S(:,:,k)=double(inverse(wishrnd(inv_SS_k,eta_k)));
    end  
    
    
    % Draw B_BAR (group mean connectivites)
    for k1=1:M
        for k2=1:M
            mu_b_bar_k1k2=zeros(Q,1);
            Sigma_sq_b_k1k2=1/(1/hyperparams.c+N_k1k2(k1,k2)/SIGMA_B(k1,k2))*eye(Q);
            for i=1:N
                n_k1=sum(Z{i}(:,k1));
                n_k2=sum(Z{i}(:,k2));
                if(and(n_k1>=1,n_k2>=1))
                    j_k1=find(Z{i}(:,k1)==1);
                    j_k2=find(Z{i}(:,k2)==1);
                    for j1=1:n_k1
                        for j2=1:n_k2
                            mu_b_bar_k1k2=mu_b_bar_k1k2+squish(B{i}(j_k1(j1),j_k2(j2),:));
                        end
                    end
                end
            end
            mu_b_bar_k1k2=Sigma_sq_b_k1k2*mu_b_bar_k1k2/SIGMA_B(k1,k2);
            B_BAR(k1,k2,:)=mvnrnd(mu_b_bar_k1k2,Sigma_sq_b_k1k2);
        end
    end
    
    % Draw SIGMA_B (group connectivity variances)
    for k1=1:M
        for k2=1:M
            p1_k1k2=hyperparams.p1+.5*Q*N_k1k2(k1,k2);
            p2_k1k2=hyperparams.p2;
            for i=1:N
                n_k1=sum(Z{i}(:,k1));
                n_k2=sum(Z{i}(:,k2));
                if(and(n_k1>=1,n_k2>=1))
                    j_k1=find(Z{i}(:,k1)==1);
                    j_k2=find(Z{i}(:,k2)==1);
                    for j1=1:n_k1
                        for j2=1:n_k2
                            p2_k1k2=p2_k1k2+.5*sum((B{i}(j_k1(j1),j_k2(j2),:)-B_BAR(k1,k2,:)).^2);
                        end
                    end
                end
            end
            SIGMA_B(k1,k2)=1/gamrnd(p1_k1k2,1/p2_k1k2);
        end
    end


    % Draw Z (indicators of group membership)   
    ind=randperm(N);
    for i=ind
        ind_i=randperm(M_i(i));
        for j=ind_i
            Log_p_ij=zeros(M,1);
            for k=1:M
                Sigma_s_ijk=SIGMA_S(:,:,k);
                s_ijk=S{i}(j,:)-S_BAR(k,:);
                log_p_s_ijk=-.5*log(det(Sigma_s_ijk))-.5*s_ijk/Sigma_s_ijk*s_ijk';
                log_p_b_ijk=0;
                for j1=1:M_i(i)
                    if j==j1
                        log_p_b_ijk=log_p_b_ijk-.5*Q*log(det(SIGMA_B(k,k)))...
                        -.5*sum((B{i}(j,j,:)-B_BAR(k,k,:)).^2)/SIGMA_B(k,k);
                    end
                    if and(not(j==j1),find(Z{i}(j1,:)==1)==k)
                        log_p_b_ijk=log_p_b_ijk-.5*Q*log(det(SIGMA_B(k,k)))...
                        -.5*sum((B{i}(j1,j1,:)-B_BAR(Z{i}(j1,:)==1,Z{i}(j1,:)==1,:)).^2)/...
                            SIGMA_B(Z{i}(j,:)==1,Z{i}(j,:)==1);                    
                    end
                    if not(find(Z{i}(j,:)==1)==find(Z{i}(j1,:)==1))
                        log_p_b_ijk=log_p_b_ijk-.5*Q*log(det(SIGMA_B(k,Z{i}(j1,:)==1)))...
                        -.5*sum((B{i}(j,j1,:)-B_BAR(k,Z{i}(j1,:)==1,:)).^2)/...
                              SIGMA_B(k,Z{i}(j1,:)==1);
                    end
                end
                for j2=1:M_i(i)
                    if not(find(Z{i}(j2,:)==1)==find(Z{i}(j,:)==1))
                        log_p_b_ijk=log_p_b_ijk-.5*Q*log(det(SIGMA_B(Z{i}(j2,:)==1,k)))...
                        -.5*sum((B{i}(j2,j,:)-B_BAR(Z{i}(j2,:)==1,k,:)).^2)/...
                                SIGMA_B(Z{i}(j2,:)==1,k);
                    end
                end
                Log_p_ij(k)=MU(k)+log_p_s_ijk+log_p_b_ijk;
            end
            Log_p_ij=Log_p_ij-max(Log_p_ij);
            P_ij=exp(Log_p_ij)/sum(exp(Log_p_ij));
            v=rand;
            Z{i}(j,:)=0*Z{i}(j,:);
            if v<=P_ij(1)
                   Z{i}(j,1)=1;
                   log_p_ij=Log_p_ij(1);
            end
            for k=1:(M-1)
               if(and(v>sum(P_ij(1:k)),v<=sum(P_ij(1:(k+1)))))
                   Z{i}(j,k+1)=1;
                   log_p_ij=Log_p_ij(k);
               end
            end
        end
    end

    % Draw MU (probabilities of clusters)   
    N_k1k2=zeros(M);
    for k1=1:M
        for k2=1:M
            for i=1:N
                n_k1=sum(Z{i}(:,k1));
                n_k2=sum(Z{i}(:,k2));
                if(and(n_k1>=1,n_k2>=1))
                    if not(k1==k2)
                        N_k1k2(k1,k2)=N_k1k2(k1,k2)+n_k1*n_k2;
                    end                
                    if k1==k2
                        N_k1k2(k1,k2)=N_k1k2(k1,k2)+n_k1;
                    end
                end
            end
        end
    end
    N_k=diag(N_k1k2);
    MU=drchrnd(1/M*ones(1,M)+N_k',1);

    
    % Posterior Draw arrays    
    if and(iter>=numBurnInSamples,thinFactor*round(iter/thinFactor)==iter)
        iter_array=iter_array+1;
        Z_array{iter_array}=Z;
        S_BAR_array{iter_array}=S_BAR;
        SIGMA_S_array{iter_array}=SIGMA_S;
        B_BAR_array{iter_array}=B_BAR;
        SIGMA_B_array{iter_array}=SIGMA_B;
        N_k_array{iter_array}=N_k;
        N_k1k2_array{iter_array}= N_k1k2;
        MU_array{iter_array}=MU;

        
        if verb
            disp(['iter: ' num2str(iter)]);
%             % compute and print out intermediate estimates of source locations
%             S_BAR_mean=S_BAR_array{1};
%             for j=2:iter_array
%                 S_BAR_mean=S_BAR_mean+S_BAR_array{j};
%             end
%             S_BAR_mean=S_BAR_mean/iter_array;
%             disp([S_BAR zeros(M,1) S_BAR_mean])
        end
        
        %[sqrt(sum(min(pdist2(S_BAR_true,S_BAR_st)').^2))/M sqrt(sum(min(pdist2(S_BAR_true,S_BAR_mean)').^2))/M]
      
        
        if verb
%             %compute and print intermediate indicator probs for random subject
%             i=discretesample(ones(N,1)/N,1);
%             Z_i_mean=Z_array{1}{i};
%             for j=2:iter_array
%                Z_i_mean=Z_i_mean+Z_array{j}{i};
%             end
%             Z_i_mean=Z_i_mean/iter_array;
% %         Z_st{i}
%            disp(Z_i_mean);
        end
    end
end


% return MCMC state
varnames = {'Z_array','S_BAR_array','SIGMA_S_array','B_BAR_array','SIGMA_B_array','MU_array','N_k1k2_array','N_k_array'};
for i = 1:length(varnames)
    % state fieldnames have trailing '_array' stripped
    fname = strrep(varnames{i},'_array','');
    vname = varnames{i};
    if appendLastState
        % append current state estimates to last state estimate
        MCMC_State.(fname) = [MCMC_State.(fname); eval(vname)];
    else
        % replace last state estimate with current state estimates 
        MCMC_State.(fname) = eval(vname); 
    end
end
MCMC_State.initstate = false;


function r = drchrnd(a,n)
    p = length(a);
    r = gamrnd(repmat(a,n,1),1,n,p);
    r = r ./ repmat(sum(r,2),1,p);

    