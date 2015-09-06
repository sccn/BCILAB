function [Conn PConn MCMC_LastState] = grp_stat_hbmm(varargin)

% if hlp_isToolboxInstalled('Parallel Computing Toolbox')
%     pardef   = 'on';
%     [tmp parprofs] = hlp_microcache('sift_domain',@defaultParallelConfig);
% else
%     pardef = 'off';
%     parprofs = {'local'};
% end

% extract some stuff from inputs for arg defaults
EEG = arg_extract(varargin,'EEG',1);

if ~isempty(EEG)
    Conn            = EEG(1).CAT.Conn;
    ConnNames       = hlp_getConnMethodNames(Conn);
    freqRangeDef    = [Conn.freqs(1) Conn.freqs(end)];
    timeRangeDef    = [Conn.erWinCenterTimes(1) Conn.erWinCenterTimes(end)];
    clear Conn;
else
    ConnNames = {''};
    [freqRangeDef, timeRangeDef] = deal([]);
end

g = arg_define([0 Inf], varargin, ...
    arg_norep({'ALLEEG','EEG'},mandatory,[],'Array of EEG datasets'), ...
    arg_nogui({'MCMC_LastState'},struct([]),[],'Structure containing state of MCMC sampler. This will be used to re-initialize the sampler'), ...
    arg({'connmethods','Estimator'},ConnNames,ConnNames,'Connectivity estimators to smooth','type','logical'), ...
    arg({'timeRange','TimeRange'},timeRangeDef,[],'[Min Max] Time range to smooth (sec). Leave blank to use all time points','shape','row','type','denserealdouble'), ...
    arg({'freqRange','FrequencyRange'},freqRangeDef,[],'[Min Max] Frequency range to smooth (Hz). Leave blank to use all frequencies','type','expression','shape','row'), ...
    arg({'collapseTime','CollapseTime'},false,[],'Average across time before smoothing'), ...
    arg({'collapseFreqs','CollapseFreqs'},true,[],'Integrate across frequencies before smoothing'), ...
    arg_sub({'smoothData','SmoothData'},{},@stat_ubsplsm,'Apply B-spline smoothing to data before group inference','cat','Smoothing'), ... 
    arg_sub({'mcmc_opts','MCMC'},{}, ...
    {...
        arg_sub({'mcmc_init','InitMCMC'},{},@grp_hbmm_init,'Initialize MCMC','suppress','verb') ...
        arg_sub({'mcmc_run','RunMCMC'},{},@grp_hbmm_mcmc,'MCMC runtime options','suppress','verb') ...
    },'Perform Markov Chain Monte Carlo inference'), ...
    arg({'verb','VerbosityLevel'},2,{int32(0) int32(1) int32(2)},'Verbosity level. 0 = no output, 1 = text, 2 = graphical') ...
    );

% arg_subtoggle({'runPll','SmoothInParallel','RunInParallel'},'off', ...
%     { ...
%         arg({'profile','ProfileName'},parprofs{1},parprofs,'Profile name'), ...
%         arg({'numWorkers','NumWorkers'},2,[1 Inf],'Number of workers') ...
%     },'Smooth connectivities in parallel. Requires Parallel Computing Toolbox.','cat','Smoothing'), ...   
% 
% if strcmp(pardef,'off') && g.runPll.arg_selection
%     fprintf('Parallel Computing Toolbox not installed. Cannot use parallel option.\n');
%     g.runPll.arg_selection = false;
% end

if ~g.collapseTime && ~g.collapseFreqs
    error('You must collapse at least one dimension (times,freqs)');
end

nsubj = length(EEG);

for si=1:nsubj
    
    Conn = EEG(si).Conn;
    
    % collapse connectivity matrices
    % -------------------------------------------------------------------------
    if g.collapseFreqs
        % collapse across freqs
        Conn = hlp_filterConns(Conn,'connmethods',g.connmethods, ...
            'method',{'freq','net'},'frange',g.freqRange,'freqdim',3,'timedim',4,'verb',g.verb);
    else
        % only select subset of freqs
        Conn = hlp_filterConns(Conn,'connmethods',g.connmethods, ...
            'method',{'freq','shrinkonly'}, ...
            'frange',g.freqRange, 'freqdim',3,'timedim',4,'verb',g.verb);
    end

    if g.collapseTime
        % collapse across time
        Conn = hlp_filterConns(Conn,'connmethods',g.connmethods, ...
            'method',{'time','mean'},'trange',g.timeRange,'freqdim',3,'timedim',4,'verb',g.verb);
    else
        % only select subset of times
        Conn = hlp_filterConns(Conn,'connmethods',g.connmethods, ...
            'method',{'time','shrinkonly'}, ...
            'trange',g.timeRange,'freqdim',3,'timedim',4,'verb',g.verb);
    end
    
    EEG(si).Conn = Conn;
end


for m=1:length(g.connmethods)
    % smooth the time-series
    if g.smoothData.arg_selection
        connmethod = g.connmethods{m};
        if g.verb
            waitbarstr = sprintf('Applying MCMC smoothing to %s',connmethod);
            multiWaitbar(waitbarstr,0,'Color',[0.8 0.0 0.1]);   
        end
        % construct YDATA cell array
        % YDATA{i} contains the [M_i x M_i x Q] connectivity matrix for
        % subject i
        YDATA = cellfun(@(EEG_i) squeeze(EEG_i.CAT.Conn.(connmethod)), ...
                          EEG,'UniformOutput',false);
        % univariate (1D) smoothing
        [fits MCMC_LastState] = stat_ubsplsm(YDATA,g.smoothData);
        
        % store the coefficients MCMC_LastState.ALPHA (?) in DATA
        
        % take mean of posterior for ALPHA coefficients and rotate to orthonormal basis
        
        
        if g.verb
            % cleanup waitbar
            multiWaitbar(waitbarstr, 'Close');
        end
    else
        % extract the connectivity data for all subjects and store in DATA
        
    end
    
    % perform MCMC inference
    
    % store result for this connmethod
end
    
    



function rot_orthn(MCMC_LastState)

% rotate and scale the basis functions for diagonal connectivities so that the resulting
% smoothing coefficents are orthonormal. This is because the mean 
% posterior estimates from smooth_mcmc are not necessarily orthonormal
Alpha_diag_mean=mean(Alpha_diag,2);  % mean of posterior estimates of smoothing coefficients
S_Alpha_diag=(Alpha_diag-repmat(Alpha_diag_mean,1,size(Alpha_diag,2)))*...
    (Alpha_diag-repmat(Alpha_diag_mean,1,size(Alpha_diag,2)))'; % sample covariance of coefs
Sigma_psi_diag=Theta_diag*S_Alpha_diag*Theta_diag';  
Sigma_psi_diag=.5*(Sigma_psi_diag+Sigma_psi_diag');   
[U D]=eig(Sigma_psi_diag);
Theta_diag_new=U(:,(K-Q+1):K)*D((K-Q+1):K,(K-Q+1):K)^.5;
Alpha_diag_new=regress(squeeze(C{1}(1,1,:)),phi_t'*Theta_diag_new);
for i=1:N
    for j=1:M_i(i)
        Alpha_diag_new=[Alpha_diag_new...
            regress(squeeze(C{i}(j,j,:)),phi_t'*Theta_diag_new)];
    end
end
Alpha_diag_new=Alpha_diag_new(:,2:size(Alpha_diag_new,2));
Theta_diag=Theta_diag_new*cov(Alpha_diag_new')^.5;
Alpha_diag=cov(Alpha_diag_new')^-.5*Alpha_diag_new;
Fit_diag=phi_t'*Theta_diag*Alpha_diag;

% do same thing for off-diagonal connectivities
Alpha_offdiag_mean=mean(Alpha_offdiag,2);
S_Alpha_offdiag=(Alpha_offdiag-repmat(Alpha_offdiag_mean,1,size(Alpha_offdiag,2)))*...
    (Alpha_offdiag-repmat(Alpha_offdiag_mean,1,size(Alpha_offdiag,2)))';
Sigma_psi_offdiag=Theta_offdiag*S_Alpha_offdiag*Theta_offdiag';
Sigma_psi_offdiag=.5*(Sigma_psi_offdiag+Sigma_psi_offdiag');
[U D]=eig(Sigma_psi_offdiag);
Theta_offdiag_new=U(:,(K-Q+1):K)*D((K-Q+1):K,(K-Q+1):K)^.5;
Alpha_offdiag_new=regress(squeeze(C{1}(1,2,:)),phi_t'*Theta_offdiag_new);
for i=1:N
    for j1=1:M_i(i)
        for j2=1:M_i(i)
            if j1~=j2
                Alpha_offdiag_new=[Alpha_offdiag_new...
                    regress(squeeze(C{i}(j1,j2,:)),phi_t'*Theta_offdiag_new)];
            end
        end
    end
end
Alpha_offdiag_new=Alpha_offdiag_new(:,2:size(Alpha_offdiag_new,2));
Theta_offdiag=Theta_offdiag_new*cov(Alpha_offdiag_new')^.5;
Alpha_offdiag=cov(Alpha_offdiag_new')^-.5*Alpha_offdiag_new;
Fit_offdiag=phi_t'*Theta_offdiag*Alpha_offdiag;