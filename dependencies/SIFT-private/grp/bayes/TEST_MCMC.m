%% convert cell array B to matrix
Bc = B;
B  = cell(1,length(Bc));
for i=1:length(Bc)
    for j=1:size(Bc{i},1)
        for k=1:size(Bc{i},2)
            B{i}(j,k,:) = Bc{i}{j,k};
        end
    end
end

%% perform mcmc
cfg = arg_guipanel('Function',@grp_stat_hbmm,'Parameters',{'EEG',{}},'PanelOnly',false);

%% compute initial state (clustering)
cfg = arg_guipanel('Function',@grp_hbmm_init,'Parameters',{},'PanelOnly',false);
MCMC_InitState = grp_hbmm_init(B,S,cfg)

%% run MCMC
cfg = arg_guipanel('Function',@grp_hbmm_mcmc,'Parameters',{},'PanelOnly',false);
MCMC_LastState = grp_hbmm_mcmc(B,S,MCMC_InitState,cfg)