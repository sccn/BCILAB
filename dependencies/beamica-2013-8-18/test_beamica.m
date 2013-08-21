% To run this you need:
% * EEGLAB
% * Guido Nolte's sourceanalysis toolbox on the path (included here)
addpath([pwd filesep 'nolte']);

% load a high-pass filtered data et
EEG = pop_loadset('demo_BF.set');
X = double(EEG.data);

% pick locations (see calc_beamformer_constraints.m for the options)
locs = {'Frontal_Sup_L', 'Frontal_Sup_R', 'Frontal_Mid_L', 'Frontal_Mid_R', 'Precentral_L', 'Precentral_R', 'Supp_Motor_Area_L', 'Supp_Motor_Area_R', 'Frontal_Sup_Medial_L', 'Frontal_Sup_Medial_R', 'Frontal_Med_Orb_L', 'Frontal_Med_Orb_R', 'Cingulum_Ant_L', 'Cingulum_Ant_R', 'Postcentral_L', 'Postcentral_R', 'Parietal_Sup_L', 'Parietal_Sup_R', 'Parietal_Inf_L', 'Parietal_Inf_R', 'Precuneus_L', 'Precuneus_R', 'Temporal_Sup_L', 'Temporal_Sup_R', 'Temporal_Mid_L', 'Temporal_Mid_R', 'Occipital_Sup_L', 'Occipital_Sup_R', 'Occipital_Mid_L', 'Occipital_Mid_R', 'Occipital_Inf_L', 'Occipital_Inf_R'};

% calc robust mean/cov
mu = median(X,2);
sig = cov_blockgeom(X');
S = inv(sqrtm(sig));

%% calc constraint set & show cortex-perpendicular projections

% calc beamformer constraint matrices (B)
[B,Winit,Chanmask,anchor_pos,~,sa] = calc_beamformer_constraints({EEG.chanlocs.labels},locs,sig,'nasion');


% plot cortex-perpendicular forward projections at anchor locations
% (note: Winit is actually not used in the subsequent ICA solutions)
figure;topoplot_grid(sig*(Winit*S)',EEG.chanlocs,'titles',locs);

%% run ICA solutions

% run 1000 iterations of vanilla Infomax, initialized with eye(32), empty constraint set
[W0,S0] = beamica(X,{},eye(32),inv(sqrtm(sig)),mu,1000,0.15);
% [W1,S1] = runica(X); -- also try runica() for comparison

% run 750 iterations of Beam-constrained Infomax, slightly lower learning rate, beam cost 5x that of Infomax cost
[W,S] = beamica(X,B,eye(32),inv(sqrtm(sig)),mu,750,0.05,0.8);

%% plot mixing matrix for vanilla ICA
figure;topoplot_grid(inv(W0*S0),EEG.chanlocs)
% plot it without assuming that W is invertible (for comparison with the BF case, which is not necessarily invertible)
% this is using Haufe's trick (should look virtually identical)
figure;topoplot_grid(sig*(W0*S0)',EEG.chanlocs)
% plot forward projection of the Beam case
figure;topoplot_grid(sig*(W*S)',EEG.chanlocs,'titles',locs)

%% also plot spatial filters designed by vanilla Beamforming (for cortex-perpendicular sources)
% versus Beam-Infomax designed filters
figure;topoplot_grid((Winit*S)',EEG.chanlocs,'titles',locs);
figure;topoplot_grid((W*S)',EEG.chanlocs,'titles',locs);


%% fit dipoles to solution maps and measure distance to anchors...
icawinv = sig*(W*S)';
numsamples = 15;
offsets = [];
dipoles = {};
residuals = [];
for k=1:size(icawinv,2)
    errors = [];
    fields = [];
    dips = {};
    fprintf('fitting dipole');
    for r=1:numsamples
        [dips{end+1},errors(end+1),x,fields(:,end+1)] = dipole_fit_field(icawinv(:,k),sa.fp,1); % one iteration takes ca. 1 second
        fprintf('.');
    end
    fprintf('\n');
    [besterr,bestidx] = min(errors);
    offsets(k) = sqrt(sum((anchor_pos(k,:)-dips{bestidx}(1:3)).^2));
    dipoles{k} = dips{bestidx};
    residuals(k) = errors(bestidx);
    fprintf('Offset is  %.1f cm (residual error = %.1f%%).\n',offsets(k),residuals(k)*100);
    % figure;topoplot_grid([fields,fields(:,bestidx),icawinv(:,k)],ica.state.root_chanlocs)
end

% calculate residual variances in dipfit style...
rv = residuals.^2;

