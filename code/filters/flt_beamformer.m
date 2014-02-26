function signal = flt_beamformer(varargin)
% Recovers activity from the given ROIs via a beamformer.
% function Signal = flt_laplace(Signal,ROIs,)
%
% In:
%   Signal : EEGLAB data set, either continuous or epoched
%
%   ROIs   : Regions of interest from which to recover signals.
%
%   FieldType : type of field to recover ('axial' or 'normal')
%
%   ReferenceType : type of referencing used before this filter ('nasion' or 'common_average')
%
%   OverrideOriginal : whether to override the original signal
%
% Out:
%   Signal : activity from the given ROIs
%
% Notes:
%   This function currently does not perform activity renormalization.
%
% Examples:
%   % recover activity from four regions
%   eeg = flt_beamformer(eeg,{'Precentral_L', 'Precentral_R', 'Frontal_Sup_L', 'Frontal_Sup_R'},'normal')
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2014-02-11

if ~exp_beginfun('filter') return; end

declare_properties('name','Beamformer', 'follows',{'flt_selchans','flt_repair_bursts'}, 'cannot_follow','flt_ica', 'independent_channels',false, 'independent_trials',true);

arg_define(varargin, ...
    arg_norep({'signal','Signal'}), ...
    arg({'roi_labels','ROIs'},{}, {'Precentral_L', 'Precentral_R', 'Frontal_Sup_L', 'Frontal_Sup_R', 'Frontal_Sup_Orb_L', 'Frontal_Sup_Orb_R', 'Frontal_Mid_L', 'Frontal_Mid_R', 'Frontal_Mid_Orb_L', 'Frontal_Mid_Orb_R', 'Frontal_Inf_Oper_L', 'Frontal_Inf_Oper_R', 'Frontal_Inf_Tri_L', 'Frontal_Inf_Tri_R', 'Frontal_Inf_Orb_L', 'Frontal_Inf_Orb_R', 'Rolandic_Oper_L', 'Rolandic_Oper_R', 'Supp_Motor_Area_L', 'Supp_Motor_Area_R', 'Frontal_Sup_Medial_L', 'Frontal_Sup_Medial_R', 'Frontal_Med_Orb_L', 'Frontal_Med_Orb_R', 'Insula_L', 'Insula_R', 'Cingulum_Ant_L', 'Cingulum_Ant_R', 'Cingulum_Mid_L', 'Cingulum_Mid_R', 'Cingulum_Post_L', 'Cingulum_Post_R', 'Hippocampus_L', 'Hippocampus_R', 'ParaHippocampal_L', 'ParaHippocampal_R', 'Calcarine_L', 'Calcarine_R', 'Cuneus_L', 'Cuneus_R', 'Lingual_L', 'Lingual_R', 'Occipital_Sup_L', 'Occipital_Sup_R', 'Occipital_Mid_L', 'Occipital_Mid_R', 'Occipital_Inf_L', 'Occipital_Inf_R', 'Fusiform_L', 'Fusiform_R', 'Postcentral_L', 'Postcentral_R', 'Parietal_Sup_L', 'Parietal_Sup_R', 'Parietal_Inf_L', 'Parietal_Inf_R', 'SupraMarginal_L', 'SupraMarginal_R', 'Angular_L', 'Angular_R', 'Precuneus_L', 'Precuneus_R', 'Paracentral_Lobule_L', 'Paracentral_Lobule_R', 'Temporal_Sup_L', 'Temporal_Sup_R', 'Temporal_Pole_Sup_L', 'Temporal_Pole_Sup_R', 'Temporal_Mid_L', 'Temporal_Mid_R', 'Temporal_Pole_Mid_L', 'Temporal_Pole_Mid_R', 'Temporal_Inf_L', 'Temporal_Inf_R', 'Olfactory_L', 'Olfactory_R', 'Rectus_L', 'Rectus_R', 'Amygdala_L', 'Amygdala_R', 'Caudate_L', 'Caudate_R', 'Thalamus_L', 'Thalamus_R', 'Heschl_L', 'Heschl_R'}, ...
        'Cortical anchor locations. List of locations to which components shall be constrained. The first k components are encouraged to lie close to the given locations, in the order of appearance. This is experimental and currently requires a) 10-20 locations and b) Guido Nolte''s source analysis toolbox (not included).','experimental',true), ...
    arg({'field_type','FieldType'},'normal',{'normal','axial'},'Regions of interest. These are the regions from which to recover signals.'), ...
    arg({'reference_type','ReferenceType'},'nasion',{'nasion','common_average'},'Referencing scheme. This is the type of re-referencing that was applied before flt_beamformer.'), ...
    arg({'override_original','OverrideOriginal'},true,[],'Override original data. If checked, the original signals will be replaced by the recovery.'), ...
    arg_norep('M'), ...
    arg_norep('channel_mask'));

% calculate spatial filter matrix, if necessary
if ~exist('M','var')
    [dummy,dummy,channel_mask,dummy,dummy,dummy,normField,leadField] = hlp_diskcache('filterdesign',@calc_beamformer_constraints,{signal.chanlocs.labels},roi_labels,eye(signal.nbchan),reference_type); %#ok<ASGLU>
    if strcmp(field_type,'normal')
        LF = normField;        
    elseif strcmp(field_type,'axial')
        LF = leadField(:,:);
    else
        error('Unsupported field type requested.');
    end
    M = LF'/cov(signal.data(channel_mask,:)'); 
end

% apply M
signal = utl_register_field(signal,'timeseries','srcpot',reshape(M*signal.data(channel_mask,:,:),[],size(signal.data,2),size(signal.data,3)));
signal.etc.roi_labels = roi_labels; 

if override_original
    % override signal.data & relabel channels
    signal.data = signal.srcpot;
    if size(signal.data,1) == signal.nbchan
        signal.chanlocs = struct('labels',roi_labels);
    else
        signal.nbchan = size(signal.data,1);
        signal.chanlocs = hlp_nanocache('cached_labels',10,@make_labels,roi_labels); 
    end
end

% append the M and ok arguments to the online expression
exp_endfun('append_online',{'M',M,'channel_mask',channel_mask});

function chanlocs = make_labels(roi_labels)
chanlocs = struct('labels',[cellfun(@(f){sprintf('%s_X',f)},roi_labels) cellfun(@(f){sprintf('%s_Y',f)},roi_labels) cellfun(@(f){sprintf('%s_Z',f)},roi_labels)]);
