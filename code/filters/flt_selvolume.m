function signal = flt_selvolume(varargin)
% Select independent components according to what brain volumes they are in.
% Signal = flt_selvolume(...)
%
% In:
%   Signal        : input data set, assumed to have an associated IC decomposition
%
%   Hemisphere    : Cell array of hemisphere names to retain (or empty to retain all) (default: {})
%
%   Lobe          : Cell array of lobe names to retain (or empty to retain all) (default: {})
%
%   Gyrus         : Cell array of gyrus names to retain (or empty to retain all) (default: {})
%
%   TransformData : Return the result as data (instead of as IC decomposition).
%                   This allows to use volume-based selection together with with several archaic 
%                   signal processing plugins (or paradigms) that are not aware of advanced 
%                   meta-data.
%
% Out:
%   Signal : input data set restricted to dipoles which lie in the respective areas
%
% Notes:
%   Requires that set_fit_dipoles (and flt_ica) have been run before.
%
% Examples:
%   % after a data set has been annotated with an ICA decomposition, and dipoles have been
%   % fit for the components, e.g., as follows:
%   eeg = set_fit_dipoles(flt_ica(eeg));
%
%   % ... retain only components that are located in the left or right Cerebrum
%   eeg = flt_selvolume(eeg,{'Left Cerebrum','Right Cerebrum'})
%
%   % ... retain only components that are located in the Occipital lobe
%   eeg = flt_selvolume(eeg,true,{'Occipital Lobe'})
%
%   % ... retain only components that are located both in the Occipital lobe and left Cerebrum
%   eeg = flt_selvolume(eeg,{'Left Cerebrum'},{'Occipital Lobe'})
%
%   % ... retain all components that are located in the Cingulate Gyrus or Anterior Cingulate
%   % passing all arguments by their name
%   eeg = flt_selvolume('Signal',eeg,'Gyrus',{'Anterior Cingulate','Cingulate Gyrus'})
%
%   % ... retain only components that are located in the left or right Cerebrum and replace
%   % the .data field of the data set with the retained component activations
%   eeg = flt_selvolume('Signal',eeg,'Hemisphere',{'Left Cerebrum','Right Cerebrum'},'TransformData',true)
%
%   % ... retain all ICs (i.e., do noting)
%   eeg = flt_selvolume(eeg)
%
% See also:
%   set_fit_dipoles, flt_ica
%
% TODO:
%   Come up with a better way to find a small number of well-dispersed components in ROIs?
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2011-01-20

% flt_selvolume_version<0.85> -- for the cache

if ~exp_beginfun('filter') return; end

% declare GUI name, etc.
declare_properties('name','VolumeSelection', 'depends','set_fit_dipoles', 'precedes',{'set_makepos','flt_project'}, 'independent_channels',false, 'independent_trials',true);

% define arguments...
arg_define(varargin, ...
    arg_norep({'signal','Signal'}), ...    
    arg({'sel_hemi','Hemisphere'},true, {'Left Cerebrum','Right Cerebrum','Left Cerebellum','Right Cerebellum','Left Brainstem','Right Brainstem','Inter-Hemispheric'}, 'Hemisphere to retain. Restrict ICs to those that fall in the specified hemispheres.'),...
    arg({'sel_lobe','Lobe'},true, {'Anterior Lobe','Frontal Lobe','Frontal-Temporal Space','Limbic Lobe','Medulla','Midbrain','Occipital Lobe','Parietal Lobe','Pons','Posterior Lobe','Sub-lobar','Temporal Lobe'}, 'Lobe to retain. Restrict ICs to those that fall in the specified lobes.'),...
    arg({'sel_gyrus','Gyrus'},true, {'Angular Gyrus','Anterior Cingulate','Caudate','Cerebellar Lingual','Cerebellar Tonsil','Cingulate Gyrus','Claustrum','Culmen','Culmen of Vermis','Cuneus','Declive','Declive of Vermis','Extra-Nuclear','Fastigium','Fourth Ventricle','Fusiform Gyrus','Inferior Frontal Gyrus','Inferior Occipital Gyrus','Inferior Parietal Lobule','Inferior Semi-Lunar Lobule','Inferior Temporal Gyrus','Insula','Lateral Ventricle','Lentiform Nucleus','Lingual Gyrus','Medial Frontal Gyrus','Middle Frontal Gyrus','Middle Occipital Gyrus','Middle Temporal Gyrus','Nodule','Orbital Gyrus','Paracentral Lobule','Parahippocampal Gyrus','Postcentral Gyrus','Posterior Cingulate','Precentral Gyrus','Precuneus','Pyramis','Pyramis of Vermis','Rectal Gyrus','Subcallosal Gyrus','Sub-Gyral','Superior Frontal Gyrus','Superior Occipital Gyrus','Superior Parietal Lobule','Superior Temporal Gyrus','Supramarginal Gyrus','Thalamus','Third Ventricle','Transverse Temporal Gyrus','Tuber','Tuber of Vermis','Uncus','Uvula','Uvula of Vermis'}, 'Gyri to retain. Restrict ICs to those that fall in the specified gyri.'),...
    arg({'probability_cutoff','ProbabilityCutoff'},0.7,[0 1],'Minimum probability. If a component has less than this probability of being in the selected structures, it will be removed.'),...
    arg({'do_transform','TransformData','transform'},false,[],'Transform the data rather than annotate. By default, ICA decompositions are maintained as annotations to the data set, but several algorithms operate by default on the raw data, and are unaware of these annotations.'),...
    arg({'max_components','MaxComponents'},Inf,[],'Retain at most this many components. This allows to prune the dimensionality of the source space.','experimental',true),...
    arg_norep('retain_ics', unassigned));

if ~exist('retain_ics','var')    
    
    if isfield(signal.dipfit,'multimodel') && length(signal.dipfit.multimodel) > length(signal.dipfit.model) %#ok<NODEF>
        % this is, among others, because .amica.W is currently a 3d array which
        % cannot be pruned separately for each model
        error('Volume selection can currently not be applied to multi-model decompositions.');
    end
    
    retain_ics = true(1,length(signal.dipfit.model));
    multiplied_probs = zeros(1,length(signal.dipfit.model));
    
    % subset ICs based on probability under retained brain volumes
    for ic=1:length(retain_ics)
        multiplied_prob = 1;
        % for each of the three partitions...
        for partition = {sel_hemi,sel_lobe,sel_gyrus}
            summed_prob = 0;
            % for each of the checked labels in the partition...
            for checked_label = partition{1}
                % if this label shows up in any chunk of the labels for the given IC, we're okay
                mask = ~cellfun('isempty',strfind(signal.dipfit.model(ic).structures,checked_label{1}));
                summed_prob = summed_prob + sum(signal.dipfit.model(ic).probabilities(mask));
            end
            multiplied_prob = multiplied_prob*summed_prob;
        end
        % if, for any of the three partitions, none of the checked labels shows up in the IC's
        % labels, then the IC is dropped
        if multiplied_prob < probability_cutoff
            retain_ics(ic) = false; end
        multiplied_probs(ic) = multiplied_prob;
    end
        
    % optionally prune excess components
    if max_components < nnz(retain_ics)
        % calculate a measure of spatial isolation for each component
        isolation_measure = zeros(1,length(signal.dipfit.model));
        for ic=find(retain_ics)
            isolation_measure(ic) = min(sqrt(sum((bsxfun(@minus,signal.dipfit.model(ic).posxyz,vertcat(signal.dipfit.model(setdiff(find(retain_ics),ic)).posxyz))).^2,2))); end
        isolation_measure(retain_ics) = log(1+isolation_measure(retain_ics)); %#ok<NASGU>
        % retain the highest-probability ones
        [dummy,order]= sort(multiplied_probs(retain_ics),'descend'); %#ok<ASGLU>
        %[dummy,order]= sort(isolation_measure(retain_ics).*multiplied_probs(retain_ics),'descend'); %#ok<ASGLU>
        retain_ics = find(retain_ics);
        retain_ics = sort(retain_ics(order(1:max_components)));
    end
end


% update the dipole model
signal.dipfit.model = signal.dipfit.model(retain_ics);

% restrict the ICs & some derived data
signal.icaweights = signal.icaweights(retain_ics,:);
signal.icawinv = signal.icawinv(:,retain_ics);
if ~isempty(signal.icaact)
    signal.icaact = signal.icaact(retain_ics,:,:); end

% and transform the data itself, if desired
if do_transform
    signal.data = (signal.icaweights*signal.icasphere)*signal.data(signal.icachansind,:);
    signal.nbchan = size(signal.data,1);
    signal.chanlocs = struct('labels',cellfun(@num2str,num2cell(1:signal.nbchan,1),'UniformOutput',false));
    signal.icaweights = [];
    signal.icasphere = [];
    signal.icawinv = [];
    signal.icachansind = [];
    signal.dipfit = [];
end

% to replicate this processing step online, we directly append the list of IC's to be retained to the arguments of the function
exp_endfun('append_online',{'retain_ics',retain_ics});

