function signal = set_fit_dipoles(varargin)
% Fit dipoles for each independent component of an IC-decomposed dataset.
%
% This function uses the Dipfit plugin for EEGLAB to automatically derive dipole locations
% for all IC's in the data set. The returned coordinates are in the MNI coordinate system.
% In addition, (probabilistic) anatomical labels can be looked up from different brain atlases [1,2].
% 
%
% In:
%   Signal           : data set with valid IC decomposition
%
%   HeadModel        : head model file (see pop_dipfit_settings) (default: standard BEM volume)
%
%   MRImage          : anatomical MR head image file (default: standard BEM image)
%
%   ChannelLocations : coregistered channel locations file (default: standard 10-20 locations)
%
%   LookupLabels     : whether to look up anatomical labels & probabilities (default: true)
%
%   ConfusionRange   : radius (in mm) of uncertainty over which to scan for labels (default: 4)
%
%   BrainAtlas       : brain atlas to use. Talairach has a larger repertoire of areas, but is (in 
%                      this version) non-probabilistic. The LONI LBPA40 atlas is a high-quality 
%                      probabilistic atlas.
%
% Out:
%   Signal : data set with added dipole annotations
%
% Notes:
%   In the default settings, a standard Eurasian BEM head model amd MR image is used, and channel
%   locations are assumed to be standard 10-20 locations if not specifically given.
%
% Examples:
%   % for a data set which is annotated with independent components, e.g. via
%   eeg = flt_ica(raw);
%
%   % ... fit dipoles for each component using default settings, and look up the nearby anatomical 
%   % structures and their hit probabilities)
%   eeg = set_fit_dipoles(eeg)
%
%   % ... fit dipoles and and look up anatomical structures within a radius of 10mm around each 
%   % dipole coordinate (passing arguments by name)
%   eeg = set_fit_dipoles('Signal',eeg,'ConfusionRange',10)
%
%   % ... fit dipoles but do not look up anatomical structures
%   eeg = set_fit_dipoles('Signal',eeg,'LookupLabels',false)
%
%   % ... fit dipoles using a (subject-) specific head model and MR file
%   eeg = set_fit_dipoles('Signal',eeg,'HeadModel','volXY.mat','MRImage','mriXY.mat')
%   
%   % ... fit dipoles using a (subject-) specific head model and MR file, as well as digitized
%   % channel locations
%   eeg = set_fit_dipoles('Signal',eeg,'HeadModel','my_vol.mat','MRImage','my_mri.mat','ChannelLocations','my_locs.mat')
%
% References:
%   [1] Lancaster JL, Woldorff MG, Parsons LM, Liotti M, Freitas CS, Rainey L, Kochunov PV, Nickerson D, Mikiten SA, Fox PT, "Automated Talairach Atlas labels for functional brain mapping". 
%       Human Brain Mapping 10:120-131, 2000
%   [2] Shattuck DW, Mirza M, Adisetiyo V, Hojatkashani C, Salamon G, Narr KL, Poldrack RA, Bilder RM, Toga AW, "Construction of a 3D probabilistic atlas of human cortical structures."
%       Neuroimage. 2008 39(3):1064-80
%
% See also:
%   flt_ica, coregister, pop_multifit
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2011-01-21

% set_fit_dipoles_version<0.91> -- for the cache

if ~exp_beginfun('filter') return; end

% declare GUI name, etc.
declare_properties('name',{'DipoleFitting','dipfit'}, 'depends','flt_ica', 'precedes','flt_project', 'independent_channels',true, 'independent_trials',true);

% define arguments...
arg_define(varargin, ...
    arg_norep({'signal','Signal'}), ...
    arg({'hdm_file','HeadModel'}, '', [],'Head model file. The BEM head model that should be used by dipfit (using a standard model if empty).'),...
    arg({'mri_file','MRImage'}, '', [],'MR head image file. The MR head image that should be used by dipfit (using a standard image if empty).'),...
    arg({'chan_file','ReferenceLocations'}, '', [],'Coregistered reference chanlocs. The coregistered electrode locations that should be used as reference for alignment warping.'),...
    arg({'lookup_labels','LookupLabels'},true,[], 'Look up dipole labels. If enabled, look up anatomical labels from the selected brain atlas.'),...
    arg({'confusion_range','ConfusionRange'}, 4, [1 4 15 30], 'Confusion Radius (mm). Assumed standard deviation of uncertainty in the dipole fits, in millimeters.'),...    
    arg({'brain_atlas','BrainAtlas'}, 'Talairach', {'Talairach','LBPA40'}, 'Brain atlas to use. Talairach has a larger repertoire of areas, but is (in this version) non-probabilistic. The LONI LBPA40 atlas is a high-quality probabilistic atlas; however, you have to to download it yourself.'),...
    arg({'var_threshold','VarianceThreshold'}, 15, [1 10 25 100], 'Residual Variance Thresold (%). The threshold in dipole fit residual variance above which dipoles models are rejected.'),...    
    arg({'mri_constraints','UseMRIConstraints'}, false, [], 'Use MRI constraints. If enabled, restricts dipole solutions to the grey-matter volume. The currently requires the SPM toolbox to be in the path.'),...
    arg({'discard_nonlocalizable','DiscardNonlocalizable'}, false, [], 'Discard non-localizable components.'),...
    arg({'verbose','VerboseOutput'}, false, [], 'Verbose output.'),...
    arg_norep('dipfit_info'));

if ~exist('dipfit_info','var')
    % use standard data if unspecified
    if isempty(hdm_file) %#ok<*NODEF>
        hdm_file = env_translatepath('resources:/standard_BEM/standard_vol.mat'); end
    if isempty(mri_file)
        if mri_constraints
            mri_file = env_translatepath('resources:/standard_BEM/standard_mri_gm.mat'); 
        else
            mri_file = '';
        end
    end
    if isempty(chan_file)
        % figure out the labeling scheme to determine the correct coregistration...
        if ~isfield(signal.chaninfo,'labelscheme')
            signal.chaninfo.labelscheme = '10-20'; end
        switch signal.chaninfo.labelscheme
            case '10-20'
                chan_file = env_translatepath('resources:/standard_BEM/elec/standard_1005.elc');
            case 'sccn_128_v1'
                chan_file = env_translatepath('resources:/sccn_BEM_coregistered_128_v1.xyz');
            case 'sccn_128_v2'
                chan_file = env_translatepath('resources:/sccn_BEM_coregistered_128_v2.xyz');
            case 'sccn_256_v1'
                chan_file = env_translatepath('resources:/sccn_BEM_coregistered_256_v1.xyz');
            case 'sccn_256_v2'
                chan_file = env_translatepath('resources:/sccn_BEM_coregistered_256_v2.xyz');
            otherwise
                error('Unknown channel labeling scheme -- please supply a chan_file argument.');
        end
    end
    
    disp(['Now fitting dipoles... (montage reference: ' chan_file ')']);
    % fit icaweights
    if verbose
        dipfit_info = hlp_diskcache('dipfits',@do_fitting,signal,mri_file,hdm_file,chan_file,lookup_labels,brain_atlas,confusion_range,var_threshold);
    else
        [text,dipfit_info] = evalc('hlp_diskcache(''dipfits'',@do_fitting,signal,mri_file,hdm_file,chan_file,lookup_labels,brain_atlas,confusion_range,var_threshold)');
    end
    try
        % check if we need to fit Amica models, too...
        if isfield(signal.etc,'amica')
            sig = signal;
            multimodel = {};
            for m=1:size(signal.etc.amica.W,3)
                sig.icaweights = signal.etc.amica.W(:,:,m);
                sig.icawinv = signal.etc.amica.A(:,:,m);
                sig.icasphere = signal.etc.amica.S;
                if verbose
                    tmp = hlp_diskcache('dipfits',@do_fitting,sig,mri_file,hdm_file,chan_file,lookup_labels,brain_atlas,confusion_range,var_threshold);
                else
                    [text,tmp] = evalc('hlp_diskcache(''dipfits'',@do_fitting,sig,mri_file,hdm_file,chan_file,lookup_labels,brain_atlas,confusion_range,var_threshold);');
                end
                multimodel{m} = tmp.model; %#ok<AGROW>
            end
            dipfit_info.multimodel = multimodel;
        end
    catch
        disp('Could not compute dipfit solutions for Amica models.');
    end
end
    
signal.dipfit = dipfit_info;
if discard_nonlocalizable
    retain_ics = ~cellfun(@isempty,{signal.dipfit.model.posxyz});
    % update the dipole model
    signal.dipfit.model = signal.dipfit.model(retain_ics);    
    % restrict the ICs & some derived data
    signal.icaweights = signal.icaweights(retain_ics,:);    
    signal.icawinv = signal.icawinv(:,retain_ics);
    if ~isempty(signal.icaact) && ~isscalar(signal.icaact)
        signal.icaact = signal.icaact(retain_ics,:,:); end
end

if isfield(signal.etc,'amica')
    signal.etc.amica.dipfit = dipfit_info; end


global tracking;
tracking.inspection.dipfit = dipfit_info;

% when applied online, include the dipfit info into the parameters (so it gets attached immediately)
exp_endfun('append_online',{'dipfit_info',dipfit_info});

function result = do_fitting(signal,mri_file,hdm_file,chan_file,lookup_labels,brain_atlas,confusion_range,var_threshold)

% coregister chanlocs to reference channels
[dummy,warping] = coregister(signal.chanlocs,chan_file,'warp','auto','manual','off'); %#ok<ASGLU>
% generate the dipfit info
%nosedir needs to be set to +X here for dipfit to work
signal.chaninfo.nosedir = '+X';
tmp = pop_multifit(pop_dipfit_settings(pop_select(signal,'channel',signal.icachansind),  ...
    'mrifile',mri_file, 'hdmfile',hdm_file,'chanfile',chan_file,'coordformat','MNI','coord_transform',warping), 1:size(signal.icaweights,1),'threshold',var_threshold);
% optionally derive anatomical labels & probabilities
if lookup_labels
    switch lower(brain_atlas)
        
        case 'talairach'
            db = org.talairach.Database;
            db.load(env_translatepath('resources:/talairach.nii'));
            for k=1:length(tmp.dipfit.model)
                try
                    p = icbm_spm2tal(tmp.dipfit.model(k).posxyz);
                    tmp.dipfit.model(k).labels = cellfun(@(d)char(d),cell(db.search_range(p(1),p(2),p(3),1.5*confusion_range)),'UniformOutput',false);
                    % and compute structure probabilities within the selected volume
                    [structures,x,idxs] = unique(hlp_split(sprintf('%s,',tmp.dipfit.model(k).labels{:}),',')); %#ok<ASGLU>
                    probabilities = mean(bsxfun(@eq,1:max(idxs),idxs'));
                    [probabilities,reindex] = sort(probabilities,'descend');
                    structures = structures(reindex);
                    mask = ~strcmp(structures,'*');
                    tmp.dipfit.model(k).structures = structures(mask);
                    tmp.dipfit.model(k).probabilities = probabilities(mask)*5; % there are 5 partitions
                    % figure('Position',[0 0 2560 900]); topoplot(tmp.icawinv(:,k),tmp.chanlocs); title([hlp_tostring(structures(mask)) 10 hlp_tostring(5*probabilities(mask))]); pop_dipplot(tmp,k);
                catch
                    tmp.dipfit.model(k).labels = {};
                    tmp.dipfit.model(k).structures = {};
                    tmp.dipfit.model(k).probabilities = [];
                end
            end
            
        case 'lbpa40'
            mask = find(~cellfun('isempty',{tmp.dipfit.model.posxyz}));
            % build query
            coords = vertcat(tmp.dipfit.model.posxyz);
            coords = [coords confusion_range * ones(size(coords,1),1)];
            % look up from atlas (note: slow!)
            [probs,labels] = label_dipoles(coords);
            % determine overlapped specialty regions to match up with Talairach's labeling scheme
            L = ~cellfun('isempty',strfind(labels,' L '));
            R = ~cellfun('isempty',strfind(labels,' R '));
            B = ~cellfun('isempty',strfind(labels,'Brainstem'));
            C = ~cellfun('isempty',strfind(labels,'Cerebellum'));
            allprobs = [sum(probs(:,L),2) sum(probs(:,R),2) sum(probs(:,C),2)*[0.5 0.5] sum(probs(:,B),2)];
            allstructs = {'Left Cerebrum' 'Right Cerebrum' 'Left Cerebellum' 'Right Cerebellum' 'Brainstem'};
            % go through all gyri and add up left & right probabilities
            gyri = unique(cellfun(@(l)l(12:end),labels(1:end-2),'UniformOutput',false))';
            allstructs = [allstructs gyri];
            for g=1:length(gyri)
                curgyrus = gyri{g};
                matches = ~cellfun('isempty',strfind(labels,curgyrus));
                allprobs = [allprobs sum(probs(:,matches),2)];
            end
            for k=1:length(mask)
                % retain only those with non-zero probability
                sel = allprobs(k,:) ~= 0;
                probabilities = allprobs(k,sel);
                structures = allstructs(sel);
                % sort probs & associated structs by descending probability
                [probabilities,reindex] = sort(probabilities,'descend');
                structures = structures(reindex);
                % store in the model...
                tmp.dipfit.model(mask(k)).structures = structures;
                tmp.dipfit.model(mask(k)).probabilities = probabilities;
                % figure('Position',[0 0 2560 900]);topoplot(tmp.icawinv(:,k),tmp.chanlocs); title([hlp_tostring(structures) 10 hlp_tostring(probabilities)]); pop_dipplot(tmp,k);
            end
    end
end
result = tmp.dipfit;
