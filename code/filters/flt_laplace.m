function signal = flt_laplace(varargin)
% Applies a simple Hjorth-style surface laplacian filter.
% function Signal = flt_laplace(Signal,NeighbourCount)
%
% The surface laplacian [1] is a spatial filter in which from the signal of each channel, the
% (spatially) averaged signals of its N closest neighbours are subtracted (N usually 4 or 8). This
% implements a spatial high-pass filter, which dampens large-scale scalp signals and amplifies
% localized signals. The primary benefit of this spatial filter over more advanced ones such as
% Independent Component Analysis (flt_ica) or Common Spatial Patterns (e.g., para_csp) is that it is
% the simplest, most robust and quickest (to compute) filter. Its major downside is that it is far
% from optimal for almost any given task, and the predictive performance of paradigms relying on it
% will therefore usually be sub-optimal.
%
% In some cases, inclusion of a surface laplacian filter has no effect on the output of a BCI,
% namely when it is directly followed by certain variants of linear adaptive mappings, such as
% unregularized LDA, CSP, PCA or ICA.
%
% In:   
%   Signal         : EEGLAB data set, either continuous or epoched
%
%   NeighbourCount : number of neighbor directions to consider for each channel (default: 8)
%
% Out:
%   Signal : laplacian-filtered EEG set
%
% Examples:
%   % use default settings
%   eeg = flt_laplace(eeg)
%
%   % use 4 neighbour channels
%   eeg = flt_laplace(eeg,4)
%
%   % pass all arguments by name
%   eeg = flt_laplace('Signal',eeg,'NeighbourCount',4)
%
% References:
%  [1] Hjorth, B. "An on-line transformation of EEG scalp potentials into orthogonal source derivations."
%      Electroencephalography and Clinical Neurophysiology, 39 (1975), 526-530.
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-06-29

% flt_laplace_version<1.0> -- for the cache

if ~exp_beginfun('filter') return; end

declare_properties('name','SurfaceLaplacian', 'follows','flt_selchans', 'cannot_follow','flt_ica', 'independent_channels',false, 'independent_trials',true);

arg_define(varargin, ...
    arg_norep({'signal','Signal'}), ...
    arg({'neighs','NeighbourCount'}, 8, [3 16], 'Number of neighbors per channel. Typical values are 4-8.'),...
    arg_norep('M',unassigned), ...
    arg_norep('ok',unassigned));

% in some cases, we may already know M (namely in the online case)
if ~exist('M','var')
    % look up standard locations if necessary
    if ~all(isfield(signal.chanlocs,{'X','Y','Z'})) || all(cellfun('isempty',{signal.chanlocs.X}))
        signal.chanlocs = pop_chanedit(signal.chanlocs,'lookup','standard-10-5-cap385.elp'); end
    
    % map admissible ("ok") channels into a spherical coordinate system (that has no seams across
    % the scalp)
    ok = find(~(cellfun('isempty',{signal.chanlocs.X}) | cellfun('isempty',{signal.chanlocs.Y}) | cellfun('isempty',{signal.chanlocs.Z})));
    [px,py] = cart2sph([signal.chanlocs(ok).Z],[signal.chanlocs(ok).X],[signal.chanlocs(ok).Y]);
    
    % compute neighbor matrix M
    M = zeros(length(ok));
    for c=1:length(ok)
        % get surface angles/distances to all other channels
        v = [px(c)-px; py(c)-py];
        [ang,dst] = deal((180/pi)*(pi+atan2(v(1,:),v(2,:))), sqrt(v(1,:).^2+v(2,:).^2));
        for s = 1:neighs
            % find all channels within the sector s (excl. c) and store the closest one
            validchns = find(within(ang, 360*(s-1.5)/neighs, 360*(s-0.5)/neighs) & c~=(1:length(ok)));
            [dummy,idx] = min(dst(validchns)); %#ok<ASGLU>
            if ~isempty(idx)
                M(c,validchns(idx)) = 1; end
        end
    end
    
    % finalize
    M = eye(length(ok)) - normrow(M);
end

% apply M
signal.data(ok,:,:) = reshape(M*reshape(signal.data(ok,:,:),length(ok),[]),length(ok),size(signal.data,2),size(signal.data,3));

% append the M and ok arguments to the online expression
exp_endfun('append_online',{'M',M,'ok',ok});



function tf = within(x,a,b)
if b<a b = b+360; end
tf = (x>=a & x<b) | ((x+360)>=a & (x+360)<b) | ((x-360)>=a & (x-360)<b);
