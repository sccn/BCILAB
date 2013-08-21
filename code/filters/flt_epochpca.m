function signal = flt_epochpca(varargin)
% Apply a principal component decomposition across time (or frequency) in an epoch.
% Signal = flt_epochpca(Signal, RetainDimensions)
%
% This allows to reduce the dimensionality of the data time course or spectrum.
%
% In:
%   Signal     : Epoched data set to be processed
%
%   RetainDimensions : Reduce the dimensionality (per channel) to this number, or [] to retain all 
%                      dimensions (default: [])
%
% Out: 
%   Signal  :   processed data set
%
% See also:
%   flt_epochica
%
% TODO: 
%   Allow robust covariance
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2011-04-27

% flt_epochpca_version<0.9> -- for the cache

if ~exp_beginfun('filter') return; end

% requires epoched data, works best on spatially filtered data
declare_properties('name','EpochPCA', 'depends','set_makepos', 'follows',{'flt_fourier'}, 'independent_channels',true, 'independent_trials',false);

arg_define(varargin,...
    arg_norep({'signal','Signal'}), ...
    arg({'reducedim','RetainDimensions'}, [], [], 'Per-channel dimensionality reduction. Reduce to this number of dimensions.'), ...
    arg_norep('mapping',unassigned));

if ~exist('mapping','var')
    if isempty(reducedim)
        reducedim = signal.pnts; end
    % compute principal eigenspace for each channel
    mapping = cell(1,signal.nbchan);
    for c = 1:signal.nbchan
        try
            [mapping{c},D] = eig(hlp_diskcache('filterdesign',@cov_shrink,squeeze(signal.data(c,:,:))'));
        catch
            mapping{c} = eye(size(signal.data,2));
        end
        mapping{c} = mapping{c}(:,1:reducedim);
    end
end

% transform each channel into the eigenspace and write back
tmp = zeros(signal.nbchan,size(mapping{1},2),signal.trials);
for c = 1:signal.nbchan
    tmp(c,:,:) = (mapping{c}' * reshape(signal.data(c,:,:),size(mapping{c},1),[])); end
signal.data = tmp;
signal.pnts = size(signal.data,2);


% keep track of the mapping decomposition
global tracking;
tracking.inspection.epochpca = mapping;

exp_endfun('append_online',{'mapping',mapping});