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
    arg({'reducedim','RetainDimensions'}, [], [0 Inf], 'Per-channel dimensionality reduction. Reduce to this number of dimensions.','shape','scalar'), ...
    arg_norep('mapping',unassigned));

if ~exist('mapping','var')
    if isempty(reducedim)
        reducedim = signal.pnts; end
    for f = utl_timeseries_fields(signal)
        X = signal.(f{1});
        if ~isempty(X)
            mapping.(f{1}) = cell(1,size(X,1));
            % compute principal components for each channel
            for c = size(X,1):-1:1
                try
                    [mapping.(f{1}){c},D] = eig(hlp_diskcache('filterdesign',@cov_shrink,squeeze(X(c,:,:))')); %#ok<NASGU>
                catch %#ok<CTCH>
                    mapping.(f{1}){c} = eye(size(X,2));
                end
                mapping.(f{1}){c} = mapping.(f{1}){c}(:,1:reducedim);
            end
        end
    end
end

% transform each channel into the eigenspace and write back
for f = utl_timeseries_fields(signal)
    X = signal.(f{1});
    if ~isempty(X)
        M = mapping.(f{1});
        tmp = zeros(size(X,1),size(M{1},2),size(X,3));
        for c = 1:size(X,1)
            tmp(c,:,:) = (M{c}' * reshape(X(c,:,:),size(M{c},1),[])); end
        signal.(f{1}) = tmp;
    end
end
signal.pnts = size(signal.data,2);
signal.etc.epochpca = mapping;

exp_endfun('append_online',{'mapping',mapping});