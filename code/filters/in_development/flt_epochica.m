function signal = flt_epochica(varargin)
% Apply an independent component decomposition across time (or frequency) in an epoch.
% Signal = flt_epochica(Signal)
%
% This is currently experimental.
%
% In:
%   Signal     : Epoched data set to be processed
%
%   RetainDimensions : Reduce the dimensionality (per channel) to this number, or [] to retain all 
%                      dimensions (default: [])
%
%   Fast : use a fast approach (instead of a slow one) (default: false)
%
% Out: 
%   Signal  :   processed data set
%
% See also:
%   flt_epochpca, flt_ica
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2011-04-27

% flt_epochica_version<0.9> -- for the cache

if ~exp_beginfun('filter') return; end

% requires epoched data, works best on spatially filtered data
declare_properties('name','EpochICA', 'experimental',true, 'depends','set_makepos', 'follows',{'flt_fourier','flt_epochpca'}, 'independent_channels',false, 'independent_trials',false);

arg_define(varargin,...
    arg_norep({'signal','Signal'}), ...
    arg({'reducedim','RetainDimensions'}, [], [], 'Dimensionality reduction. Reduce to this number of dimensions.','shape','scalar'), ...
    arg({'fastmode','FastApproach'}, false, [], 'Use FastICA. As a faster alternative to the (high-quality) amica solution.'), ...
    arg_norep('mapping',unassigned), ...
    arg_norep('scaling',unassigned));

if ~exist('mapping','var')
    % mapping needs to be created first
    [C,S,T] = size(signal.data);
    
    % concatenate the data across epochs & channels
    tmp = zeros(S,C*T);
    for s = 1:S
        vec = squeeze(signal.data(:,s,:));
        tmp(s,:) = vec(:); 
    end

    % standardize the data
    scaling = hlp_findscaling(tmp','std');
    tmp = hlp_applyscaling(tmp',scaling)';
    
    if fastmode
        % do a FastICA
        [mixing,mapping] = fastica(tmp,'maxNumIterations',1000,'approach','symm','numOfIC',reducedim,'g','tanh','finetune','tanh','a1',1,'a2',1,'mu',1,...
            'stabilization','off','epsilon',0.001,'maxFinetune',100,'sampleSize',1,'initGuess',[],'verbose','off','displayMode','off','displayInterval',1,'firstEig',1,'lastEig',[]); %#ok<ASGLU>
    else
        % do an AMICA
        tmpset = exp_eval(set_new('data',tmp,'srate',C/S * signal.srate));
        tmpdecomp = exp_eval(flt_ica(tmpset,{'amica','useqsub','on','num_models',1,'max_iter',250},'clean','hardcore_nochans'));
        mapping = tmpdecomp.icaweights*tmpdecomp.icasphere;
    end
end

% transform each channel into the spectral eigenspace and write back
tmp = zeros(signal.nbchan,size(mapping,1),signal.trials);
for c = 1:signal.nbchan
    tmp(c,:,:) = (mapping * hlp_applyscaling(squeeze(signal.data(c,:,:))',scaling)'); end

% mapping is not yet defined: obtain it via ICA
signal.data = tmp;
signal.pnts = size(signal.data,2);
% make sure that the mapping gets overridden by the precomputed one online
exp_endfun('append_online',{'mapping',mapping,'scaling',scaling});


% keep track of the mapping decomposition
global tracking;
tracking.inspection.epochica = mapping;
