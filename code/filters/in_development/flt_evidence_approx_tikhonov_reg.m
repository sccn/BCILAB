function [signal, state] = flt_evidence_approx_tikhonov_reg(varargin)
% Return the current source density for a given head model and data using
% the cortically-constrained LORETA (low resolution electrical
% tomographic analysis) with a Bayesian update scheme for hyperparameters.
% The reconstructed CSD time-series (or source potential maps) will be 
% stored in signal.srcpot. This matrix has dimension [num_voxels x num_samples].
% 
% Author: Tim Mullen, Jan 2013, SCCN/INC/UCSD
%         Alejandro Ojeda, Jan 2013, SCCN/INC/UCSD
%         Christian Kothe, Jan 2013, SCCN/INC/UCSD


if ~exp_beginfun('filter'), return; end

declare_properties('name','LORETA', 'experimental',true, 'independent_channels',false, 'independent_trials',false);

arg_define(varargin, ...
    arg_norep({'signal','Signal'}), ...
    arg_nogui({'K','ForwardModel'},[],[],'Forward model (matrix)','shape','matrix'), ...
    arg_nogui({'L','LaplacianOperator'},[],[],'Laplacian operator. Sparse matrix of N sources x N sources, this is matrix is used as the square root of the precision matrix of the sources.'), ...
    arg_sub({'options','LoretaOptions'},{},...
        { ...
        arg({'maxTol','MaxTolerance'},1e-6,[0 Inf],'Tolerance for hyperparameter update loop','cat','Loreta Options'), ...
        arg({'maxIter','MaxIterations'},100,[1 Inf],'Maximum iterations for hyperparameter update loop','cat','Loreta Options'), ...
        arg({'verbose','VerboseOutput'},false,[],'Verbosity','cat','Loreta Options'), ...
        arg({'block_size','BlockSize'},5, [], 'Block granularity for processing. The inverse operator will be updated using blocks of this many samples. This assumes that the inverse solution is spatially stationary over this many samples.'), ...
        arg({'skipFactor','SkipFactor'},0,[0 Inf],'Number of blocks to skip'), ...
        arg_deprecated({'maxblocks','MaxBlocks'},Inf,[0 Inf],'Maximum number of blocks'), ...
        arg({'standardize','Standardize'},'all',{'none','channels','time','all'},'Rescale data to unit variance. If ''channels'', standardization is carried out across channels for each time point. If ''all'' each data sample is normalized by the standard deviation taken over all data.'), ...
        arg({'useGPU','UseGPU'},false,[],'Use GPU to accelerate computation.'), ...
        },'Additional options for Loreta function'), ...
    arg({'verb','Verbosity'},0,[],'Verbosity level.','typecheck',false), ...
    arg_nogui({'state','State'},[],[],'State object. When provided, hyperparameters will be estimated adaptively from prior state'));
if verb
    fprintf('Estimating current source density using the Evidence Approximation (%s)\n',mfilename); 
end

[nchs, npnts, ntrs] = size(signal.data);
if isempty(options.block_size) || options.block_size > npnts
    options.block_size = npnts;
end
numsplits    = floor(npnts/options.block_size);

% if necessary, cast to double-precision
if ~strcmpi(class(signal.data),'double')
    signal.data = double(signal.data);
end
    
% options the data
if ~strcmpi(options.standardize,'none')
    switch options.standardize
        case 'channels'
            scale = std(signal.data,[],1);
        case 'time'
            scale = std(signal.data,[],2);
        case 'all'
            scale = std(signal.data(:));
    end
    data = bsxfun(@rdivide,signal.data,scale+eps);
else
    data = signal.data;
end
data(isnan(data)) = 0;
data(isinf(data)) = 0;

if isempty(state) || ~isfield(state,'iLV') || isempty(state.iLV)
    if verb
        fprintf('...computing SVD of LFM.\n');
    end
    % mode is offline or we are initializing online filter
    % perform one-time SVD for faster computation.
    [U,S,V]      = svd(K/L,'econ');
    state.iLV    = L\V;
    state.s2     = diag(S).^2; %s^2
    state.Ut     = U';
    state.alpha = [];
    state.beta  = [];
    options.maxTol(options.maxTol<1e-6) = 1e-6;
end
 
if npnts == 0
    % no data
    signal.srcpot    = [];
    state.srcweights = [];
    exp_endfun; return;
end

signal.srcpot    = zeros([size(K,2), npnts, ntrs]);
state.srcweights = zeros(size(L,1),nchs);
sum_srcweights   = zeros(size(L,1),nchs);
signal.loretaHistory = struct([]);
tmpHist = [];
if verb
    fprintf('...assuming %d stationary blocks of length %d\n',numsplits,options.block_size);
end


% loop over all trials
for tr=1:ntrs
    if verb && ntrs > 1
        fprintf('\nTrial (%d\\%d).',tr,ntrs);
    end
    k = 0;
    % loop over sub-blocks and estimate CSD for each block
    for i=0:options.skipFactor+1:numsplits-1
        if verb
            if i+1 >= floor(numsplits*(k+1)/10)
                k = k + 1;
                fprintf('%0.3g%%...',round((i/numsplits)*100));
            end
        end
        range = 1+floor(i*npnts/numsplits) : min(npnts,floor((i+1)*npnts/numsplits));
        % call (dynamic bayesian) loreta estimator
        [signal.srcpot(:,range,tr), state.alpha, state.beta, state.srcweights] = ...
            evidence_approx_tikhonov_reg_svd(data(:,range,tr), state.Ut, state.s2, state.iLV,state.alpha, state.beta,options);
        if ~isempty(tmpHist)
            signal.loretaHistory{tr} = [signal.loretaHistory{tr},tmpHist]; 
        end

        if options.skipFactor > 0
            % estimate CSD for samples between blocks using current inverse operator
            range = 1+floor((i+1)*npnts/numsplits) : min(npnts,floor((i+options.skipFactor+1)*npnts/numsplits));
            signal.srcpot(:,range,tr) = state.srcweights*signal.data(:,range,tr);
        end
        
        % running sum
        sum_srcweights = sum_srcweights + state.srcweights;
    end
end


if numsplits > 1
    % store the mean inverse operator over all splits          
    state.srcweights = sum_srcweights/ceil(numsplits*ntrs/(options.skipFactor+1));
end

if ~strcmpi(options.standardize,'none')
    % recale data to original units
    signal.srcpot = bsxfun(@times,signal.srcpot,scale);
end
if verb, fprintf('done.\n');end

exp_endfun;
