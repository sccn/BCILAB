function [signal, state] = flt_beamforming(varargin)
% Return the current source density for a given head model and data using
% a simple beamforming method.
% The reconstructed CSD time-series (or source potential maps) will be 
% stored in signal.srcpot. This matrix has dimension [num_voxels x num_samples].
% 
% Author: Tim Mullen, Jan 2013, SCCN/INC/UCSD
%         Alejandro Ojeda, Jan 2013, SCCN/INC/UCSD
%         Christian Kothe, Jan 2013, SCCN/INC/UCSD


if ~exp_beginfun('filter'), return; end

declare_properties('name','Beamforming', 'experimental',true, 'independent_channels',false, 'independent_trials',false);

arg_define(varargin, ...
    arg_norep({'signal','Signal'}), ...
    arg_nogui({'lead_field','LeadField','K','ForwardModel'},[],[],'Forward model (matrix)','shape','matrix'), ...
    arg({'cov_range','CovRange'},'whole',{'whole','epoch'},'Covariance data range. This is the range over which the covariance shall be calculated.'), ...
    arg({'cov_shrinkage','CovShrinkage'},0,[],'Covariance shrinkage parameter. For better conditioning.'), ...
    arg({'cov_robust','CovRobust'},false,[],'Use robust covariance estimate.'), ...
    arg({'rescale_sources','RescaleSources','Rescale'},false,[],'Whether to rescale the solution. If false, the source activity scale will be incorrect, but downstream methods may be scale-invariant.'), ...
    arg({'standardize_trials','StandardizeTrials'},false,[],'Perform per-trial standardization'), ...
    arg_nogui({'state','State'},[],[],'State object. When provided, hyperparameters will be estimated adaptively from prior state'));

[dummy,S,T] = size(signal.data); %#ok<ASGLU>

if strcmp(cov_range,'whole')
    % use whole-data covariance matrix
    if isempty(state)
        LF = lead_field;
        if cov_robust
            C = cov_blockgeom(signal.data(:,:)');
        else
            C = cov(signal.data(:,:)');
        end
        C = (1-cov_shrinkage)*C + cov_shrinkage*mean(trace(C))*eye(length(C));
        state.srcweights = LF'/C;
        if rescale_sources
            if ndims(LF) == 3
                error('Rescaling not yet implemented for vectorial lead-field matrices.'); end
            for k=size(LF,2):-1:1
                scales(k) = 1./((LF(:,k)'/C)*LF(:,k)); end
            state.srcweights = bsxfun(@times,scales(:),state.srcweights); 
        end
    end
    signal.srcpot = reshape(state.srcweights*signal.data(:,:),[],S,T);
elseif strcmp(cov_range,'epoch')
    % use per-epoch covariance matrix
    if isempty(signal.epoch)
        error('Trying to use epoch-wise covariance matrix, but data not epoched.'); end
    for t=T:-1:1
        LF = lead_field;
        if cov_robust
            C = cov_blockgeom(signal.data(:,:,t)');
        else
            C = cov(signal.data(:,:,t)');
        end
        C = (1-cov_shrinkage)*C + cov_shrinkage*mean(trace(C))*eye(length(C));
        srcweights = LF'/C;
        if rescale_sources
            if ndims(LF) == 3
                error('Rescaling not yet implemented for vectorial lead-field matrices.'); end
            for k=size(LF,2):-1:1
                scales(k) = 1./((LF(:,k)'/C)*LF(:,k)); end
            srcweights = bsxfun(@times,scales(:),srcweights); 
        end        
        signal.srcpot(:,:,t) = srcweights*signal.data(:,:,t); 
    end
end
if standardize_trials
    signal.srcpot = bsxfun(@minus,signal.srcpot,mean(signal.srcpot,2));
    signal.srcpot = bsxfun(@times,signal.srcpot,1./(std(signal.srcpot,[],2)+eps));
end

exp_endfun;