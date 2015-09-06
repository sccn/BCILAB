function [AR PE out] = mvar_scsa_em(varargin)

% Algorithm: SCSA EM
%
% Description:
%
% This implements the Sparsely Connected
% Sources Analysis (SCSA) algorithm of
% Haufe et al, 2010 [1]. Optimization is
% performed using an Expectation-
% Maximization (EM) algorithm.
%
%
% Author Credits:
%
% The code for SCSA was contributed by 
% Stefan Haufe, Ph.D [2]
%
% References and Code:
%
% [1] Haufe, S., Tomioka, R., & Nolte, G. (2010). Modeling sparse 
%     connectivity between underlying brain sources for EEG/MEG. 
%     Biomedical Engineering, (c), 1?10.
% [2] http://doc.ml.tu-berlin.de/publications/
%
% Dependencies: scsa_em()
%
% ------------------------------------------------------------------------
%
%
% See Also: est_fitMVAR(), mvar_dalSCSA(), mvar_ridge(), mvar_vieiramorf(),
%           est_fitMVARKalman()
%
% References:
%
% Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual.
%
% Author: Stefan Haufe, 2013
%         Tim Mullen,   2013, SCCN/INC, UCSD
% Email:  stefan.haufe@tu-berlin.de
%         tim@sccn.ucsd.edu

% This function is part of the Source Information Flow Toolbox (SIFT)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

persistent initAR;

verb = arg_extract(varargin,{'verb','VerbosityLevel'},[],2);

g = arg_define([0 1],varargin, ...
                arg_norep({'data','Data'},mandatory,[],'Data Matrix. Dimensions are [nchs x npnts].'), ...
                arg({'morder','ModelOrder','p'},10,[],'VAR Model order'), ...
                arg_subtoggle({'warmStart','WarmStart'},[], ...
                {...
                    arg({'initState','InitialState'},[],[],'Initial SCSA state. This is a structure with fields ''B'' and ''H'', which represent, respectively, the initial [M x M] ICA demixing matrix and the [M x M x morder] MVAR coefficient matrix.') ...
                },'Warm start. The previously estimated state will be used as a starting estimate for successive operations. Alternately, you may provide an non-empty initial state structure via the ''InitialState'' argument.'), ...
                arg_subtoggle({'pca','PCA'},'on',...
                {...
                    arg({'dim','NumComps'},[],[],'Number of PCs to retain. If empty, retain all components. If dim < 1 this is considered the proportion of signal variance explained to retain'), ...
                    arg_nogui({'initState','InitialState'},[],[],'PCA filter structure. This contains fields ''pcaweights'' and ''pcawinv'' containing PCA filter weights and patterns (loadings), respectively'), ...
                    arg({'iscentered','DataIsCentered'},false,[],'Set to true if data is zero mean. Otherwise mean will be removed from each channel/trial') ...
                },'Reduce dimensionality with PCA'), ...
                arg_sub({'scsa_opts','SCSA_Options'},[], ...
                {...
                	arg({'lambda','ReguParamLambda','RegularizationParam'},0.2,[0 Inf],'Regularization parameter (lambda)'), ...
                    arg({'nstart','NumRestarts'},1,[1 Inf],'Number of EM restarts'), ...
                    arg({'emiter','MaxIterations'},5,[1 Inf],'Maximum EM iterations'), ...
                    arg({'emstop','Tolerance'},1e-6,[0 Inf],'EM stopping tolerance'), ...
                },'Options for ADMM algorithm'), ...
                arg({'verb','VerbosityLevel'},verb,{int32(0) int32(1) int32(2)},'Verbosity level. 0 = no output, 1 = text, 2 = graphical') ...
                );
                
[nchs npnts ntr] = size(g.data);


% perform dimensionality reduction to largest principal components
% --------------------------------------------------------------------------------------------------
if g.pca.arg_selection
    if g.verb
        fprintf('Applying PCA...');
    end
    
    pcadim = g.pca.dim;
    if isempty(pcadim)
        pcadim = nchs; end
    
    if ~g.pca.iscentered
        % center each channel/trial
        g.data = bsxfun(@minus, g.data, sum(g.data,2)/npnts);
        % indicate that the data have been made zero-mean before running PCA and
        % SCSA. The same must be done if one wants to use pcaweights or 
        % scsafilt to recover the sources
        g.pca.iscentered = true;
    end
    
    if ~isempty(g.pca.initState)
        pcaweights  = g.pca.initState.pcaweights;
        pcawinv     = g.pca.initState.pcawinv;
        pcadim      = size(pcaweights, 1);
        g.data      = reshape(pcaweights*g.data(:,:), [pcadim, npnts, ntr]);
    else
        % Perform PCA
        % first compute average data covariance matrix across trials...
        % ...concat trials (reshape to [nchs x npnts*ntr])
        cdata = g.data(:,:);
        % ...compute (unbiased) mean covariance
        C     = cdata*cdata'/(ntr*(npnts-1));
        % compute eigendecomposition
        [V D] = eig(C);
        d     = diag(D);
        % compute variance explained
        ve = cumsum(d(end:-1:1))./sum(d);
        % determine dim reduction
        if g.pca.dim < 1
            pcadim = nnz(ve<=pcadim); end
        varexplained = ve(pcadim);
        % compute pruned principal components
        cps = length(d)-pcadim+1:length(d);
        d   = d(cps);
        V   = V(:,cps);
        pcaweights   = diag(1./sqrt(d))*V';
        pcawinv      = V*diag(sqrt(d));
        g.data       = reshape(pcaweights*cdata, [pcadim, npnts, ntr]);
        clear cdata;
    end
    
    if g.verb
        fprintf('Dimensionality reduced to %d components\n',pcadim);
    end
end

% Initialize CSA state (warm-start)
% --------------------------------------------------------------------------------------------------
if g.warmStart.arg_selection
    initAR = g.warmStart.initState;
    runCSA = false;
else
    % reset initAR
    if g.pca.arg_selection
        nvars    = pcadim;
    else
        nvars = nchs;
    end
    initAR.B = eye(nvars);
    initAR.H = zeros(nvars, nvars, g.morder);
    runCSA = true;
end

% Run SCSA-EM algorithm
% --------------------------------------------------------------------------------------------------
if g.verb
    fprintf('Running SCSA...\n'); end
scsa_opts = catstruct(g.scsa_opts,struct('start', initAR,'runCSA', runCSA,'verb',g.verb));
[scsaweights, AR, cost, srcact] = scsa_em(g.data, g.morder, scsa_opts);


% Prepare outputs
% --------------------------------------------------------------------------------------------------

% estimate noise covariance matrix
if nargout>1
    res = est_mvarResiduals(srcact,AR,zeros(1,size(srcact,1)));
    res = res(:,:);
    PE  = cov(res',1);
end

if nargout>2
    % compute filters
    scsawinv = inv(scsaweights);
    if g.pca.arg_selection
      icawinv    = pcawinv/scsaweights;
      icaweights = scsaweights*pcaweights;
    else
      icawinv    = scsawinv;
      icaweights = scsaweights;
    end

    % store outputs
    if g.pca.arg_selection
        out.pca = struct('pcaweights',  pcaweights, ...
                         'pcawinv'   ,  pcawinv,    ...
                         'varexplained',varexplained,...
                         'pcadim',   pcadim,        ...
                         'iscentered',  g.pca.iscentered);
    else
        out.pca = struct([]);
    end
    out.scsa    = struct('scsaweights',scsaweights, ...
                         'scsawinv'   ,scsawinv,    ...
                         'finalcost'  ,cost);
    out.icaact      = srcact;
    out.icaweights  = icaweights;
    out.icawinv     = icawinv;
end


