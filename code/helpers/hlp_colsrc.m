function srcmat_out = hlp_colsrc(srcmat,roiVertices,rule)
% collapse source matrix across groups of vertices using some rule.
% rule may be one of {'mean','sum','max','maxmag','median'}
% Author: Tim Mullen, SCCN/INC/UCSD, 2013

% arg_define([0 Inf],varargin, ...
%            arg_norep('srcmat',[],[],'Source matrix. [num_sources x dim1 x ... x dim N]'), ...
%            arg_norep('roiVertices',[],[],'ROI Vertices. Cell array of vertices for each group to collapse over.'), ...
%            arg({'rule','CollapseRule'},'mean',{'mean','sum','max','maxmag','median'},'Collapse rule'));

nrois = length(roiVertices);
szact = size(srcmat);

% initialization
srcmat_out = zeros([nrois,szact(2:end)]);

switch rule
    case 'mean'
        % average CSD for each ROI
        for k=1:nrois
            x = roiVertices{k};
            srcmat_out(k,:,:) = mean(srcmat(x,:,:),1);
        end
    case 'sum'
        % sum CSD for each ROI
        for k=1:nrois
            x = roiVertices{k};
            srcmat_out(k,:,:) = sum(srcmat(x,:,:),1);
        end
    case 'max'
        % get max CSD for each ROI
        for k=1:nrois
            x = roiVertices{k};
            srcmat_out(k,:,:) = max(srcmat(x,:,:),[],1);
        end
    case 'maxmag'
        % get CSD with maximum magnitude for each ROI
        for k=1:nrois
            x = roiVertices{k};
            Q = srcmat(x,:,:);
            [dummy idx] = max(abs(Q),[],1);
            srcmat_out(k,:,:) = Q((0:numel(idx)-1).*size(Q,1) + idx);
        end
    case 'median'
        % get median CSD for each ROI
        for k=1:nrois
            x = roiVertices{k};
            if ~isempty(x)
                srcmat_out(k,:,:) = median(srcmat(x,:,:),1); end
        end
    otherwise
        error('Unsupported collapse rule: %s',rule);
end


% OLD METHOD (simulation shows loop is faster)
% switch colRoiCsd
%     case 'mean'
%         % average CSD for each ROI
%         roifunc_csd = @(x) mean(srcmat(x,:),1);
%         roifunc_inv = @(x) mean(srcweights(x,:),1);
%     case 'sum'
%         % sum CSD for each ROIs
%         roifunc_csd = @(x) sum(srcmat(x,:),1);
%         roifunc_inv = @(x) sum(srcweights(x,:),1);
%     case 'max'
%         % get max CSD for each ROI
%         roifunc_csd = @(x) max(srcmat(x,:),[],1);
%         roifunc_inv = @(x) max(srcweights(x,:),[],1);
% end
%
% if ~strcmp(colRoiCsd,'none')
%     if verb
%         fprintf('Computing %s CSD for each ROI \n', colRoiCsd);
%     end
%     roicsd = cellfun(roifunc_csd,state.roiVertices,'UniformOutput',false);
%     invcsd = cellfun(roifunc_inv,state.roiVertices,'UniformOutput',false);
%     srcmat         = cell2mat(roicsd');
%     srcweights     = cell2mat(invcsd');
% end
