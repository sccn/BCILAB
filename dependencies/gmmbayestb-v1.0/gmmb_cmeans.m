% GMMB_CMEANS  simple c-means clustering
%
% T = CMEANS(data, nclust, count)
% [T, CLUST] = CMEANS(...)
%
% data    input data, N x D matrix
% nclust  number of clusters
% count   number of iterations
%
% T       output, data labels, 1 x N vector
% CL      output, cluster centers, nclust x D matrix
%
% Author: Jarmo Ilonen
% Editor: Pekka Paalanen
%
% $Name:  $ $Id: gmmb_cmeans.m,v 1.1 2004/11/02 08:32:22 paalanen Exp $

function [pclass, clust]=gmmb_cmeans(pdata,nclust,count);

rp = randperm(size(pdata,1));
clust = pdata(rp(1:nclust),:);

for kierros=1:count,
	% compute squared distance from every point to every cluster center.
	for i=1:nclust,
		vd = pdata - repmat(clust(i,:),size(pdata,1),1);
		cet(:,i) = sum(abs(vd).^2, 2);
	end;

	% compute new cluster centers
	[a, pclass]=min(cet');

	for i=1:nclust,
		clust(i,:) = mean( pdata(find(pclass==i), :) );
	end;
end;
