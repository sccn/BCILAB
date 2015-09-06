% gl_softth - soft threshold function for grouped L1 regularization
%
% Copyright(c) 2009 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt

function [vv,ss] = gl_softth(vv,lambda,info)

% init ss
ss=zeros(length(info.blks),1);
% for each block group...
for k=1:length(info.blkgrp)
    kk = info.blkvec{k};
    I = info.blkgrp{k};
    % compute the row-wise norms of the k'th shard of vv
    ssn = sqrt(sum(vv(I)'.^2));
    % threshold it
    ssk = max(ssn-lambda,0);
    % scale vv
    vv(I) = bsxfun(@times,vv(I),(ssk./ssn)');
    % update ss(kk)
    ss(kk) = ssk;
end
vv=vv(:);
ss=ss(:);
