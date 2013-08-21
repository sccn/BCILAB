% ds_softth - soft threshold function for DS regularization
%
% Copyright(c) 2009 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt

function [vv,ss,info]=ds_softth(vv,lambda,info)
% soft-threshold operation for a dual-spectral norm (also returning some misc stuff)
% returns the thresholded matrix, as well as the regularizer spectrum

ss=zeros(sum(min(info.blks,[],2)),1);

ixs=0;
ixv=0;

% block-shapedness...
for kk=1:size(info.blks,1)
    blk=info.blks(kk,:);
    I=ixv+(1:blk(1)*blk(2));
    ixv=I(end);
    J=ixs+(1:min(blk));
    ixs=J(end);
    vvk = reshape(vv(I),blk);
    
    % apparently, this thresholds the (matrix-shaped gradient) in a singular-value thresholding
    % manner
    
    % note: this is part of the problem with the Hessian: if min(size(vvk))>500, the lansvd is used
    %       which returns only those SV's that are larger than lambda... so U,S,V are tall, small, and fat, respectively...
    [U,S,V] = svdmaj(vvk, lambda, info.nsv(kk)+1,5); 
    
    % get the diagonal SV's
    dd=diag(S);
    
    % find the singular values beyond the threshold
    K=find(dd>lambda);
    % subtract lambda
    ssk=dd(K)-lambda;
    % write'em back (note: the rest is already initialized to zero)
    ss(J(1:length(ssk)))=ssk;
    %reconstruct only the non-thresholded parts of VV
    vvk=U(:,K)*diag(ssk)*V(:,K)';
    vv(I)=vvk(:);
    
    % this stuff is later re-used for other parts (e.g. hessian computation)
    info.nsv(kk)=length(ssk);
    info.U{kk}=U;
    info.S{kk}=S;
    info.V{kk}=V;
end
