function X = tensorsum(A,B)
%TENSORSUM
%
% Modification of Laurent Sorber's KRON

[I J] = size(A);
[K L] = size(B);
if ~issparse(A) && ~issparse(B)
    A = reshape(A,[1 I 1 J]);
    B = reshape(B,[K 1 L 1]);
    X = reshape(bsxfun(@plus,A,B),[I*K J*L]);
else
    
    [ia,ja,sa] = find(A); ia=ia(:); ja=ja(:); sa=sa(:);
    [ib,jb,sb] = find(B); ib=ib(:); jb=jb(:); sb=sb(:);
    
    ix = bsxfun(@plus,K*(ia-1).',ib);
    jx = bsxfun(@plus,L*(ja-1).',jb);
    

    if ~(isdouble(sa)||islogical(sa)); sa=double(sa); end 
    if ~(isdouble(sb)||islogical(sb)); sb=double(sb); end 
     
    X = sparse(ix,jx,bsxfun(@plus,sb,sa.'),I*K,J*L);

end