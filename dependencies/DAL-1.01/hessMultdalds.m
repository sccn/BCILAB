% hessMultdalds - function that computes H*x for DAL with the
%                 dual spectral norm (trace norm) regularization
%
% Copyright(c) 2009 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt

function bb = hessMultdalds(aa, A, AT, eta, Hinfo)

info=Hinfo.info;
hloss=Hinfo.hloss;
lambda=Hinfo.lambda;

nn=sum(prod(info.blks,2));

bb = hloss*aa;


ixv=0;                      % index vector before the next block begins
nblks=size(info.blks,1);    % number of blocks
% for each block...
for kk=1:nblks
    blk=info.blks(kk,:);    % current block dimensions
    bsz=blk(1)*blk(2);      % block size, numel(block)
    I=ixv+(1:bsz);          % index range into aa that indexes the block
    ixv=I(end);
    
    nsv=info.nsv(kk);
    if nsv==0
        continue;
    end
    
    if nblks>1              % note: this is a good chunk of the total compute time...
        AI = A(sparse(I,1:bsz,ones(bsz,1),nn,bsz));
    end
    
    U=info.U{kk};
    V=info.V{kk};
    S=info.S{kk};
    ss=diag(S);
    
    dd=min(blk);
    
    if nblks==1
        Q = U'*reshape(AT(aa), blk)*V;
    else
        Q = U'*reshape(AI'*aa, blk)*V;
    end
    
    % make sure that Q and ss are at least dd-sized
    if length(Q) < dd
        Q(dd,dd) = 0;  end
    if length(ss) < dd
        ss(dd) = 0; end
    
    ds=diag(Q);
    
    % get the [dd,dd]-sized principal minor of Z & W right
    Qpm = Q(1:dd,1:dd);
    ssp = ss.^2;
    M = bsxfun(@minus,ssp',ssp);
    M(1:size(M,1)+1:numel(M)) = Inf;
    Qss = bsxfun(@times,Qpm,ss');
    Z = (Qss + Qss') ./ M;
    Qss = bsxfun(@times,Qpm,ss);
    W = -(Qss + Qss')./ M';
    
    % zero-extend Z & W, if necessary
    if any(size(Z) < blk(1))
        Z(blk(1),blk(1)) = 0; end
    if any(size(W) < blk(2))
        W(blk(2),blk(2)) = 0; end
    
    % this part operates on the edges of Q (beyond the principal minor)
    if blk(1)<blk(2)
        W(dd+1:end,1:dd)=Q(1:dd,dd+1:end)'/diag(ss);
        W(1:dd,dd+1:end)=-W(dd+1:end,1:dd)';
    else
        Z(dd+1:end,1:dd)=Q(dd+1:end,1:dd)/diag(ss);
        Z(1:dd,dd+1:end)=-Z(dd+1:end,1:dd)';
    end
    
    ssth = ss(1:nsv) - lambda;
    yy=U*(Z(:,1:nsv)*spdiag(ssth)*V(:,1:nsv)')...
        +U(:,1:nsv)*spdiag(ds(1:nsv))*V(:,1:nsv)'...
        +(U(:,1:nsv)*spdiag(ssth)*W(:,1:nsv)')*V';
    
    if nblks==1
        bb = bb + eta*(A(yy(:)));
    else
        bb = bb + eta*(AI*yy(:));
    end
end

B=Hinfo.B;
if ~isempty(B)
    bb = bb + eta*(B'*aa)*B;
end
