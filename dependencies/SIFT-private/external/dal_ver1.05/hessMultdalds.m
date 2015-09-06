% hessMultdalds - function that computes H*x for DAL with the
%                 dual spectral norm (trace norm) regularization
%
% Copyright(c) 2009 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt

function bb = hessMultdalds(aa, A, eta, Hinfo)

info=Hinfo.info;
hloss=Hinfo.hloss;
lambda=Hinfo.lambda;

nn=sum(prod(info.blks,2));

bb = hloss*aa;

ixv=0;
nblks=size(info.blks,1);
for kk=1:nblks
  blk=info.blks(kk,:);
  bsz=blk(1)*blk(2);
  I=ixv+(1:bsz);
  ixv=I(end);

  nsv=info.nsv(kk);
  if nsv==0
    continue;
  end

  if nblks>1
    AI=A.slice(I);
  end
  
  U=info.U{kk};
  V=info.V{kk};
  S=info.S{kk};
  ss=diag(S);
  
  Z=zeros(blk(1));
  W=zeros(blk(2));
  
  dd=min(blk);
  if nblks==1
    Q = U'*reshape(A.Ttimes(aa), blk)*V;
  else
    Q = U'*reshape(AI'*aa, blk)*V;
  end
  ds=diag(Q);
  for ii=1:dd-1, for jj=ii+1:dd
      Z(ii,jj)=(ss(jj)*Q(ii,jj)+ss(ii)*Q(jj,ii))/(ss(jj)^2-ss(ii)^2);
      W(jj,ii)=(ss(ii)*Q(ii,jj)+ss(jj)*Q(jj,ii))/(ss(ii)^2-ss(jj)^2);
      Z(jj,ii)=-Z(ii,jj);
      W(ii,jj)=-W(jj,ii);
    end
  end
  if blk(1)<blk(2)
    W(dd+1:end,1:dd)=Q(1:dd,dd+1:end)'/diag(ss);
    W(1:dd,dd+1:end)=-W(dd+1:end,1:dd)';
  else
    Z(dd+1:end,1:dd)=Q(dd+1:end,1:dd)/diag(ss);
    Z(1:dd,dd+1:end)=-Z(dd+1:end,1:dd)';
  end
  
  ssth=ss(1:nsv)-lambda;
  
  yy=U*(Z(:,1:nsv)*spdiag(ssth)*V(:,1:nsv)')...
      +U(:,1:nsv)*spdiag(ds(1:nsv))*V(:,1:nsv)'...
      +(U(:,1:nsv)*spdiag(ssth)*W(:,1:nsv)')*V';

  if nblks==1
    bb = bb + eta(1)*(A.times(yy(:)));
  else
    bb = bb + eta(1)*(AI*yy(:));
  end
end

B=Hinfo.B;
if ~isempty(B)
  bb = bb + eta(2)*(B'*aa)*B;
end
