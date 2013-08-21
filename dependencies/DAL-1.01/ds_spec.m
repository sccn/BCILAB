function nm=ds_spec(ww,blks)
% spectrum function of the regularizer
% note that this function is ALSO block-aware, so it can be sped up...

nn=sum(min(blks,[],2));

nm=zeros(nn,1);

ix0=0;
ixw=0;
for kk=1:size(blks,1)
  blk=blks(kk,:);
  I=ixw+(1:blk(1)*blk(2));
  ixw=I(end);
  J=ix0+(1:min(blk));
  ix0=J(end);
  
  % just get the SVs here - simple!
  nm(J)=svd(reshape(ww(I),blk));
end
