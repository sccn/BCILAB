% gl_spec - spectrum function for teh grouped L1 regularizer
% 
% Copyright(c) 2009 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt

function nm=gl_spec(ww,blks)

nn=length(blks);
nm=zeros(nn,1);

ixw=0;
for kk=1:length(blks)
  I=ixw+(1:blks(kk));
  ixw=I(end);
  nm(kk)=norm(ww(I));
end
