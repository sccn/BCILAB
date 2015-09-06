% gl_spec - spectrum function for the grouped L1 regularizer
% 
% Copyright(c) 2009 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt

function nm=gl_spec(ww,info)



% old code (slow)
% blks = info.blks;
% nn=length(blks);
% nm=zeros(nn,1);
% ixw=0;
% for kk=1:length(blks)
%   I=ixw+(1:blks(kk));
%   ixw=I(end);
%   nm(kk)=norm(ww(I));
% end

nn=length(info.blks);
nm=zeros(nn,1);
for k=1:length(info.blkgrp)
    kk = info.blkvec{k};
    I = info.blkgrp{k};    
    nm(kk)=sqrt(sum(ww(I)'.^2));
end
