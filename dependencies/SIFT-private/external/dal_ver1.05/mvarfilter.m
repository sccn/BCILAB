% mvarfilter - Multivariate AR filter
%
% Example:
%  H=randmvar(20,3,10);
%  Z=mvarfilter(H0, 2/pi*log(tan(pi*rand(N,M)/2)))';
%
% Copyright(c) 2011 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt

function Y=mvarfilter(A, X)

[M1,M2,P]=size(A);
[N,M]=size(X);

if M1~=M2 || M1~=M
  error('Matrix sizes don''t match.');
end

X=X';
Y=zeros(M,N);
zz=zeros(M*P,1);
for ii=1:N
  zz1=A(:,:)*zz+X(:,ii);
  zz=[zz1; zz(1:end-M)];
  Y(:,ii)=zz1;
end
Y=Y';