function H = randmvar(M,P,k)

H=zeros(M*M,P);
I=randperm(M*M); I=I(1:k);
H(I,:)=0.1*randn(k,P);
H=reshape(H,[M,M,P]);
for ii=1:P
  H(:,:,ii)=H(:,:,ii)+0.01*diag(randn(M,1));
end

