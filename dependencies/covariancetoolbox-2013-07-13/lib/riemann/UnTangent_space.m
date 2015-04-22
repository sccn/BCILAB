function COV = UnTangent_space(T,C)
NTrial = size(T,2);
N_elec = (sqrt(1+8*size(T,1))-1)/2;
COV = zeros(N_elec,N_elec,NTrial);

if nargin<2
    C = riemann_mean(COV);
end

index = reshape(triu(ones(N_elec)),N_elec*N_elec,1)==0;

Out = zeros(N_elec*N_elec,NTrial);

Out(not(index),:) = T;
P = C^0.5;
for i=1:NTrial
  tmp = reshape(Out(:,i),N_elec,N_elec,[]);
  tmp =  diag(diag(tmp))+triu(tmp,1)/sqrt(2) + triu(tmp,1)'/sqrt(2);
  tmp = P*tmp*P;
  COV(:,:,i) = RiemannExpMap(C,tmp);
end
