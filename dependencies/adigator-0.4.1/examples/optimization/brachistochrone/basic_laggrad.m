function GL = basic_laggrad(z,lambda,probinfo)
% This function builds the lagrangian gradient

Z.f = z;
Z.dz = ones(length(z),1);
GO = sparse(1,probinfo.tfind,1,1,probinfo.tfind); % Obj Grad

[~,Ceq] = basic_cons_z(Z,probinfo);
JCeq = sparse(Ceq.dz_location(:,1),Ceq.dz_location(:,2),Ceq.dz,...
  Ceq.dz_size(1),Ceq.dz_size(2)); % Cons Jac
GL = GO + lambda.'*JCeq; % Lag Grad