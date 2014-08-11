function [C,Ceq,JC,JCeq] = basic_conswrap(z,probinfo)


if nargout == 2
  [C,Ceq] = basic_cons(z,probinfo);
else
  Z.f = z;
  Z.dz = ones(size(z));
  [~,Ceq] = basic_cons_z(Z,probinfo);
  JC = [];
  C = [];
  JCeq = sparse(Ceq.dz_location(:,1),Ceq.dz_location(:,2),Ceq.dz,...
    Ceq.dz_size(1),Ceq.dz_size(2)).';
  Ceq = Ceq.f;
end