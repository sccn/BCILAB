function [C,Ceq,JC,JCeq] = vect_conswrap(z,probinfo)

C = [];
JC = [];
if nargout > 2
% Create Vectorized Inputs
Xvec.f = z(probinfo.xind);
Xvec.dY = ones(size(Xvec.f));
Uvec.f  = z(probinfo.uind);
Uvec.dY = ones(size(Uvec.f));

% Get Vectorized Derivs
Fvec = dynamics_Y(Xvec,Uvec);

% Create Non-Vectorized Inputs
X.f  = Xvec.f;
X.dz = Xvec.dY(:);
F.f = Fvec.f;
F.dz = Fvec.dY(probinfo.dFind);
tf.f = z(probinfo.tfind);
tf.dz = 1;

Ceq = vect_cons_z(X,F,tf,probinfo);
JCeq = sparse(Ceq.dz_location(:,1),Ceq.dz_location(:,2),Ceq.dz,...
  Ceq.dz_size(1),Ceq.dz_size(2)).';
Ceq = Ceq.f;

else
  X = z(probinfo.xind);
  U = z(probinfo.uind);
  tf = z(probinfo.tfind);
  F = dynamics(X,U);
  Ceq = vect_cons(X,F,tf,probinfo);
end