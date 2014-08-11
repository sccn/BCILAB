function GL = vect_laggrad(Xf,dX,Ff,dF,tff,dtf,lambda,probinfo)

% Make inputs for vect_cons_z
x.f  = Xf;
x.dz = dX;
f.f  = Ff;
f.dz = dF;
tf.f  = tff;
tf.dz = dtf;

% Build Constraint Jacobian
C = vect_cons_z(x,f,tf,probinfo);
JC = sparse(C.dz_location(:,1),C.dz_location(:,2),...
  C.dz,C.dz_size(1),C.dz_size(2));

% Build Objective Gradient
GO = sparse(1,probinfo.tfind,1,1,probinfo.tfind);

% Build Lagrangian Gradient
GL = GO + lambda.'*JC;