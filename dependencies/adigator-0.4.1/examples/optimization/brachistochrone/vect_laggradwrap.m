function [H,G] = vect_laggradwrap(z,lambda,probinfo)


% Create Vectorized Inputs
Xvec.f = z(probinfo.xind);
Xvec.dY = ones(size(Xvec.f));
Uvec.f  = z(probinfo.uind);
Uvec.dY = ones(size(Uvec.f));

% Get Vectorized Derivs
Fvec = dynamics_YY(Xvec,Uvec);



% Create Non-Vectorized Inputs
X.f  = Xvec.f;
X.dz = Xvec.dY(:);
dX   = X.dz;


F.f = Fvec.f;
F.dz = Fvec.dY(probinfo.dFind);
dF.f = F.dz;
dF.dz = Fvec.dYdY(probinfo.dF2ind);

tf.f  = z(probinfo.tfind);
tf.dz = 1;
dtf   = tf.dz;


GL = vect_laggrad_z(X,dX,F,dF,tf,dtf,lambda.eqnonlin,probinfo);

H = sparse(GL.dz_location(:,1),GL.dz_location(:,2),GL.dz,...
    GL.dz_size(1),GL.dz_size(2));
  
G = GL.f;