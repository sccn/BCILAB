function [H,G] = basic_laggradwrap(z,lambda,probinfo)

Z.f = z;
Z.dz = ones(size(z));
GL = basic_laggrad_z(Z,lambda.eqnonlin,probinfo);
H = sparse(GL.dz_location(:,1),GL.dz_location(:,2),GL.dz,...
    GL.dz_size(1),GL.dz_size(2));
  
  G = GL.f;