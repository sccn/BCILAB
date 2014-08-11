function [tf,G] = basic_objwrap(z,probinfo)

tf = z(probinfo.tfind);
if nargout > 1
  G = sparse(1,probinfo.tfind,1,1,probinfo.tfind);
end