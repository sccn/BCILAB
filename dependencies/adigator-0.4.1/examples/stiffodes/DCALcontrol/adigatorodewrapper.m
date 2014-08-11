function [J f] = adigatorodewrapper(t,x)
X.f = x;
X.dx = ones(size(x));
dXdt = TwoLinkSys_x(t,X);
f = dXdt.f;
J = sparse(dXdt.dx_location(:,1),dXdt.dx_location(:,2),...
  dXdt.dx,dXdt.dx_size(1),dXdt.dx_size(2));
end