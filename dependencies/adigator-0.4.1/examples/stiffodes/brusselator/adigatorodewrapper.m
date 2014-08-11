function [J f] = adigatorodewrapper(t,y,N)
Y.f = y;
Y.dy = ones(size(y));
dYdt = brussderiv(t,Y,N);
f = dYdt.f;
J = sparse(dYdt.dy_location(:,1),dYdt.dy_location(:,2),...
  dYdt.dy,dYdt.dy_size(1),dYdt.dy_size(2));
end