function y = normpdfg(x,scale,shape,center,max)
% y = normpdfg(x,scale,shape,center,max)
%
% Evalute the generalized gaussian pdf at point(s) x
% the pdf is defined by the given scale (2*variance) and shape and mean
% (center)
%
% Note the variance of the distribution is scale/2
%
% As shape --> Inf we get a (step function) uniform distribution over [center-scale, center+scale]
% As shape --> 0 we get a delta function
% Beta = 1 is the Laplace 
% Beta = 2 is gaussian
%
% (C) Tim Mullen, 2011. SCCN/INC UCSD

A = shape./(2.*scale.*gamma(1/shape));
C = max./A;

y = C.*A.*exp(-((abs(x-center)/scale).^shape));