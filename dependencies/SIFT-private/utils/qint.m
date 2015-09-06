function [p,y,a] = qint(ym1,y0,yp1) 
%QINT   Quadratic interpolation of 3 uniformly spaced samples
%
%               [p,y,a] = qint(ym1,y0,yp1) 
%
%       returns extremum-location p, height y, and curvature a
%       of a parabolic fit through three points. 
%       The parabola is given by y(x) = a*(x-p)^2+b, 
%       where y(-1)=ym1, y(0)=y0, y(1)=yp1. 

   p = (yp1 - ym1)/(2*(2*y0 - yp1 - ym1)); 
   if nargout>1
     y = y0 - 0.25*(ym1-yp1)*p;
   end;
   if nargout>2
     a = 0.5*(ym1 - 2*y0 + yp1);
   end;