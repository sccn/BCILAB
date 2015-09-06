function [xi,yi,curvatures] = maxr(a) 
%MAXR   Find interpolated maximizer(s) and max value(s) 
%       for (each column of) a.
%
%               [xi,yi,curvatures] = maxr(a) 
%
%       Calls max() followed by qint() for quadratic interpolation.
%  
   [m,n] = size(a); 
   if m==1, a=a'; t=m; m=n; n=t; end; 
   [y,x] = max(a); 
   xi=x;    % vector of maximizer locations, one per col of a
   yi=y;    % vector of maximum values, one per column of a
   if nargout>2, curvatures = zeros(1,n); end
   for j=1:n,   % loop over columns
     if x(j)>1, % only support interior maxima
       if x(j)<m, 
         [xdelta,yij,cj] = qint(a(x(j)-1,j),y(j),a(x(j)+1,j)); 
         xi(j) = x(j) + xdelta;
         if nargout>2, curvatures(j) = cj; end
         if (nargout>1), yi(j) = yij; end
       end; 
     end; 
   end;