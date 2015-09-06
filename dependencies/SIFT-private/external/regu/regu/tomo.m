function [A,b,x] = tomo(N,f)
%TOMO  Create a 2D tomography test problem
%
% [A,b,x] = tomo(N,f);
%
% This function creates a simple two-dimensional tomography test problem.
% A 2D domain [0,N] x [0,N] is divided into N^2 cells of unit size, and a
% total of round(f*N^2) rays in random directions penetrate this domain.
% The default value is f = 1.
%
% Each cell is assigned a value (stored in the vector x), and for each
% ray the corresponding element in the right-hand side b is the line
% integral along the ray, i.e.
%    sum_{cells in ray}  x_{cell j} * length_{cell j}
% where length_{cell j} is the length of the ray in the j-th cell.
%
% The matrix A is sparse, and each row (corresponding to a ray) holds
% the value length_{cell j} in the j-th position.  Hence:
%    b = A*x .
% Once a solution x_reg has been computed, it can be visualized by means
% of imagesc(reshape(x_reg,N,N)).
%
% The exact solution, reshape(x,N,N), is identical to the exact image in
% the function blur.

% Per Christian Hansen, IMM, DTU; March 22, 2007.

% Default values.
if nargin==1, f = 1; end

% Allocate space for sparse matrix.
A = spalloc(N,N,2*N);

% Prepare for loop over rows (i.e., rays).
x = (0:N)'; y = x;

% Loop over all rays.
for i=1:round(f*N^2)
    % Random line through the domain (y = a*x + b).
    x0 = N*rand; x1 = N*rand;
    y0 = N*rand; y1 = N*rand;
    a = (y1-y0)/(x1-x0);
    b = y0 - a*x0;

    % All points where the line intersects the grid.
    yp = a*x + b;
    xp = (y - b)/a;
    xp = [x;xp];
    yp = [yp;y];

    % Sort them, and skip those outside the box.
    [xp,I] = sort(xp);
    yp = yp(I);
    I = find(xp >= 0 & xp <= N & yp >= 0 & yp <= N);
    xp = xp(I);
    yp = yp(I);
    % Skip double points.
    I = find(diff(xp)==0);
    xp(I) = []; yp(I) = [];

    % Lengths within cells.
    d = sqrt( diff(xp).^2 + diff(yp).^2 );

    % Column indices to cells penetrated by the ray.
    xm = 0.5*(xp(1:end-1)+xp(2:end));
    ym = 0.5*(yp(1:end-1)+yp(2:end));
    j = ( floor(xm) )*N + floor(ym) + 1;
    
    % Store lengths in i-th row.
    A(i,j) = d';
    
end

if nargout>1

  % Start with an image of all zeros. 
  x = zeros(N,N); 
  N2 = round(N/2); 
  N3= round(N/3); 
  N6 = round(N/6); 
  N12 = round(N/12); 
 
  % Add a large ellipse. 
  T = zeros(N6,N3); 
  for i=1:N6
    for j=1:N3 
      if ( (i/N6)^2 + (j/N3)^2 < 1 ), T(i,j) = 1; end 
    end
  end
  T = [fliplr(T),T]; 
  T = [flipud(T);T]; 
  x(2+(1:2*N6),N3-1+(1:2*N3)) = T; 
 
  % Add a smaller ellipse. 
  T = zeros(N6,N3); 
  for i=1:N6
    for j=1:N3 
      if ( (i/N6)^2 + (j/N3)^2 < 0.6 ), T(i,j) = 1; end 
    end
  end
  T = [fliplr(T),T];
  T = [flipud(T);T];
  x(N6+(1:2*N6),N3-1+(1:2*N3)) = x(N6+(1:2*N6),N3-1+(1:2*N3)) + 2*T;
  % Correct for overlap.
  f = find(x==3);
  x(f) = 2*ones(size(f));

  % Add a triangle.
  T = triu(ones(N3,N3));
  [mT,nT] = size(T);
  x(N3+N12+(1:nT),1+(1:mT)) = 3*T;

  % Add a cross.
  T = zeros(2*N6+1,2*N6+1);
  [mT,nT] = size(T);
  T(N6+1,1:nT) = ones(1,nT);
  T(1:mT,N6+1) = ones(mT,1);
  x(N2+N12+(1:mT),N2+(1:nT)) = 4*T;

  % Make sure x is N-times-N, and stack the columns of x.
  x = reshape(x(1:N,1:N),N^2,1);
  
end

if nargout==3, b = A*x; end