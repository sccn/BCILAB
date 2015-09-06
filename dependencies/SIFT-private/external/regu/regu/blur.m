function [A,b,x] = blur(N,band,sigma)
%BLUR Test problem: digital image deblurring.
%
% function [A,b,x] = blur(N,band,sigma)
%
% The matrix A is an N*N-by-N*N symmetric, doubly block Toeplitz matrix that
% models blurring of an N-by-N image by a Gaussian point spread function.
% It is stored in sparse matrix format.
%
% In each Toeplitz block, only matrix elements within a distance band-1
% from the diagonal are nonzero (i.e., band is the half-bandwidth).
% If band is not specified, band = 3 is used.
%
% The parameter sigma controls the width of the Gaussian point spread
% function and thus the amount of smoothing (the larger the sigma, the wider
% the function and the more ill posed the problem).  If sigma is not
% specified, sigma = 0.7 is used.
%
% The vector x is a columnwise stacked version of a simple test image, while
% b holds a columnwise stacked version of the blurrred image; i.e, b = A*x.

% Per Christian Hansen, IMM, 11/11/97.

% Initialization.
if (nargin < 2), band = 3; end
band = min(band,N);
if (nargin < 3), sigma = 0.7; end

% Construct the matrix as a Kronecker product.
z = [exp(-((0:band-1).^2)/(2*sigma^2)),zeros(1,N-band)];
A = toeplitz(z);
A = sparse(A);
A = (1/(2*pi*sigma^2))*kron(A,A);

% Generate x and b, if required.
if (nargout > 1)

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

  % Compute the blurred image.
  b = A*x;

end