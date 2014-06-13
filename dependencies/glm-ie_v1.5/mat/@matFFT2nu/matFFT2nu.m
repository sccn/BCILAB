% Non-uniformly spaced 2d FFT matrix: constructor
%
% We approximate the matrix A by P*F*diag(d) where d is a (real) scaling vector,
% P is a highly sparse (complex) matrix containting interpolation weights and 
% F is the Fourier transform matrix twice as large as the original image due 
% to overgridding.
%
% The implementation is based on Jeffrey Fessler's code: 
%                              http://www.eecs.umich.edu/~fessler/irt/irt/nufft/
%
% A = matFFT2nu(imsz,k,complex,kbwsz,ogf)
%
%  imsz    [1,2]: size of the image to transform
%  k       [m,1]: normalized k-space coordinates i.e. real&imag from[-1,1]/2
%  complex [1,2]: complex number treatment A:IN^n->OUT^m,          default [1,1]
%                 first number specifies output field OUT^m and the second the
%                 input field IN^m of the matrix A
%                 1 -  real(z) + 1i*imag(z)    IC:   native complex numbers
%                 2 - [real(z); imag(z)]       IR^2: stacked real numbers
%                 3    real(z)                 IR:   real part only
%  kbwsz   [1,2]: Kaiser-Bessel interpolation window size,         default [3,3]
%  ogf     [1,1]: overgridding factor                                [default 2]
% 
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 December 9

function  A = matFFT2nu(imsz,k,complex,kbwsz,ogf)

  if numel(imsz)==1, imsz = imsz*[1,1]; end, A.imsz = imsz;         % image size
  sz = [numel(k),prod(imsz)];                                      % matrix size
  if nargin<3, complex = 1; end
  if numel(complex)==1, complex = complex*[1,1]; end
  if nargin>3, A.kbwsz = max(fix(kbwsz(:)'),3); else A.kbwsz = [3,3]; end
  if nargin<5, ogf = 2; end                   % default overgridding factor of 2
  A.fftsz = fix(ogf*imsz);
  [A.d, A.P] = nufft_init(k, A.imsz, A.fftsz, A.kbwsz);   % nufft precomputation
  A = class(A,mfilename,mat(sz,complex));


% Initialize structure for NUFFT using KB interpolator,
% particularly the interpolation matrix in sparse format.
% caution: this routine can require a lot of memory!
% INPUT
%	k	      k-space coordinates 
%	imsz 		image size
% OUTPUT
%	P		[M, fftsz]	sparse interpolation matrix
%	d		[(imsz)]		scaling factors
%
% Like fft(), the NUFFT expects the signals to be x(0,0), ...
%
% (c) by Jeff Fessler, The University of Michigan

function [d,P] = nufft_init(k, imsz, fftsz, kbwsz)

  %% 1) (real) scaling factors d: outer product of vectors
  [d1, alpha{1}] = nufft_scale(imsz(1), kbwsz(1), fftsz(1));
  [d2, alpha{2}] = nufft_scale(imsz(2), kbwsz(2), fftsz(2));
  d = d1*d2';
  %% 2) interpolation coefficient  matrix P: huge sparse matrix
  om = [-imag(k(:)), real(k(:))];             % separate real and imaginary part
  for id=1:2
    N = imsz(id); J = kbwsz(id); K = fftsz(id); al = alpha{id};

    % compute matrix T = [C' S S' C]\inv used in NUFFT, independent of frequency
    L = length(al) - 1;	cssc = zeros(J,J);
    [j1 j2] = ndgrid(1:J, 1:J);
    for l1 = -L:L
      for l2 = -L:L
        alf1 = al(abs(l1)+1);
        if l1 < 0, alf1 = conj(alf1); end
        alf2 = al(abs(l2)+1);
        if l2 < 0, alf2 = conj(alf2); end

        tmp = j2 - j1 + l1 - l2;
        tmp = nufft_sinc(tmp/(K/N));
        cssc = cssc + alf1 * conj(alf2) * tmp;
      end
    end
    T = inv(cssc);

    % compute NUFFT r vector
    M = length(om(:,id));
    % offset for NUFFT
    if mod(J,2)	% odd J
      koff = round(K*om(:,id)) - (J+1)/2;
    else		% even J
      koff = floor(K*om(:,id)) - J/2;
    end
    dk = K*om(:,id) - koff;
    arg = outer_sum(-(1:J)', dk');
    r = zeros(J,M);
    for l1 = -L:L
      alf = al(abs(l1)+1);
      if l1<0, alf = conj(alf); end
      r1 = nufft_sinc( (arg+l1)/(K/N) );
      r = r + alf * r1;
    end

    ud{id} = exp(2i*pi*(N-1)/(2*K) * arg) .* (T*r);	% linear phase
    % indices into oversampled FFT components
    kd{id} = mod(outer_sum((1:J)', koff'), K) + 1;
    if id > 1	% trick: pre-convert these indices into offsets!
      kd{id} = (kd{id}-1) * prod(fftsz(1:(id-1)));
    end
  end

  % build sparse matrix with entries per frequency point
  kk = kd{1};	uu = ud{1};
  kk = block_outer_sum(kk, kd{2});	% outer sum of indices
  kk = reshape(kk, prod(kbwsz), M);
  uu = block_outer_prod(uu, ud{2});	% outer product of coefficients
  uu = reshape(uu, prod(kbwsz), M);
  mm = (1:M); mm = mm(ones(prod(kbwsz),1),:);                 % pre-do transpose
  P = sparse(mm(:), kk(:), conj(uu(:)), M, prod(fftsz));	       % sparse matrix


% scaling factors for 1D NUFFT (from Fourier series coefficients)
function [d,al] = nufft_scale(N, J, K)                                % real out
  % Compute scaling factors using the "do no harm" strategy.
  ssn = 1 ./ kaiser_bessel_ft( J, (1/2-N/2:N/2-1/2)'/K );
  % use regression to match NUFFT with best Kaiser scalings
  nn = -(N-1)/2:(N-1)/2;  L = min(ceil(N/3),13);
  X = cos(2*pi*nn'*(0:L)/K);
  al = (X \ ssn)';	al(2:end) = al(2:end)/2; % LS regress
  L = length(al)-1;
  d = exp( 1i*2*pi/K * nn'*(-L:L) ) * [fliplr(al(2:end)),al]';
  d = real(d);

% generalized Kaiser-Bessel function for x in support [-J/2,J/2]
function [alpha,kb] = kaiser_bessel(J,x)                      % real in real out
  o2 = [1 1];                        % hardwired, because it was nearly the best
  zn = [2.5 2.27 2.31 2.34 2.32*o2 2.35 2.34*o2 2.35 2.34 2.35*[1 1 1] 2.33]';
  alpha = J*zn(J==(2:16)');
  if nargout>1 && nargin>1
    ii = abs(x) < J/2;
    f = sqrt(1 - (x(ii)/(J/2)).^2);
    kb = zeros(size(x));
    kb(ii) = real( besseli(0, alpha*f) / besseli(0,alpha) );
  end

% Fourier transform of generalized Kaiser-Bessel function
function y = kaiser_bessel_ft(J, x)                           % real in real out
  alpha = kaiser_bessel(J);
  z = sqrt( (2*pi*(J/2)*x).^2 - alpha^2 );
  y = real( sqrt(2*pi)*J/2 * besselj(1/2,z) ./ besseli(0,alpha) ./ sqrt(z) );

% sinc
function y = nufft_sinc(x)
  iz = find(x == 0);  x(iz) = 1;  y = sin(pi*x) ./ (pi*x);  y(iz) = 1;

function y = block_outer_sum(x1, x2)
  x1 = reshape(x1, size(x1,1), 1, []);
  x2 = reshape(x2, 1, size(x2,1), []);
  y = x1(:,ones(size(x2,2),1),:) + x2(ones(size(x1,1),1),:,:);

function y = block_outer_prod(x1, x2)
  x1 = reshape(x1, size(x1,1), 1, []);
  x2 = reshape(x2, 1, size(x2,1), []);
  y = x1(:,ones(size(x2,2),1),:) .* x2(ones(size(x1,1),1),:,:);

function s = outer_sum(x,y)
	x = x(:); y = y(:)'; s = x(:, ones(numel(y),1)) + y(ones(numel(x),1),:);