% Non-uniformly spaced 2d/3d FFT matrix: constructor
%
% A = matFFTnu(sz,k,type)
%
% sz:   size of the array; either 2d or 3d
% k:    normalized k-space coordinates i.e. in [-1,1]/2
%       In case of 2d input, either a [m,1] array where real and complex part
%       contain first and second cooordinate or an [m,2] array containing the
%       locations in k-space. For 3 input, simply an [m,3] array.
%
% type: string indicating the nufft interpolation type,
%       we support 'dense/full', 'nearest', 'linear', 'cubic',
%                  'fessler' and 'ferrara'                      default 'linear'
%  'dense' or 'full'
%    We return the exact dense matrix. This requires quite some memory, hence it
%    is only possible for small sz.
%  'nearest', 'linear', 'cubic'
%    We use nearest neighbor, linear or cubic interpolation after oversampling.
%  'fessler'
%    Code inspired by Jeff Fessler's nufft library [2], implementing [1].
%    The package is very general, we use Kaiser-Bessel windows of size 3.
%    The matrix is approximated by P*F*diag(d) where d is a (real) scaling
%    vector, P is a highly sparse (complex) matrix containting interpolation
%    weights and F is the Fourier transform matrix twice as large as the
%    original image due to overgridding.
%  'ferrara'
%    Code inspired by Matthew Ferrara's nufft package [4], implementing [3].
%    The code uses Gaussian windows which speeds up computations and makes sure
%    that the nufft and fftn coincide if the points are placed on a grid.
%    We use overgridding with a factor of 2.
%
% [1] Jeff Fessler and Brad P. Sutton, Nonuniform fast Fourier transforms using
%     min-max interpolation, IEEE-SP, 2003, Vol. 51, pp. 560-574.
% [2] http://www.eecs.umich.edu/~fessler/irt/irt/nufft/
% [3] Leslie Greengard and June-Yub Lee, Accelerating the nonuniform fast 
%     Fourier transform, SIAM Review, Vol. 46, pp. 443-454.
% [4] http://www.mathworks.com/matlabcentral/fileexchange/25135-nufft-nfft-usfft
%
% (c) by Hannes Nickisch & Alexander Loktyushin, 
%                                  MPI for Biological Cybernetics, 2011 March 03

function  A = matFFTnu(sz,k,type,complex)

  if numel(sz)==3 && sz(3)==1, sz = sz(1:2); k = k(:,1:2); end % pseudo 3d image
  A.imsz = sz;
  A.ndims = numel(sz);                            % Are we doing 2d or 3d stuff?
  if A.ndims<2 || A.ndims>3, error('We support only 2d and 3d nufft.'), end
  sz = [size(k,1),prod(A.imsz)];  
  if nargin<2, error('We need at least two input arguments.'), end
  if nargin<3, A.type = 'ferrara'; else A.type = type; end         % set default

  if any(strfind(A.type,'full')) || any(strfind(A.type,'dense'))   % full matrix
    c = ceil((A.imsz-1)/2)+1;                                           % center
    gvec = cell(A.ndims,1); for i=1:A.ndims, gvec{i} = (1:A.imsz(i))-c(i); end
    [gvec{:}] = ndgrid(gvec{:}); g = reshape(cat(A.ndims+1,gvec{:}),[],A.ndims);
    A = exp(-2i*pi*conj(k)*g')/sqrt(sz(2)); return % g is integer Cartesian grid
  end

  A.ogf = 2;                                               % overgridding factor
  A.kbwsz=[]; A.fftsz=[]; A.accuracy=[]; A.k=[]; A.w=[];
  A.f=[]; A.P=[]; A.id=[]; A.in=[]; A.dP=cell(A.ndims,1); A.di=0; A.F=[]; % init
  if     any(strfind(A.type,'fessler'))                        % Fessler's nufft
    kbwsz = 3; A.k = k;                              % Kaiser-Bessel window size
    A.kbwsz = kbwsz*ones(1, A.ndims);
    A.fftsz = fix(A.ogf*A.imsz);
    [A.f, A.P] = nufft_init(A.k, A.imsz, A.fftsz, A.kbwsz, 0);      % precompute
  elseif any(strfind(A.type,'ferrara'))                        % Ferrara's nufft
    A.accuracy = 2;                                % 1..6, high = slow & precise
    if ~all(mod(A.imsz,2)), error('Method ferrara requires odd sizes.'), end
    % remember outside grid points
    kk = bsxfun(@times,k,A.imsz);
    kk = bsxfun(@plus,kk,ceil((A.imsz-1)/2)+1);
    A.in =    0<kk(:,1) & kk(:,1)<=A.imsz(1) & 0<kk(:,2) & kk(:,2)<=A.imsz(2);
    if A.ndims>2, A.in = A.in & 0<kk(:,3) & kk(:,3)<=A.imsz(3); end% inside grid
    k = k(A.in,:); clear  kk
    if A.ndims==3                                       % shrinkage compensation
      k = k.*repmat(((A.imsz)./(A.imsz-1)), size(k,1),1);
    end
    A.id = abs(k(:,1))<=1/2;  % remove points outside image area after expansion
    for dim=2:A.ndims, A.id = A.id & abs(k(:,dim))<=1/2; end
    A.k = -k(A.id,:);
  elseif any([strfind(lower(A.type),'near'), strfind(lower(A.type),'lin'), ...
                                             strfind(lower(A.type),'cub')])
     A.k = 2*bsxfun(@times, k', A.imsz(:) );         % get back scale, padding 2
     psz = fix(A.imsz/2);                                         % padding size
     A.k = bsxfun(@plus, A.k, A.imsz(:)+1 );                   % add back center
     A.F = matFFTN(2*A.imsz)*matZpad(A.imsz,psz,A.imsz-psz);    % padded Fourier
     A.f = 2;                                    % scaling factor due to padding
     
     if any(strfind(lower(A.type),'mex')), mode = 'mex'; else mode = 'mat'; end
     A.P = matResample(2*A.imsz,A.k,A.type,mode);
  else
    error('Unknown nufft type.')
  end

  if nargin<4, complex = 1; end
  if numel(complex)==1, complex = complex*[1,1]; end
  A = class(A,mfilename,mat(sz,complex));