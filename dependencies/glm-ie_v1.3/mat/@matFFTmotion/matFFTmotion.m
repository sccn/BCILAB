% Two and threedimensional (centered) Fourier transformation
% where the input can be subject to rigid body motion
%
% A = matFFTmotion(sz,trans,rot,complex,type)
%
% sz: size of the array; either 2d or 3d
%
% 2d i.e. numel(sz)=2
%  trans: [2,1] or [2,sz(2)], the translation in pixels,           default [0;0]
%  rot:   [1,1] or [1,sz(2)], the rotation angle in radian,            default 0
%
% 3d i.e. numel(sz)=3
%  trans: [3,1] or [3,sz(2),sz(3)], the translation in pixels,   default [0;0;0]
%  rot:   [3,1] or [3,sz(2),sz(3)], the rotation axis, whose norm equals the 
%                                    rotation angle in radian,   default [0;0;0]
% Note that by convention, translation is done after rotation.
% Note further that the matrix relies on the fact that the pixels/voxels are
% cubes i.e. we have the same resolution along every dimension.
%
% complex: [1,2] complex number treatment A:IN^n->OUT^m,           default [1,1]
%                 first number specifies output field OUT^m and the second the
%                 input field IN^m of the matrix A
%                 1 -  real(z) + 1i*imag(z)    IC:   native complex numbers
%                 2 - [real(z); imag(z)]       IR^2: stacked real numbers
%                 3    real(z)                 IR:   real part only
%
% type: string indicating the nufft interpolation type,
%       we support 'nearest', 'linear', 'cubic',
%                  'fessler' and 'ferrara'                      default 'linear'
%  'nearest', 'linear', 'cubic'
%    We use nearest neighbor, linear or cubic interpolation after oversampling.
%    One can use type{mex,mat} to indicate whether the mex or the mat option is
%    chosen in the matResample class. E.g 'cubicmex' is possible.
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
%                                  MPI for Biological Cybernetics, 2011 March 02

function A = matFFTmotion(sz, trans, rot, complex, type)

  %% general info
  A.ndims = numel(sz);                            % Are we doing 2d or 3d stuff?
  if A.ndims<2 || A.ndims>3, error('We support only 2d and 3d.'), end
  A.imsz = sz; sz = prod(A.imsz)*[1,1];
  if nargin<4, complex = 1; end, if numel(complex)==0, complex = 1; end
  if numel(complex)==1, complex = complex*[1,1]; end
  center = ceil((A.imsz-1)/2)+1; A.center = center;          % center in k-space
  if nargin<5, type = 'linear'; end, A.type = type;   % set (default) nufft type

  %% derivatives
  A.dkind = 'none';                    % kind of derivative, either rot or trans
  A.di = 0;               % which derivative, 0 means all of the respective kind

  %% rotation
  % The outcome of the computations below is
  %  1) a matrix A.RF implementing the Fourier transform followed by a rotation,
  %     which is done by means of a non-uniform Fourier transform with griding
  %     interpolation and oversampling to reduce artefacts;
  %  2) a bunch of vectors A.k indicating the rotated k-space points, some of
  %     which can be placed outside the image rectangle;
  %  3) a set of bunches of vectors A.dk containing the perturbed versions of
  %     A.k needed to compute the rotational derivatives;
  %  4) an array of booleans indicating whether the vectors in A.k are placed
  %     inside the image rectangle.
  nrot = A.ndims*(A.ndims-1)/2;                % number of rotational parameters
  if nargin<3,        rot = zeros(nrot,1); end
  if numel(rot)<nrot, rot = zeros(nrot,1); end, A.rot = rot(:);
  if numel(A.rot)>nrot                                 % bring into proper shape
    A.rot = reshape(A.rot,[nrot,prod(A.imsz(2:end))]);
  end
  if any(strfind(lower(A.type),'fessler'))
    ctrans  = zeros(A.ndims,size(A.rot,2)); % center translation due to rotation
    dctrans = zeros(A.ndims,size(A.rot,2),nrot);
  end
  sc = A.imsz/A.imsz(1);             % correction for A.imsz(1) ~= A.imsz(2:end)
  k = zeros(A.ndims,prod(A.imsz));                             % allocate memory
  A.dk = repmat(k,[1,1,nrot]);                    % init with no rotation at all
  T = size(A.rot,2);                             % number of different rotations
  dR = cell(1,nrot);                      % derivatives of the rotation matrices
  if A.ndims==2                                          % construct target grid
    [i2,i1]=meshgrid(1:A.imsz(2),1:A.imsz(1)); ii=[i1(:)';i2(:)']; clear i1 i2
  else
    [i2,i1,i3] = meshgrid(1:A.imsz(2),1:A.imsz(1),1:A.imsz(3));
    ii = [i1(:)'; i2(:)'; i3(:)']; clear i1 i2 i3
  end
  for t = 1:T                                              % loop over scanlines
    if T==1
      i = 1:prod(A.imsz);                                  % all indices at once
    else
      i = (t-1)*A.imsz(1) + (1:A.imsz(1));      % indices of individual scanline
    end
    [R,dR{:}] = rotmatsc(A.rot(:,t),sc);   % rotation matrix and its derivatives
    k(:,i) = bsxfun(@minus,ii(:,i),center(:));                   % centered grid
    for j=1:nrot, A.dk(:,i,j) = dR{j}'*k(:,i); end       % dk/dx grid derivative
    k(:,i) = R'*k(:,i);                                   % non-integer src grid
    if any(strfind(lower(A.type),'fessler'))                        % correction
      ctrans(:,t) = -diag(sc)*R*diag(1./sc)*A.imsz(:)/2;   % off-center-rotation
      for j=1:nrot
        dctrans(:,t,j) = -diag(sc)*dR{j}*diag(1./sc)*A.imsz(:)/2;
      end
    end
  end
  A.k = bsxfun(@times,k,1./A.imsz(:)); A.dk = bsxfun(@times,A.dk,1./A.imsz(:));
  A.RF = matFFTnu(A.imsz,A.k',A.type); % Fourier followed by rot, using gridding

  %% translation
  % The outcome of the computations below is
  %  1) a vector A.t implementing the translation;
  %  2) a set of vectors A.dt needed for the translation derivatives;
  %  3) a vector A.tc needed to correct for the rotation that is not exactly
  %     done around the center.
  if nargin<2,             trans = zeros(A.ndims,1); end
  if numel(trans)<A.ndims, trans = zeros(A.ndims,1); end, A.trans = trans(:);
  if numel(A.trans)>A.ndims                            % bring into proper shape
    A.trans = reshape(A.trans,[A.ndims,A.imsz(2:end)]);
  end
  [A.t,A.dt] = trans_diag(A.imsz,A.trans);  % diag translation matrix and derivs
  
  % compute translational correction for not rotating around the origin
  if any(strfind(lower(A.type),'fessler'))
    if numel(ctrans)>A.ndims                           % bring into proper shape
      ctrans = reshape(ctrans,[A.ndims,A.imsz(2:end)]);
    end
    [A.tc,dtc] = trans_diag(A.imsz,ctrans);      % diagonal trans matrix, derivs
    A.dtc = zeros(numel(A.tc),nrot); dtc = reshape(dtc,[A.imsz,A.ndims]);
    for j=1:nrot
      dctrans_j = dctrans(:,:,j).';
      if numel(dctrans_j)>A.ndims
        dctrans_j = reshape(dctrans_j,[1,A.imsz(2:end),A.ndims]);
      else
        dctrans_j = reshape(dctrans_j,[ones(1,A.ndims),A.ndims]);
      end
      A.dtc(:,j) = reshape(sum(bsxfun(@times,dtc,dctrans_j),A.ndims+1),[],1);
    end
  else
    A.tc = 1; A.dtc = zeros(1,nrot);
  end

  %% construct matrix class, inherit from mat
  A = class(A,mfilename,mat(sz,complex));


% Compute diagonal translation equivalent in Fourier space
function [t,dt] = trans_diag(sz,tr)

  ndims = numel(sz); N = prod(sz);      % number of dimensions, number of pixels
  deriv = nargout>1;               % flag saying, whether we compute derivatives
  t = ones(sz); if deriv, dt = ones(N,ndims); end              % allocate memory
  for d=1:ndims
    rp = ((1:sz(d))-fix(1+sz(d)/2))/sz(d);    % ramp in direction of dimension d
    rp = -2i*pi*reshape(rp,[ones(1,d-1),sz(d),ones(1,ndims-d)]);
    erp = exp(bsxfun(@times,rp,tr(d,:,:)));           % tr is singleton in dim 1
    t = bsxfun(@times, t, erp);
    if deriv
      for dd=1:ndims
        dtdd = bsxfun(@times, reshape(dt(:,dd),sz), erp);
        if dd==d
          dtdd = bsxfun(@times, dtdd, rp);
        end
        dt(:,dd) = reshape(dtdd,N,1);
      end
    end
  end


% Compute scaled rotation matrices and their derivatives
function [R,varargout] = rotmatsc(rot,sc)

  R = diag(1./sc)*rotmat(rot)*diag(sc);                 % scaled rotation matrix
  if nargout>1
    varargout = cell(nargout-1, 1);  % allocate right number of output arguments
    [varargout{:}] = drotmat(rot);
    for i=1:nargout-1
      varargout{i} = diag(1./sc)*varargout{i}*diag(sc);      % scaled derivative
    end
  end
