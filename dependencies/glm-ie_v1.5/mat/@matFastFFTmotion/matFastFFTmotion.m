% Two and threedimensional (centered) Fourier transformation
% where the input can be subject to rigid body motion
%
% A = matFastFFTmotion(sz,trans,rot,complex,type)
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
%    We use cubic interpolation after oversampling. One can set type{mex,gpu}
%    to indicate whether CPU or GPU implementation of resampling is used.
%
% (c) Alexander Loktyushin,     MPI for Biological Cybernetics, 2011

function A = matFastFFTmotion(sz, trans, rot, complex, type)

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
%{
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
  end
%}
  setenv('OMP_NUM_THREADS', num2str(feature('numCores')));

  dkR = cell(1,nrot);                      % derivatives of the rotation matrices
  [k, dkR{:}] = grid_create(ii, A.rot, T, A.imsz, center, sc, nrot);

  for i = 1:nrot
    A.dk(:,:,i) = dkR{i};
  end;

  A.k = bsxfun(@times,k,1./A.imsz(:)); A.dk = bsxfun(@times,A.dk,1./A.imsz(:));
    
  % matFFTnu supplement
  k = A.k'; nusz = A.imsz;
  if numel(nusz)==3 && nusz(3)==1, nusz = nusz(1:2); k = k(:,1:2); end % pseudo 3d image
  A.nuimsz = nusz;
  A.nundims = numel(nusz);                            % Are we doing 2d or 3d stuff?
  
  A.ogf = 2;                                               % overgridding factor
  A.kbwsz=[]; A.fftsz=[]; A.accuracy=[]; A.k=[]; A.w=[];
  A.f=[]; A.P=[]; A.id=[]; A.in=[]; A.dP=cell(A.nundims,1); A.di=0; A.F=[]; % init
  
  A.k = 2*bsxfun(@times, k', A.nuimsz(:) );         % get back scale, padding 2
  psz = fix(A.nuimsz/2);                                         % padding size
  A.k = bsxfun(@plus, A.k, A.nuimsz(:)+1 );                   % add back center
  A.F = matFFTN(2*A.nuimsz)*matZpad(A.nuimsz,psz,A.nuimsz-psz);    % padded Fourier
  A.f = 2;                                    % scaling factor due to padding
     
  if any(strfind(lower(A.type),'mex'))
      mode = 'mex';
  elseif any(strfind(lower(A.type),'gpu'))
      mode = 'mex';                           % !! disable GPU capability for this class for now, otherwise set mode = 'gpu';
  else mode = 'mat'; 
  end
  
  A.mode = mode;
   

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
  A.tc = 1; A.dtc = zeros(1,nrot);

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

