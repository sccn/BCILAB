% Translation only motion
%
% A = matFFTtrans(sz,trans,rot,complex,type)
%
%
% (c) Alexander Loktyushin,     MPI for Biological Cybernetics, 2012

function A = matFFTtrans(sz, trans, complex)

  %% general info
  A.ndims = numel(sz);                            % Are we doing 2d or 3d stuff?
  if A.ndims<2 || A.ndims>3, error('We support only 2d and 3d.'), end
  A.imsz = sz; sz = prod(A.imsz)*[1,1];
  if nargin<3, complex = 1; end, if numel(complex)==0, complex = 1; end
  if numel(complex)==1, complex = complex*[1,1]; end
  center = ceil((A.imsz-1)/2)+1; A.center = center;          % center in k-space
    
  %% derivatives
  A.dkind = 'none';                    % kind of derivative, either rot or trans
  A.di = 0;               % which derivative, 0 means all of the respective kind

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
  
  A.F = matFFTN(A.imsz); 
  
  % compute translational correction for not rotating around the origin
  A.tc = 1;

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

