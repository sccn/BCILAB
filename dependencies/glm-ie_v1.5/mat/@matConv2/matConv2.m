% matrix object for 2d convolutions in the frequency domain.
%
% Inputs:
%   f        the left-hand-side of the convolution
%   sx       the size of the right-hand-side of the convolution
%   shape    one of 'circ', 'same' or 'valid'                   [default 'circ']
%   type     one of 'conv' or 'corr'                            [default 'conv']
%
% Example of usage:
%   F = matConv2(f, sx, 'same');
%   y = F*x;    % same as y = imfilter(x, f, 'same', 'conv');
%
% Michael Hirsch * 20/07/2010
% wrapped by Hannes Nickisch, MPI for Biological Cybernetics, 2011 October 20

function A = matConv2(f, sx, shape, type, complex)

  sf = size(f);        % size of the psf, not of the matrix

  if any(sf < sx) && any(sf > sx)
    error('[%s.m] size missmatch', mfilename)
  end

  sF = max(sx, sf);    % the bigger one
  sw = min(sx, sf);    % the smaller one

  if nargin<3, shape = 'circ'; end
  switch shape
    case 'circ'
      sy = sF;
    case 'valid'
      sy = sF - sw + 1;
    case 'same'
      sy = sF;
      sF = sx + sf - 1;
    otherwise
      error('[matConv2.m] Unknown shape. Choose between ''circ'', ''same'' and ''valid''');
  end

  if nargin<4, type = 'conv'; end
  switch type
    case 'conv'
    case 'corr'
      f = rot90(f,2);
    otherwise
      error('[matConv2.m] Unknown type. Choose between ''conv'' and ''corr''');
  end

  F = fftn(f, sF);                                    % precompute fft of filter
  sz = [prod(sy), prod(sx)];                                % calculate its size

  % store all stuff in the structure
  A.shape   = shape;    % shape of the convolution
  A.type    = type;     % conv or corr
  A.f       = f;        % PSF
  A.F       = F;        % precomputed FFT
  A.sx      = sx;       % input dimension
  A.sy      = sy;       % output dimension
  A.sf      = sf;       % size of the PSF
  A.sF      = sF;       % size of the FFT
  A.sw      = sw;       % size of the PSF

  if nargin<5, complex = 1; end
  A.complex = complex;            % just memorise the complex type for later use
  if numel(complex)==1, complex = complex*[1,1]; end
  A = class(A,mfilename,mat(sz,complex));