% matrix object for 2d convolutions in the frequency domain.
%
% Inputs:
%   f        the left-hand-side of the convolution
%   sx       the size of the right-hand-side of the convolution
%   shape    one of 'circ', 'same' or 'valid'                   [default 'circ']
%
% Example of usage:
%   F = matConv2(f, sx, 'same');
%   y = F*x;    % same as y = imfilter(x, f, 'same');
%
% Michael Hirsch * 20/07/2010
% wrapped by Hannes Nickisch, MPI for Biological Cybernetics, 2010 August 23

function A = matConv2(f, sx, shape, complex)

  sf = size(f);        % size of the psf, not of the matrix

  if any(sf < sx) && any(sf > sx)
    error('[%s.m] size missmatch', mfilename);
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

  F = fftn(f, sF);                                    % precompute fft of filter
  sz = [prod(sy), prod(sx)];                                % calculate its size

  % store all stuff in the structure
  A.shape = shape;    % shape of the convolution
  A.f     = f;        % PSF
  A.F     = F;        % precomputed FFT
  A.sx    = sx;       % input dimension
  A.sy    = sy;       % output dimension
  A.sf    = sf;       % size of the PSF
  A.sF    = sF;       % size of the FFT
  A.sw    = sw;       % size of the PSF

  if nargin<4, complex = 1; end
  if numel(complex)==1, complex = complex*[1,1]; end
  A = class(A,mfilename,mat(sz,complex));