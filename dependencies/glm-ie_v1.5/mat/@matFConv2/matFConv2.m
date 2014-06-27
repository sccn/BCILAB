% FFT-free 2D Convolution
%
% Inputs:
%   f        the left-hand-side of the convolution
%   sx       the size of the right-hand-side of the convolution
%   shape    one of 'circ', 'same' or 'valid'                   [default 'circ']
%
% Example of usage:
%   F = matFConv2(f, sx, 'same');
%   y = F*x;    % same as y = imfilter(x, f, 'same');
%

function A = matFConv2(f, sx, shape, complex)

  sf = size(f);        % size of the psf, not of the matrix

  if any(sf < sx) && any(sf > sx)
    error('[%s.m] size missmatch', mfilename);
  end

  sF = max(sx, sf);    % the bigger one
  sw = min(sx, sf);    % the smaller one

  sy = sF;
  sF = sx;

  F = fftshift(fftn(f, sF));                                    % precompute fft of filter
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