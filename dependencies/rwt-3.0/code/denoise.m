function [xd,xn,option] = denoise(x,h,type,option)
%    [xd,xn,option] = denoise(x,h,type,option); 
%
%    DENOISE is a generic program for wavelet based denoising.
%    The program will denoise the signal x using the 2-band wavelet
%    system described by the filter h using either the traditional 
%    discrete wavelet transform (DWT) or the linear shift invariant 
%    discrete wavelet transform (also known as the undecimated DWT
%    (UDWT)). 
%
%    Input:  
%       x         : 1D or 2D signal to be denoised
%       h         : Scaling filter to be applied
%       type      : Type of transform (Default: type = 0)
%                   0 --> Discrete wavelet transform (DWT)
%                   1 --> Undecimated DWT (UDWT)
%       option    : Default settings is marked with '*':
%                   *type = 0 --> option = [0 3.0 0 0 0 0]
%                   type = 1 --> option = [0 3.6 0 1 0 0]
%       option(1) : Whether to threshold low-pass part
%                   0 --> Don't threshold low pass component 
%                   1 --> Threshold low pass component
%       option(2) : Threshold multiplier, c. The threshold is
%                   computed as: 
%                     thld = c*MAD(noise_estimate)). 
%                   The default values are:
%                     c = 3.0 for the DWT based denoising
%                     c = 3.6 for the UDWT based denoising
%       option(3) : Type of variance estimator
%                   0 --> MAD (mean absolute deviation)
%                   1 --> STD (classical numerical std estimate)
%       option(4) : Type of thresholding
%                   2 --> Soft thresholding
%                   1 --> Hard thresholding
%       option(5) : Number of levels, L, in wavelet decomposition. By
%                   setting this to the default value '0' a maximal
%                   decomposition is used.
%       option(6) : Actual threshold to use (setting this to
%                   anything but 0 will mean that option(3)
%                   is ignored)
%
%    Output: 
%       xd     : Estimate of noise free signal 
%       xn     : The estimated noise signal (x-xd)
%       option : A vector of actual parameters used by the
%                program. The vector is configured the same way as
%                the input option vector with one added element
%                option(7) = type.
%
%  HERE'S AN EASY WAY TO RUN THE EXAMPLES:
%  Cut-and-paste the example you want to run to a new file 
%  called ex.m, for example. Delete out the % at the beginning 
%  of each line in ex.m (Can use search-and-replace in your editor
%  to replace it with a space). Type 'ex' in matlab and hit return.
%
%    Example 1: 
%       h = daubcqf(6); [s,N] = makesig('Doppler'); n = randn(1,N);
%       x = s + n/10;     % (approximately 10dB SNR)
%       figure;plot(x);hold on;plot(s,'r');
%
%       %Denoise x with the default method based on the DWT
%       [xd,xn,opt1] = denoise(x,h);
%       figure;plot(xd);hold on;plot(s,'r');
%
%       %Denoise x using the undecimated (LSI) wavelet transform
%       [yd,yn,opt2] = denoise(x,h,1);
%       figure;plot(yd);hold on;plot(s,'r');
%
% Example 2: (on an image)  
%      h = daubcqf(6);  load lena; 
%      noisyLena = lena + 25 * randn(size(lena));
%      figure; colormap(gray); imagesc(lena); title('Original Image');
%       figure; colormap(gray); imagesc(noisyLena); title('Noisy Image'); 
%       Denoise lena with the default method based on the DWT
%      [denoisedLena,xn,opt1] = denoise(noisyLena,h);
%      figure; colormap(gray); imagesc(denoisedLena); title('denoised Image');
%       
%
%    See also: mdwt, midwt, mrdwt, mirdwt, SoftTh, HardTh, setopt
%
%Author: Jan Erik Odegard  <odegard@ece.rice.edu>

if(nargin < 2)
  error('You need to provide at least 2 inputs: x and h');
end;
if(nargin < 3),
  type = 0;
  option = [];
elseif(nargin < 4)
  option = [];
end;
if(isempty(type)),
  type = 0;
end;
if(type == 0),
  default_opt = [0 3.0 0 2 0 0];
elseif(type == 1),
  default_opt = [0 3.6 0 1 0 0];
else
  error(['Unknown denoising method',10,...
	  'If it is any good we need to have a serious talk :-)']);
end;
option = setopt(option,default_opt);
[mx,nx] = size(x);
dim = min(mx,nx);
if(dim == 1),
  n = max(mx,nx);
else
  n = dim;
end;
if(option(5) == 0),
  L = floor(log2(n));
else
  L = option(5);
end;
if(type == 0), 			% Denoising by DWT
  xd = mdwt(x,h,L);
  if (option(6) == 0),
    tmp = xd(floor(mx/2)+1:mx,floor(nx/2)+1:nx);
    if(option(3) == 0),
      thld = option(2)*median(abs(tmp(:)))/.67;
    elseif(option(3) == 1),
      thld = option(2)*std(tmp(:));
    else
      error('Unknown threshold estimator, Use either MAD or STD');
    end;
  else
    thld = option(6);
  end;
  if(dim == 1)
    ix = 1:n/(2^L);
    ykeep = xd(ix);
  else
    ix = 1:mx/(2^L);
    jx = 1:nx/(2^L);
    ykeep = xd(ix,jx);
  end;
  if(option(4) == 2),
    xd = SoftTh(xd,thld);
  elseif(option(4) == 1),
    xd = HardTh(xd,thld);
  else
    error('Unknown threshold rule. Use either Soft (2) or Hard (1)');
  end;
  if (option(1) == 0),
    if(dim == 1),
      xd(ix) = ykeep;
    else
      xd(ix,jx) = ykeep;
    end;
  end;
  xd = midwt(xd,h,L);
elseif(type == 1), 			% Denoising by UDWT
  [xl,xh] = mrdwt(x,h,L);
  if(dim == 1),
    c_offset = 1;
  else
    c_offset = 2*nx + 1;
  end;
  if (option(6) == 0),
    tmp = xh(:,c_offset:c_offset+nx-1);
    if(option(3) == 0),
      thld = option(2)*median(abs(tmp(:)))/.67;
    elseif(option(3) == 1),
      thld = option(2)*std(tmp(:));
    else
      error('Unknown threshold estimator, Use either MAD or STD');
    end;
  else
    thld = option(6);
  end;
  if(option(4) == 2),
    xh = SoftTh(xh,thld);
    if(option(1) == 1),
      xl = SoftTh(xl,thld);
    end;
  elseif(option(4) == 1),
    xh = HardTh(xh,thld);
    if(option(1) == 1),
      xl = HardTh(xl,thld);
    end;
  else
    error('Unknown threshold rule. Use either Soft (2) or Hard (1)');
  end;
  xd = mirdwt(xl,xh,h,L);
else 					% Denoising by unknown method
  error(['Unknown denoising method',10,...
         'If it is any good we need to have a serious talk :-)']);
end;
option(6) = thld;
option(7) = type;
xn = x - xd; 
