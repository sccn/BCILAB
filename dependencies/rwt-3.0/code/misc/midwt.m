function [y,L] = midwt(x,h,L)
%    [x,L] = midwt(y,h,L);
% 
%    Function computes the inverse discrete wavelet transform x for a 1D or
%    2D input signal y using the scaling filter h.
%
%    Input:
%	y : finite length 1D or 2D input signal (implicitly periodized)
%           (see function mdwt to find the structure of y)
%       h : scaling filter
%       L : number of levels. In the case of a 1D signal, length(x) must be
%           divisible by 2^L; in the case of a 2D signal, the row and the
%           column dimension must be divisible by 2^L.  If no argument is
%           specified, a full inverse DWT is returned for maximal possible
%           L.
%
%    Output:
%       x : periodic reconstructed signal
%       L : number of decomposition levels
%
%    1D Example:
%       xin = makesig('LinChirp',8);
%       h = daubcqf(4,'min');
%       L = 1;
%       [y,L] = mdwt(xin,h,L);
%       [x,L] = midwt(y,h,L)
%
%    1D Example's  output:
%
%       x = 0.0491 0.1951 0.4276 0.7071 0.9415 0.9808 0.6716 0.0000
%       L = 1
%
%    See also: mdwt, mrdwt, mirdwt
%
%Author: Markus Lang  <lang@jazz.rice.edu>
