function [y,L] = mdwt(x,h,L)
%    [y,L] = mdwt(x,h,L);
%
%    Function computes the discrete wavelet transform y for a 1D or 2D input
%    signal x using the scaling filter h.
%
%    Input:
%	x : finite length 1D or 2D signal (implicitly periodized)
%       h : scaling filter
%       L : number of levels. In the case of a 1D signal, length(x) must be
%           divisible by 2^L; in the case of a 2D signal, the row and the
%           column dimension must be divisible by 2^L. If no argument is
%           specified, a full DWT is returned for maximal possible L.
%
%    Output:
%       y : the wavelet transform of the signal 
%           (see example to understand the coefficients)
%       L : number of decomposition levels
%
%    1D Example:
%       x = makesig('LinChirp',8);
%       h = daubcqf(4,'min');
%       L = 2;
%       [y,L] = mdwt(x,h,L)
%
%    1D Example's  output and explanation:
%
%       y = [1.1097 0.8767 0.8204 -0.5201 -0.0339 0.1001 0.2201 -0.1401]
%       L = 2
%
%    The coefficients in output y are arranged as follows
%
%       y(1) and y(2) : Scaling coefficients (lowest frequency)
%       y(3) and y(4) : Band pass wavelet coefficients
%       y(5) to y(8)  : Finest scale wavelet coefficients (highest frequency)
%
%    2D Example:
%
%       load test_image        
%       h = daubcqf(4,'min');
%       L = 1;
%       [y,L] = mdwt(test_image,h,L);
%
%    2D Example's  output and explanation:
%
%       The coefficients in y are arranged as follows.
%
%              .------------------.
%              |         |        |
%              |    4    |   2    |
%              |         |        |
%              |   L,L   |   H,L  |
%              |         |        |
%              --------------------
%              |         |        |
%              |    3    |   1    |
%              |         |        |
%              |   L,H   |  H,H   |
%              |         |        |
%              `------------------'
%       
%       where 
%            1 : High pass vertically and high pass horizontally
%            2 : Low pass vertically and high pass horizontally
%            3 : High pass vertically and low  pass horizontally
%            4 : Low pass vertically and Low pass horizontally 
%                (scaling coefficients)
%
%
%
%
%    See also: midwt, mrdwt, mirdwt
%
%Author: Markus Lang  <lang@jazz.rice.edu>
