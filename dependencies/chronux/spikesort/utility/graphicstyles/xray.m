function map = xray(m)
%XRAY              Inverted gray-scale color map.
%   XRAY(M) returns an M-by-3 matrix containing an inverted gray-scale
%   colormap (from white to black, rather than black-to-white).
%   XRAY, by itself, is the same length as the current colormap.

if (nargin < 1), m = size(colormap,1);  end;
map = flipud(gray(m));