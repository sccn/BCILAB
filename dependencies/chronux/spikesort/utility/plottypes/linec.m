function h = linec(x,y,z,c)
%LINEC             Plots a line with color varying along its length.
%   LINEC(X,Y,C) draws the parametric plot Y(k) versus vector X(k) for
%   k = 1:length(X) with colors at each vertex given by the vector C
%   (scaled into the figure colormap).  The color between vertices is
%   interpolated.
%
%   LINEC(Y,C) is equivalent to LINEC(1:length(Y), Y, C).
%
%   LINEC(X,Y,Z,C) plots [X(k), Y(k), Z(k)] with colors given by C(k).
%
%   A patch object is used to actually render the line.  Note that a flaw
%   in the Matlab 'painters' renderer causes extremely slow execution on
%   some machines when zooming in on these patch 'lines'; for this reason,
%   LINEC sets the renderer of the current figure to 'OpenGL'. 
%
%   H = LINEC(...) returns a handle to the line.  The 'line' is actually
%   a patch object, so the standard parameter/value pairs used with
%   Matlab's builtin PLOT command are not supported.
%
%   Matrix arguments for X,Y are currently not supported.

%%%%%%%%%%%%%%%%% ARGUMENT CHECKING %%%%%%%%%%%%%%%%%
if (nargin > 4),        error('Only (x,y,z,c) arguments supported.');
elseif (nargin == 4),
elseif (nargin == 3),   c = z;  z = repmat(0,[1,0]);
elseif (nargin == 2),   c = y;  y = x;  x = (1:length(y))';  z = repmat(0,[1,0]);
else                    error('At least two arguments are required.');
end

xlen = length(x);  ylen = length(y);  zlen = length(z);  clen = length(c);

if (isempty(z))
	if (~isequal(xlen,ylen,clen)), error('Vector dimensions must agree.'); end;
else
	if (~isequal(xlen,ylen,zlen,clen)), error('Vector dimensions must agree.');  end;
end
	
if ((ndims(x) > 2) || (~any(size(x) == 1)) || (ndims(y) > 2) || (~any(size(y) == 1)) || ...
	(ndims(z) > 2) || (~any(size(z) == 1)) || (ndims(c) > 2) || (~any(size(c) == 1)))
    error('Arguments must be row or column vectors.');
end

x = x(:); y = y(:);  z = z(:);  c = c(:);    % Force column vectors

%%%%%%%%%%%%%%%%%%% DO THE WORK %%%%%%%%%%%%%%%%%%%%%
xdata = [x; flipud(x)];
ydata = [y; flipud(y)];
cdata = [c; flipud(c)];
zdata = [z; flipud(z)];

if (isempty(zdata))
	h = patch('XData', xdata, 'YData', ydata, 'CData', cdata);  view(2);
else
	h = patch('XData', xdata, 'YData', ydata, 'ZData', zdata, 'CData', cdata');   view(3);
end
set(h, 'EdgeColor', 'interp', 'LineWidth', 2);
set(gcf, 'Renderer', 'OpenGL');

if (nargout == 0)
    clear h
end
