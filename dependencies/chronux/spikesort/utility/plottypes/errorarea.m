function [linehandles, patchhandles] = errorarea(x,y,l,u)
%ERRORAREA         Plot with confidence region.
%   ERRORAREA(X,Y,L,U) is similar to ERRORBAR(X,Y,L,U), except that the
%   confidence bounds are drawn as a shaded patch above and below the
%   plot.  ERRORAREA(X,Y,E) and ERRORAREA(Y,E) are similarly allowed and
%   are analagous to their ERRORBAR equivalents.
%
%   [LINEHANDLE,PATCHHANDLE] = ERRORAREA(...) returns handles to the line
%   and patch objects.
%
%   See also ERRORBAR.

% Parse the data arguments
if (min(size(x)) == 1),  x = x(:);  end;
[npts,nlns] = size(x);
switch (nargin),
	case 4,   % do nothing
	case 3,   u = l;
	case 2,   l = y;  u = y;  y = x;  x = (1:npts)' * ones(1,nlns);
	otherwise,  error('Invalid syntax.');
end
if (nlns == 1),  y = y(:); l = l(:); u = u(:);  end;
if (~isequal(size(x), size(y), size(u), size(l)))
	error('The sizes of X, Y, L and U must be the same.');
end


% Plot the main line
linehandles = plot(x,y);

% Plot the error patches ...
patchhandles = zeros(1,nlns);
for ln = 1:nlns
	xwrap = [x(:,ln)',fliplr(x(:,ln)')];
	ywrap = [(y(:,ln)+u(:,ln))', flipud(y(:,ln)-l(:,ln))'];
	patchhandles(ln) = patch(xwrap, ywrap, get(linehandles(ln), 'Color'), ...
		                     'EdgeColor', 'none', 'FaceAlpha', 0.25);
end

if (nargout == 0)
	clear linehandles patchhandles
end