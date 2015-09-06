function points = spleval(f)
%SPLEVAL Evaluation of a spline or spline curve.
%
% points = spleval(f)
%
% Computes points on the given spline or spline curve f between
% its extreme breaks.

% Original routine fnplt by C. de Boor / latest change: Oct. 25, 1997
% Simplified by Per Christian Hansen, IMM, 04/16/98.

% Set default number of points.
npoints = 300;

if (f.form(1)=='B'), f = sp2pp(f); end

[breaks,coefs,l,k,d] = ppbrk(f);
x = breaks(1) + (0:npoints)*((breaks(l+1)-breaks(1))/npoints);
v=ppual(f,x);

if (d==1), points=[x;v]; else points = v; end