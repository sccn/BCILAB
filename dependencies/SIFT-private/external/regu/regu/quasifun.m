function Q = quasifun(lambda,s,xi,dsvd)

% Auxiliary routine for quasiopt.  PCH, IMM, 12/29/97.

if (nargin==3)
   f = (s.^2)./(s.^2 + lambda^2);
else
   f = s./(s + lambda);
end

Q = norm((1 - f).*f.*xi);