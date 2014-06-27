function K = compute_cubic_spline_kernel(x,y);
nx = size(x,1);
ny = size(y,1);
Kmax = max(repmat(abs(x),1,ny),repmat(abs(y)',nx,1) );
Kmin = min(repmat(abs(x),1,ny),repmat(abs(y)',nx,1) );
K = Kmin .* Kmin .* ( 3 * Kmax - Kmin );
K = K .* (  x*y' > 0 );
 
