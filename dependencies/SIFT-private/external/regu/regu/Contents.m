% Regularization Tools.
% Version 4.1  9-march-08.
% Copyright (c) 1993 and 1998 by Per Christian Hansen and IMM.
%
% Demonstration.
%   regudemo  - Tutorial introduction to Regularization Tools.
%
% Test problems.
%   baart     - Fredholm integral equation of the first kind.
%   blur      - Image deblurring test problem with structured matrix.
%   deriv2    - Computation of the second derivative.
%   foxgood   - Severely ill-posed problem.
%   gravity   - One-dimensional gravity surveying problem.
%   heat      - Inverse heat equation.
%   i_laplace - Inverse Laplace transformation.
%   parallax  - Stellar parallax problem with 28 fixed observations.
%   phillips  - Philips' "famous" test problem.
%   shaw      - One-dimensional image restoration problem.
%   spikes    - Test problem with a "spiky" solution.
%   tomo      - Two-dimensional tomography problem with sparse matrix.
%   ursell    - Integral equation with no square integrable solution.
%   wing      - Test problem with a discontinuous solution.
%
% SVD- and GSVD-based regularization routines.
%   discrep   - Minimizes the solution (semi-)norm subject to an upper
%               bound on the residual norm (discrepancy principle).
%   dsvd      - Computes the damped SVD/GSVD solution.
%   lsqi      - Minimizes the residual norm subject to an upper bound
%               on the (semi-)norm of the solution.
%   mtsvd     - Computes the modified TSVD solution.
%   tgsvd     - Computes the truncated GSVD solution.
%   tikhonov  - Computes the Tikhonov regularized solution.
%   tsvd      - Computes the truncated SVD solution.
%   ttls      - Computes the truncated TLS solution.
%
% Iterative regularization routines.
%   art       - Algebraic reconstruction technique (Kaczmarz's method).
%   cgls      - Computes the least squares solution based on k steps
%               of the conjugate gradient algorithm.
%   lsqr_b    - Computes the least squares solution based on k steps
%               of the LSQR algorithm (Lanczos bidiagonalization).
%   maxent    - Computes the maximum entropy regularized solution.
%   mr2       - Variant of MINRES with starting vector Ab.
%   nu        - Computes the solution based on k steps of Brakhage's
%               iterative nu-method.
%   pcgls     - Same as cgls, but for general-form regularization.
%   plsqr_b   - Same as lsqr, but for general-form regularization.
%   pmr2      - Same as mr2, but for general-form regularization.
%   pnu       - Same as nu, but for general-form regularization.
%   prrgmres  - Same as rrgmres, but for general-form regularization.
%   rrgmres   - Variant of GMRES with starting vector Ab.
%   splsqr    - Computes an approximate Tikhonov solution via the
%               subspace preconditioned LSQR algorithm.
%
% Analysis routines.
%   corner    - Locates the corner of a discrete L-curve.
%   fil_fac   - Computes filter factors for some regularization methods.
%   gcv       - Plots the GCV function and computes its minimum.
%   l_corner  - Locates the L-shaped corner of the L-curve.
%   l_curve   - Computes the L-curve, plots it, and computes its corner.
%   lagrange  - Plots the Lagrange function ||Ax-b||^2 + lambda^2*||Lx||^2,
%               and its derivative.
%   ncp       - Plots normalized cumulative periodograms (NCPs) and finds
%               the one closest to a straight line.
%   picard    - Plots the (generalized) singular values, the Fourier
%               coefficient for the right-hand side, and a (smoothed curve
%               of) the solution's Fourier-coefficients. 
%   plot_lc   - Plots an L-curve.
%   quasiopt  - Plots the quasi-optimality function and computes its minimum.
%
% Routines for transforming a problem in general form into one in
% standard form, and back again.
%   gen_form  - Transforms a standard-form solution back into the
%               general-form setting.
%   std_form  - Transforms a general-form problem into one in
%               standard form.
% 
% Utility routines.
%   bidiag    - Bidiagonalization of a matrix by Householder transformations.
%   cgsvd     - Computes the compact generalized SVD of a matrix pair.
%   csvd      - Computes the compact SVD of an m-by-n matrix.
%   get_l     - Produces a p-by-n matrix which is the discrete
%               approximation to the d'th order derivative operator.
%   lanc_b    - Performs k steps of the Lanczos bidiagonalization
%               process with/without reorthogonalization.
%   regutm    - Generates random test matrices for regularization methods.
%
% Auxiliary routines required by some of the above routines.
%   app_hh    - Applies a Householder transformation from the left.
%   gen_hh    - Generates a Householder transformation.
%   lsolve    - Inversion with A-weighted generalized inverse of L.
%   ltsolve   - Inversion with transposed A-weighted inverse of L.
%   pinit     - Initialization for treating general-form problems.
%   spleval   - Computes points on a spline or spline curve.
 
% The following four routines are not documented, since they are only used
% internally by gcv, l_corner, and quasiopt, respectively.  They cannot be
% located as private functions.
%   gcvfun    - Computes the GCV function
%   lcfun     - Computes the curvature of the L-curve
%   ncpfun    - Computes the NCP's distance to a straight line.
%   quasifun  - Computes the quasi-optimality function.
%
% For backward compatibility, the function l_corner uses the Spline
% Toolbox when available, otherwise is used the new function corner.