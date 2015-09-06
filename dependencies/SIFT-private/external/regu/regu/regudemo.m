%REGUDEMO Tutorial script for Regularization Tools.

% Per Christian Hansen, IMM, Feb. 21, 2001.

echo on
clf

% Part 1.  The discrete Picard condition
% --------------------------------------
%
% First generate a "pure" test problem where only rounding
% errors are present.  Then generate another "noisy" test
% problem by adding white noise to the right-hand side.
%
% Next compute the SVD of the coefficient matrix A.
%
% Finally, check the Picard condition for both test problems
% graphically.  Notice that for both problems the condition is
% indeed satisfied for the coefficients corresponding to the
% larger singular values, while the noise eventually starts to
% dominate.

[A,b_bar,x] = shaw(32);
randn('seed',41997);
e = 1e-3*randn(size(b_bar)); b = b_bar + e;
[U,s,V] = csvd(A);
subplot(2,1,1); picard(U,s,b_bar);
subplot(2,1,2); picard(U,s,b);
pause, clf

% Part 2.  Filter factors
% -----------------------
%
% Compute regularized solutions to the "noisy" problem from Part 1 
% by means of Tikhonov's method and LSQR without reorthogonalization.
% Also, compute the corresponding filter factors.
%
% A surface (or mesh) plot of the solutions clearly shows their dependence
% on the regularization parameter (lambda or the iteration number).

lambda = [1,3e-1,1e-1,3e-2,1e-2,3e-3,1e-3,3e-4,1e-4,3e-5];
X_tikh = tikhonov(U,s,V,b,lambda);
F_tikh = fil_fac(s,lambda);
iter = 30; reorth = 0;
[X_lsqr,rho,eta,F_lsqr] = lsqr_b(A,b,iter,reorth,s);
subplot(2,2,1); surf(X_tikh), title('Tikhonov solutions'), axis('ij')
subplot(2,2,2); surf(log10(F_tikh)), axis('ij')
                title('Tikh filter factors, log scale')
subplot(2,2,3); surf(X_lsqr(:,1:17)), title('LSQR solutions'), axis('ij')
subplot(2,2,4); surf(log10(F_lsqr(:,1:17))), axis('ij')
title('LSQR filter factors, log scale')
pause, clf

% Part 3.  The L-curve
% --------------------
%
% Plot the L-curves for Tikhonov regularization and for
% LSQR for the "noisy" test problem from Part 1.
%
% Notice the similarity between the two L-curves and thus,
% in turn, by the two methods.

subplot(1,2,1); l_curve(U,s,b); axis([1e-3,1,1,1e3])
subplot(1,2,2); plot_lc(rho,eta,'o'); axis([1e-3,1,1,1e3])
pause, clf

% Part 4.  Regularization parameters
% ----------------------------------
%
% Use the L-curve criterion and GCV to determine the regularization
% parameters for Tikhonov regularization and truncated SVD.
%
% Then compute the relative errors for the four solutions.

lambda_l = l_curve(U,s,b);   axis([1e-3,1,1,1e3]),      pause, clf
k_l = l_curve(U,s,b,'tsvd'); axis([1e-3,1,1,1e3]),      pause, clf
lambda_gcv = gcv(U,s,b);     axis([1e-6,1,1e-9,1e-1]),  pause, clf
k_gcv = gcv(U,s,b,'tsvd');   axis([0,20,1e-9,1e-1]),    pause, clf

x_tikh_l   = tikhonov(U,s,V,b,lambda_l);
x_tikh_gcv = tikhonov(U,s,V,b,lambda_gcv);
if isnan(k_l)
  x_tsvd_l = zeros(32,1); % Spline Toolbox not available.
else
  x_tsvd_l = tsvd(U,s,V,b,k_l);
end
x_tsvd_gcv = tsvd(U,s,V,b,k_gcv);
disp([norm(x-x_tikh_l),norm(x-x_tikh_gcv),...
 norm(x-x_tsvd_l),norm(x-x_tsvd_gcv)]/norm(x))
pause, clf

% Part 5.  Standard form versus general form
% ------------------------------------------
%
% Generate a new test problem: inverse Laplace transformation
% with white noise in the right-hand side.
%
% For the general-form regularization, choose minimization of
% the first derivative.
%
% First display some left singular vectors of SVD and GSVD; then
% compare truncated SVD solutions with truncated GSVD solutions.
% Notice that TSVD cannot reproduce the asymptotic part of the
% solution in the right part of the figure.

n = 16; [A,b,x] = i_laplace(n,2);
b = b + 1e-4*randn(size(b));
L = get_l(n,1);
[U,s,V] = csvd(A); [UU,sm,XX] = cgsvd(A,L);
I = 1;
for i=[3,6,9,12]
  subplot(2,2,I); plot(1:n,V(:,i)); axis([1,n,-1,1])
  xlabel(['i = ',num2str(i)]), I = I + 1;
end
subplot(2,2,1), text(12,1.2,'Right singular vectors V(:,i)'), pause
clf
I = 1;
for i=[n-2,n-5,n-8,n-11]
  subplot(2,2,I); plot(1:n,XX(:,i)), axis([1,n,-1,1]);
  xlabel(['i = ',num2str(i)]), I = I + 1;
end
subplot(2,2,1)
text(10,1.2,'Right generalized singular vectors XX(:,i)')
pause, clf

k_tsvd = 7; k_tgsvd = 6;
X_I = tsvd(U,s,V,b,1:k_tsvd);
X_L = tgsvd(UU,sm,XX,b,1:k_tgsvd);
subplot(2,1,1);
  plot(1:n,X_I,1:n,x,'x'), axis([1,n,0,1.2]), xlabel('L = I')
subplot(2,1,2);
plot(1:n,X_L,1:n,x,'x'), axis([1,n,0,1.2]), xlabel('L \neq I')
pause, clf

% Part 6.  No square integrable solution
% --------------------------------------
%
% In the last example there is no square integrable solution to
% the underlying integral equation (NB: no noise is added).
%
% Notice that the discrete Picard condition does not seem to
% be satisfied, which indicates trouble!

[A,b] = ursell(32); [U,s,V] = csvd(A);
picard(U,s,b); pause

% This concludes the demo.
echo off