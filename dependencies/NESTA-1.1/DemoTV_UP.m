%  DemoTV_UP.m
%
%  This is a short script that performs TV minimization recovery 
%   from partial Fourier measurements 
% 
%  min Lambda ||x||_TV +1/2||b - A x||_2
%
%  Here A*A is assumed to be an orthogonal projector
%
%  Parameters to set : 1) n : the size of the image to reconstruct (default: 128)
%                      2) Dyna : the dynamic range of the original image (in dB) (default : 40)
%                      3) Lambda : value of the Lagrange multiplier (default: 0.05)
%
%  Created entities : I : original image
%                     Xnesta : image recovered with Nesta
%                     NA_nesta = number of calls of A or A* for the method
%
% Written by: Jerome Bobin, Caltech
% Email: bobin@acm.caltech.edu
% Created: May 2009
%
% NESTA Version 1.1
%   See also NESTA_UP and Core_Nesterov_UP
clear; clc;
n = 128;   %--- The data are n*n images
Dyna = 40; %--- Dynamic range in dB
Lambda = 0.05; 
Setup_Nesta % -- add the solvers to the path
%% Setting up the experiment

N = n*n;
I=MakeRDSquares(n,7,Dyna);
x = reshape(I,n*n,1);
L = floor(55*(n/256)); %--- so as to get a compression ratio ~10
[M,Mh,mi,mhi] = LineMask(L,n);
OMEGA = mhi;
OMEGA = [OMEGA];
K = length(OMEGA);
A = @(z) A_fhp(z,OMEGA);
At = @(z) At_fhp(z,OMEGA,n);
b0 = A(x);
sigma = 0.1;
noise = sigma*randn(size(b0));
b = b0+noise;

fprintf('##############################################\n\n');
fprintf('NESTA: Total Variation minimization experiment\n\n');
fprintf('##############################################\n\n');
fprintf('Image size = %g x %g / Dynamic range : %g / Noise level : %g\n\n',n,n,Dyna,sigma);

%% Setting up the experiment for NESTA

fprintf(' NESTA starts\n\n');
U = @(z) z;
Ut = @(z) z;
mu = 0.2; %--- can be chosen to be small
opts = [];
opts.maxintiter = 5;
opts.TOlVar = 1e-4;
opts.verbose = 25;
opts.maxiter = 5000;
opts.U = U;
opts.Ut = Ut;
opts.stoptest = 1;  
opts.typemin = 'tv';
counter();
Ac = @(z) counter(A,z);
Atc = @(z) counter(At,z);
tic;
[x_nesta,niter,resid,err] = NESTA_UP(Ac,Atc,b,Lambda,1,mu,opts);
t.NESTA = toc;
NA_nesta = counter();
tvnesta = calctv(n,n,reshape(x_nesta,n,n));
Xnesta = reshape(x_nesta,n,n);


% -- Print results

fprintf('\nNESTA: TV: %g -- ||b-Ax|| = %g -- # calls : %g\n\n',tvnesta,norm(b - A(Xnesta(:))),NA_nesta);
fprintf('relative l2 error compared to original signal: %.2e\n',norm(Xnesta-I,'fro')/norm(I,'fro') );
%% Plot results
figure(1); clf;

subplot(1,2,1);
image( I );  title('Original image');
map = colormap(hsv);
map(1,:) = [1,1,1];
colormap(map);
% colorbar('westoutside');
axis square; axis off

subplot(1,2,2);
image(Xnesta); title(sprintf(...
    'Reconstruction via NESTA\n%.1f sec',t.NESTA));
axis square; axis off