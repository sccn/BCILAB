%  DemoTV.m
%
%  This is a short script that performs comparisons between NESTA and RecPF
%  in the scope of TV minimization recovery from partial Fourier measurements 
% 
%  min ||x||_TV s.t. ||b - A x||_2 <= epsilon
%
%  Here A*A is assumed to be an orthogonal projector
%
%  Parameters to set : 1) n : the size of the image to reconstruct (default: 128)
%                      2) Dyna : the dynamic range of the original image (in dB) (default : 40)
%
%  Created entities : I : original image
%                     Xnesta : image recovered with Nesta
%                     Xrpf_default : image recovered with RecPF using the default parameters
%                     Xrpf : image recovered with RecPF using the tuned parameters
%                     NA_X = number of calls of A or A* for the method (X = nesta, nesta_chg, rpf, rpf_default)
%
% Written by: Jerome Bobin, Caltech
% Email: bobin@acm.caltech.edu
% Created: May 2009
%
% NESTA Version 1.1
%   See also NESTA and Core_Nesterov

n = 128;   %--- The data are n*n images
Dyna = 40; %--- Dynamic range in dB
lambda = 0.05; %--- Lagrange multiplier in lambda ||x||_TV + 1/2||b-Ax|| 

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

%% Running RecPF for a fixed value of lambda to get epsilon(lambda)

fprintf('1) run RecPF to get epsilon = || b - A x_lambda || as a function of lambda (= %g)\n\n',lambda);
Psi = @(x) x; PsiT = Psi;
ropts = [];
ropts.mit_inn = 10000;ropts.mit_out=10;ropts.tol_rel_inn=1e-5;ropts.tol_rel_out=1e-2;
B = fft2(reshape(At(b),n,n)); %--- observation in the Fourier domain
OMEGA = [1;OMEGA];
B = B(OMEGA);
B = B(:);
X_rpf = RecPF(n,n,lambda,0,OMEGA,B,PsiT,Psi,ropts,[]);

%% Setting up the experiment for Nesta

fprintf('2) run NESTA\n\n');
U = @(z) z;
Ut = @(z) z;
mu = 0.2; %--- can be chosen to be small
opts = [];
opts.maxintiter = 5;
opts.TOlVar = 1e-5;
opts.verbose = 0;
opts.maxiter = 5000;
opts.U = U;
opts.Ut = Ut;
opts.stoptest = 1;  
opts.typemin = 'tv';
delta = norm(b - A(X_rpf(:)));
counter();
Ac = @(z) counter(A,z);
Atc = @(z) counter(At,z);
tic;
[x_nesta,niter,resid,err] = NESTA(Ac,Atc,b,mu,delta,opts);
t.NESTA = toc;
NA_nesta = counter();
NA_nesta_chg = 2*niter; %-- with the change of variable simplification, NA = 2*niter
tvnesta = calctv(n,n,reshape(x_nesta,n,n));
Xnesta = reshape(x_nesta,n,n);


%% Running RecPF for a fixed value of lambda and the new stopping criterion
%%% // default parameters
fprintf('3) run RecPF (default parameters) with the stopping criterion\n\n');

ropts = [];
F = @(z) fft2(z);iF = @(z) ifft2(z);
Fc = @(z) counter(F,z);iFc = @(z) counter(iF,z);
tic;
x_rpf_default = RecPF_Modified(n,n,lambda,0,OMEGA,B,PsiT,Psi,ropts,[],Fc,iFc,tvnesta,delta,b,Ac);
t.RecPF1 = toc;
NA_rpf_default = counter();
tvrpf_default = calctv(n,n,x_rpf_default);
Xrpf_default = x_rpf_default;

%% Running RecPF for a fixed value of lambda and the new stopping criterion
%%% // tuned parameters

fprintf('4) run RecPF with the modified stopping criterion and the tuned parameters\n\n');

ropts = [];
ropts.mit_inn = 10000;ropts.mit_out=10;ropts.tol_rel_out=1e-2;
ropts.tol_rel_inn=1e-5; % to refine the TV approximation
tic
x_rpf = RecPF_Modified(n,n,lambda,0,OMEGA,B,PsiT,Psi,ropts,[],Fc,iFc,tvnesta,delta,b,Ac);
t.RecPF2 = toc;
NA_rpf = counter();
tvrpf = calctv(n,n,x_rpf);
Xrpf = x_rpf;

%% Print results

fprintf('NESTA: TV: %g -- ||b-Ax|| = %g -- # calls : %g -- with chg of var. : %g\n',tvnesta,norm(b - A(Xnesta(:))),NA_nesta,NA_nesta_chg);
fprintf('RecPF, def. param. : TV: %g -- ||b-Ax|| = %g -- # calls : %g\n',tvrpf_default,norm(b - A(Xrpf_default(:))),NA_rpf_default);
fprintf('RecPF, tuned param. : TV: %g -- ||b-Ax|| = %g -- # calls : %g\n',tvrpf,norm(b - A(Xrpf(:))),NA_rpf);

%% Plot results
figure(1); clf;

subplot(2,2,1);
image( I );  title('Original image');
map = colormap(hsv);
map(1,:) = [1,1,1];
colormap(map);
% colorbar('westoutside');
axis square; axis off

subplot(2,2,2);
image(Xnesta); title(sprintf(...
    'Reconstruction via NESTA\n%.1f sec',t.NESTA));
axis square; axis off

subplot(2,2,3);
image( Xrpf_default ); title(sprintf(...
    'Reconstruction via RecPF\ndefault parameters\n%.1f sec',t.RecPF1));
axis square; axis off

subplot(2,2,4);
image( Xrpf ); title(sprintf(...
    'Reconstruction via RecPF\ntuned parameters\n%.1f sec',t.RecPF2));
axis square; axis off

%% To see what the measurements looked like (in Fourier domain)
figure(2); clf; set(gcf,'name','Measurements','toolbar','none','menubar','none');
imshow(fftshift(fftshift(M,1),2))
title(sprintf('Measurements in Fourier Domain\nsampling only %.1f%% of 2D Fourier Coefficients',...
    100*length(OMEGA)/N ) );
% This usually renders small, so make it bigger by stretching width and
% height
rect = get(gcf,'position');  
set(gcf,'position',[rect(1),rect(2),1.5*rect(3),1.5*rect(4)] );