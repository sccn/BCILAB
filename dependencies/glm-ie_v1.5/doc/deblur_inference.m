% Inference example doing sparse regression for image deblurring using glm-ie.
%
% We do nonblind deconvolution to recover the image u that has been blurred by
% an optical device with known point spread function f.
% The measurements are y = conv2(u,f) + e, e~N(0,s2).
%
% (c) by Hannes Nickisch, Philips Research, 2013 August 30
clear all, close all

load deblur
ut = u; clear u, f = f(4:end,:); f = f/sum(f(:));    % image u and blur filter f
su = size(ut); nu = prod(su); sf = size(f);                              % sizes

% blur operator and noisy measurement
X  = matConv2(f, su, 'circ'); sy = su;
s2 = 1e-5;                                                            % variance
y  = X*ut; y = y + sqrt(s2)*randn(size(y));                        % observation

% prior
B = matFD2(su);                           % total variation = finite diff matrix
tau = 15; t = 0;                                                   % prior scale

pot = @potLaplace;

%% Inference outer loop variances can be computed using Lanczos or Monte Carlo
% opts.outerMethod = 'lanczos';         % Lancos
% opts.outerMethod = 'sample';          % Monte Carlo
opts.outerMethod = 'factorial';       % factorial posterior approximation

opts.innerOutput = 1;
opts.outerOutput = 1;
opts.outerNiter = 4;                           % Number of outer loop iterations

opts.outerVarOpts = struct('MVM',170,                    ... % lanczos variances
                           'NSamples',20,'Ncg',20);           % sample variances
opts.innerMVM   = 100;        % Number of CG steps in inner loop
opts.innerVBpls = 'plsLBFGS'; % PLS algorithm for inner loop
[uinf,ga,b,z,zu,nlZ] = dli(X,y,s2,B,t,pot,tau,opts); uinf = reshape(uinf,su);
sdinf = sqrt(reshape(zu,su));

%% Estimation
opt.nMVM = 250; opt.output = 1;
u0 = zeros(nu,1);
uest = reshape( feval(opts.innerVBpls,u0,X,y,B,t,opt,s2,'penVB',pot,tau), su);

% accuracy
err_estimation = norm(uest(:)-ut(:))
err_inference  = norm(uinf(:)-ut(:))

% di-f1.png and di-f2.png
sz = [800,300]; figure('Position',[50,50,sz])
subplot('Position',[0/3,0,.32,.87])
  imagesc(abs(uest-ut),[0,.2]),  axis off, title('|true-mode|','FontSize',16)
subplot('Position',[1/3,0,.32,.87])
  imagesc(abs(uinf-ut),[0,.2]),  axis off, title('|true-mean|','FontSize',16)
  subplot('Position',[2/3,0,.31,.87])
  mn = min(sdinf(:)); mx = max(sdinf(:));
  clim = [mn, mx/2]; if mn>mx/2, clim = [mn,mx]; end
  imagesc(sdinf,clim), axis off, title('standard dev.','FontSize',16), colorbar
colormap(gray)