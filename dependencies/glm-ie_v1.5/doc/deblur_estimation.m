% Estimation example doing sparse regression for image deblurring using glm-ie.
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
W = matWav(su); D = matFD2(su);       % construct wavelet and finite diff matrix
B = [W; D];                                          % sparsity transform matrix
tau = 15;  t = 0;                                                  % prior scale
lam = s2*tau;
dof = 2; pen = @(s) penLogSmooth(s,dof);

% optimisation parameters
opt.nMVM = 50; opt.output = 1; 
opt.LBFGSnonneg = 1;                        % non-negativity condition for LBFGS
u0 = zeros(nu,1);

%% Estimation: apply different PLS schemes
plsList = {'TN','CGBT','CG','BB'}; % also LBFGS is possible
for i=1:length(plsList)
  fprintf('PLS optimisation using %s.\n',plsList{i})
  pls = ['pls',plsList{i}];
  tic, [u{i},phi(i)] = feval(pls,u0,X,y,B,t,opt,lam,pen); tt(i) = toc;
  d(i) = norm(u{i}-ut(:));
end

fprintf('objective function values\n  ')
for i=1:length(plsList), fprintf('%s %1.4e, ',plsList{i},phi(i)), end
fprintf('\b\b\n')

fprintf('running times\n  ')
for i=1:length(plsList), fprintf('%s %1.2fs, ',plsList{i},tt(i)), end
fprintf('\b\b\n')

fprintf('accuracies in %%\n  ')
for i=1:length(plsList), fprintf('%s %1.2f%%, ',plsList{i},d(i)), end
fprintf('\b\b\n')

% de-f1.png
sz = [800,300]; figure('Position',[50,50,sz])
subplot('Position',[0/3,0,.32,.87])
  imagesc(reshape(y,sy)),      axis off, title('blurry',   'FontSize',16)
subplot('Position',[1/3,0,.32,.87])
  imagesc(reshape(ut,su)),     axis off, title('original', 'FontSize',16)
  subplot('Position',[2/3,0,.32,.87])
  imagesc(reshape(u{end},su)), axis off, title('deblurred','FontSize',16)
colormap(gray)