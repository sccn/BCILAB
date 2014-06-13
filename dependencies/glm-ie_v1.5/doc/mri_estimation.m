% Estimation example doing MRI reconstruction using glm-ie.
%
% (c) by Hannes Nickisch, Philips Research, 2013 August 30
clear all, close all

type_str = {'line', 'mask', 'spiral'};
type = 3;
complex = [2,2];

if type==1                                                        % line example 
  load mri_line
  mfull = 1;
  X = matFFT2line(size(utrue),id,complex,mfull);
  nMVM = 75;
elseif type==2                                                    % mask example
  load mri_mask
  mfull = 1;
  X = matFFTNmask(mask,complex,mfull);
  nMVM = 100;
else                                                            % spiral example 
  load mri_spiral, type = 3;
  X = matFFT2nu(size(utrue),k,complex);
  nMVM = 100;
end

y = cx2re(y(:)); su = size(utrue); nu = numel(utrue);
fprintf('MRI reconstruction from %s Fourier acquisition\n', type_str{type})
uls = reshape(re2cx(X'*y),su);                    % least squares reconstruction
B = [matWav(su,[],[],complex); matFD2(su,'circ',complex)]; % wavelet + fin-diff.
pen = @penAbs; t = 0;                                        % Laplace potential

s2  = 4e-4;                                                     % noise variance
opt.nMVM = nMVM; opt.output = 1;                       % optimisation parameters

if type<3                   % plsSB is only fast for matFFT2line and matFFTNmask
  pls1 = 'SB'; opt.SBouter = 75; opt.SBinner = 3; opt.SBga = 5;
else
  pls1 = 'CGBT'; 
  uls = uls*norm(utrue(:))/norm(uls(:));                          % adjust scale
end

plsList = {pls1,'BB'};             % Estimation: apply two different PLS schemes
for i=1:length(plsList)
  fprintf('PLS optimisation using %s.\n',plsList{i})
  pls = ['pls',plsList{i}];
  tic, [u{i},phi(i)] = feval(pls,cx2re(0*uls(:)),X,y,B,t,opt,s2,pen);
  tt(i) = toc;
  u{i} = abs(re2cx(u{i}));
  d(i) = 100*norm(u{i}-utrue(:))/norm(utrue(:));
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

% me-f1.png, me-f2.png, me-f3.png
sz = [800,300]; figure('Position',[50,50,sz])
subplot('Position',[0/3,0,.32,.87])
  imagesc(reshape(abs(uls),su)), axis off
  err = norm(abs(utrue(:))-abs(uls(:)))/norm(abs(utrue(:)));
  str = sprintf('u_{LS}, err=%1.3f',err);
  title(str,'FontSize',16)
subplot('Position',[1/3,0,.32,.87])
  if type==1
    idx = zeros(su(1),1); idx(id) = 1; imagesc(idx)
    r = 100*sum(numel(id)*su(2))/nu;
  elseif type==2
    imagesc(mask)
    r = 100*sum(mask(:))/nu;
  else
    plot(k,'.-')
    r = 100*numel(k)/nu;
  end
  axis off, title(sprintf('Pattern; %1.2f%% covered',r),'FontSize',16)

subplot('Position',[2/3,0,.32,.87])
  imagesc(reshape(u{end},su)), axis off
  err = norm(abs(utrue(:))-abs(u{end}(:)))/norm(abs(utrue(:)));
  str = sprintf('u_{MAP}, err=%1.3f',err);
  title(str,'FontSize',16)
colormap(bone)