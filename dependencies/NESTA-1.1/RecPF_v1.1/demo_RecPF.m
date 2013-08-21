function demo_RecPF

close all
clear all

addpath('utilities');
addpath('solver');
addpath('images');

% load images
Ims = cell(6,1);
for k = 1:6
    filename = ['ortho' int2str(k) '.jpg'];
    I = single(imread(filename))/255;
    if ndims(I) == 3
        Ims{k} = rgb2gray(I);
    end
end

% call tester
tester(Ims);

function tester(Ims)

%
% min aTV*TV(u) + aL1*||Phi*u||_1 + 0.5*||F_p*u - f_p||_2^2
%
aTV = 1.e-6;
aL1 = 0.e-6;

h = figure(1);
set(h,'units','normalized','outerposition',[0 .4 .9 .6]);

for nimg = 1:length(Ims)
    
    % load an image
    I = single(Ims{nimg});
    [m n] = size(I);
    N = m*n;
    
    % generate mask: two methods
    mth = 1;
    if mth == 1
        acc = 3.6;
        picks = selectPF(n/acc,n);
        K = length(picks);
    else
        Ls = 45;
        picks = fftshift(MRImask(n,Ls));
        picks = union(find(picks~=0),1);
        K = length(picks);
    end
    
    % image information display
    fprintf('Image %i: size %3i by %3i, N/M = %6.2f\n',nimg,m,n,N/K);
    snr(2*I,I);
    subplot(2,6,nimg); imshow(I,[]);
    title(['Original ' int2str(m) ' x ' int2str(n)]);
    drawnow;
    
    % add noise
    sigma = 0.01;
    noise = sigma*(randn(K,1) + sqrt(-1)*randn(K,1));
    
    % DWT and IDWT
    if exist('midwt','file')
        wav = daubcqf(2);
        Psi = @(x) midwt(x,wav);
        PsiT = @(x) mdwt(x,wav);
    elseif exist('wavedec2','file')
        Psi = @(x) Wavedb1Phi(x,1);
        PsiT = @(x) Wavedb1Phi(x,0);
    else
        Psi = @(x) x; PsiT = Psi;
    end
    
    
    % observation B
    FI = fft2(I);
    B = FI(picks) + noise;
    
    % run RecPF
    opts = [];
    t = cputime;
    U = RecPF(m,n,aTV,aL1,picks,B,PsiT,Psi,opts);
    t = cputime - t;
    
    % plot reconstructed images
    subplot(2,6,6+nimg);
    imshow(U,[]);
    title(sprintf('%4.2fdB, %4.2fs',snr(U,I),t));
    drawnow;
    
end
