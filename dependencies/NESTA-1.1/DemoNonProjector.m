%{ 
DemoNonProjector.m, based off DemoAnalysis.m

 This is a short script that shows how to use NESTA
  along with an analysis/synthesis operator (e.g. frame operator)

  In addition, we assume that the measurement operator A is
  has much smaller dimensions than the frame operator W
  so that the dominant cost is calculating W(y) and W^*(x),
  and hence we can afford to compute one O(n^3) computation of A
    (e.g. chol(A), svd(A), inv(A*A'), ...)
  Or, we can solve for inv(A*A') via a Krylov subspace method.

This is the problem we are solving:

 min ||W^*x||_l1 s.t. ||b - A x||_2 < delta (ANALYSIS FORMULATION)

 A'*A is NO LONGER assumed to be an orthogonal projector
 i.e. A*A' is not necessarily the identity

There are 4 ways we can handle this:

Case 1: noiseless data (delta=0)
    - reorthogonalize A (and b) using QR
    This works for delta > 0, but the noise model is now
    slightly different, so the constraints may be off.

Case 2: noiseless data (delta=0)
    - pass in a routine to calculate inv(A*A')
    This can be done with either the "inv" function or a SVD factorization,
    but we recommend the Cholesky factorization.  See code below
    for how to do this.  You could explicitly calculate inv(AA'), but this 
    is slower and less numerically stable than the Cholesky factorization.
    Calculating the SVD works ( if U*diag(s)*V' = A, then inv(AA')=U*diag(s.^(-2))*U'  ),
    but computing svd(A) is much slower than chol(A*A'), so this is not recommended
    either.

Case 3: noisy data (delta>0)
    - pass in the SVD factors of A.  This is slower than case 2.
    If you have noisy data, this is your only option!

Case 4: noiseless data (delta=0)
    - when A is large but has a fast transform, we cab calculate 
    inv(A*A') via a Krylov Subspace method, e.g. conjugate-gradient.


Note: this script requires the Signal Toolbox for some of the plotting
tools.  If this is not available, then the reconstruction will
still work but some plots will not show up.

Written by: Stephen Becker, Caltech
Email: srbecker@acm.caltech.edu
Created: Oct 2009
Modified: Nov 2009

NESTA Version 1.1
  See also NESTA and Core_Nesterov
%}

clear all;clc;
Setup_Nesta         %-- setup the path for the solvers
addpath Analysis    % -- setup path for the analysis stuff

fprintf('###############################################\n\n');
fprintf('NESTA: l1 analysis, non projector case \n\n');
fprintf('###############################################\n\n');


HAS_SIGNAL_TOOLBOX = license('test','signal_toolbox');
%% LOAD A SIGNAL
HANDEL = load('handel');        % this comes with Matlab
x_exact = HANDEL.y;             % It's a clip from Handel's "Hallelujah"
FS = HANDEL.Fs;                 % The sampling frequency

% for non-projector case, we need to do an SVD, so make the size more
% manageable:
LENGTH = 2^13;  % you can change this

offset = 1500;
x_exact = x_exact(offset+(1:LENGTH));

PLAY = input('Play the original music clip? ','s');
if strcmpi(PLAY,'y') || strcmpi(PLAY,'yes')
    sound( x_exact );
end
N = length(x_exact);
% The Hadamard code used below requires inputs of 2^k
% So, pad the signal with zeros at the end
k = nextpow2(N);
x_exact = [x_exact; zeros(2^k-N,1) ];
N_short = N;
N = 2^k;

%% Take some measurements
% Here, we'll use a subsampled
% Hadamard, with randomly permuted columns.

M = round(N/6);  % number of measurements
seed = 1234;     % make it reproducible
randn('state',seed); rand('state',seed);
fprintf('Subsampling the signal by a factor of %d\n', round(N/M) );

% --- Randomly permute the columns ---
perm = randperm(N);      % pick a permuation
permute_cols = @(x) x(perm,:);
Sperm.type = '()'; Sperm.subs{1} = perm; Sperm.subs{2} = ':';
i_permute_cols = @(x) subsasgn( zeros(size(x)),Sperm,x);

% --- Randomly subsample the rows
ROWS = randperm(N); ROWS = ROWS(1:M);   % pick some rows
downsample = @(x) x(ROWS,:);
SS.type = '()'; SS.subs{1} = ROWS; SS.subs{2} = ':';
upsample = @(x) subsasgn( zeros(N,size(x,2)),SS,x);

% -- And make the operator (which is not a projection)
Amatrix = randn(M,N);
% if we wanted to matrix this orthonormal, just do:
% Amatrix = orth(Amatrix')';

A = @(x) Amatrix*x;
At= @(x) Amatrix'*x;

%% SETUP AN ANALYSIS OPERATOR
% The PsiWFF and PsiTransposeWFF code is a Gabor frame
% (i.e. a short-time Fourier transform)
% written by Peter Stobbe
% 
% PsiWFF is the synthesis operator, and acts on coefficients
% PsiTransposeWFF is the analysis operator, and acts on signals
gMax = 0;
gLevels = 1;
tRedundancy = 1;
fRedundancy = 0;
gWindow = 'isine';

gabor.gMax = gMax; gabor.gLevels = gLevels;
gabor.tRedundancy = tRedundancy; gabor.fRedundancy = fRedundancy;
param.gabor = gabor;

logN = log2(N);
psi_A = @(y) PsiWFF(y,gWindow,logN,logN-gLevels,logN+gMax,tRedundancy,fRedundancy);
psiT_A = @(x) PsiTransposeWFF(x,gWindow,logN,logN-gLevels,logN+gMax,tRedundancy,fRedundancy);

N_Gabor = length( psiT_A( ones(N,1) ) );
param.N_Gabor = N_Gabor;

psi = @(y) counter( psi_A, y);
psiT= @(x) counter( psiT_A,x);

%
x = x_exact;
y = psiT(x);
x2 = psi(y);
% The frame is tight, so the psuedo-inverse is just the transpose:
fprintf('Error in PSI(PSI^* x ) - x is %e\n', norm(x-x2));
%% How compressible is this signal?
% semilogy( sort(abs(y),'descend' ) ); title('Coefficients for signal');
% We can make a sparse version of this signal by truncating
% the coefficients to keep only the ones that constitute 90% of the power
cutoff = .90*(norm(y)^2);
ySorted = sort(abs(y),'descend');
threshold = find( cumsum( ySorted.^2 ) > cutoff );
fprintf('Can represent 90%% of signal power with only %.1f%% nonzeros\n', threshold(1)/length(y)*100 );
threshold = ySorted(threshold(1));

% y_sparse = y.*( abs(y) > threshold );
% x_sparse = psi(y_sparse);
% fprintf('Error in x_sparse - x is %e\n', norm(x-x_sparse)/norm(x));
%% RECONSTRUCT via NESTA, using "ANALYSIS"

opts = [];
opts.U = psiT;
opts.Ut = psi;
opts.Verbose = 10;
opts.TolVar =1e-4;
b = A(x_exact);     % the data
% b = A(x_sparse);
muf = 1e-3;
disp(' ');
disp('Here are four ways to deal with A when A''*A is not a projector');
disp('Type 1: orthognalize A.  Noiseless data only');
disp('Type 2: use Cholesky decomposition to find inv(AA'').  Noiseless data only');
disp('Type 3: use SVD decomposition of A.  Noisy data is OK');
disp('Type 4: use CG or other Krylov method to find inv(AA''). Noiseless data only');
disp(' ');
% For the non-projector case, we need more info:

Type = 2;  % change this as desired

fprintf('Current type selected: type %d\n',Type);
tic; disp('Precalculation...');
AA = A; AAt = At; bb = b;
switch Type
    case 0
        % Assume A has orthogonal rows. Verify the assumption:
        if norm( A( At(b) ) - b ) < 1e-10
            disp('Measurement matrix has orthogonal rows');
        else
            error('Measurement matrix does not have orthogonal rows');
        end
        delta = 1e-2;
    case 1
        % orthogonalize A with a QR.  This works for noiseless data
        % It's meant for the case where A is a matrix, not a function
        delta = 0;
        [Q,R] = qr(Amatrix',0);  % A = (Q*R)'= R'*Q', R is triangular
        bb = (R')\b;
        AA = @(x) (R')\A(x);  % this is now Q'
        AAt= @(x) At(R\x);    % this is now Q
    case 2
        % find a cholesky decomposition of A. Noiseless data only
        delta = 0;
        R = chol( Amatrix*Amatrix' );
        opts.AAtinv = @(x) R\( R'\x );
    case 3
        % Find the SVD of A.  Noisy data are OK.
        if N_short <= 4096
            delta = 1e-5;
            [U,S,V] = svd(Amatrix,'econ');
            opts.USV.U=U;
            opts.USV.S=S;
            opts.USV.V=V;
        else
            disp('This method not recommended when A is large');
            disp('Using Type 4 instead');
            Type = 4;
        end
    case 4
        % Use a Krylov Subspace method to solve for inv(AA')
        % Noiseless data only.  Here, we'll use Conjugate-Gradients
        % (call CGwrapper, which calls MATLAB's "pcg" function)
        delta = 0;
        A_function = @(x) A(At(x));
        cg_tol = 1e-6; cg_maxit = 40;
        CGwrapper(); % (first, zero-out the CGwrapper counters)
        opts.AAtinv = @(b) CGwrapper(A_function,b,cg_tol,cg_maxit);
    case 5
        % This is an example of what NOT TO DO. 
        disp('This option definitely not recommended');
        delta = 0;  % noiseless data only
        opts.AAtinv = @(x) inv( Amatrix*Amatrix' )*x;
    case 6
        % This is much better than case 5 because the inverse is only
        % computed once, but it's not as numericaly robust (or as fast)
        % as the Cholesky decomposition.  So, use case 2 instead of this.
        disp('This option not recommended');
        delta = 0;
        Ainv = inv( Amatrix*Amatrix' );
        opts.AAtinv = @(x) Ainv*x;  
    otherwise
        disp('Invalid option for "Type" variable');
end
disp('Finished'); toc

disp('------- Reconstruction via Analysis -----------');
tic
[xk,niter,resid,outData] = NESTA(AA,AAt,bb,muf,delta,opts);
time.analysis = toc;
fprintf('l2 error: %.2e\n',norm(xk - x_exact )/norm(x_exact));
%% PLOT IN TIME
figure(1); clf;
subplot(3,1,1);
coeff = psiT( xk );
semilogy( sort(abs(coeff),'descend') )
hold all
semilogy( sort(abs(psiT(x_exact)),'descend') );
title('Dictionary coefficients of recovered signal');
xlabel('Sorted coefficients');
ylabel('Magnitude');
legend('recovered via analysis',...
    'original signal');

% cutoff = 1e-2;
% line( [0,length(coeff)], cutoff*[1,1] );
% xk2 = psi( coeff.*( abs(coeff)>cutoff ) );

subplot(3,1,2);
plot( (1:N)/FS, xk );
xlabel('time (in seconds)');
ylim([-.8,.8]);
title('As recovered via analysis');

subplot(3,1,3);
plot( (1:N)/FS, x_exact,'r' );
xlabel('time (in seconds)');
ylim([-.8,.8]);
title('Original signal');
%% PLOT IN FREQUENCY
if HAS_SIGNAL_TOOLBOX
    window = blackmanharris( max( 512, round(N_short/4) ));
    window=[window(1:end/2);ones(N_short-length(window),1);window(end/2+1:end)];
    figure(2); clf;
    subplot(2,1,1);
    periodogram( xk(1:N_short) ,window,[],FS);
    title('As recovered via analysis');
    xlim([0,1.5])
    
    subplot(2,1,2);
    periodogram( x_exact(1:N_short) ,window,[],FS);
    title('Original signal');
    xlim([0,1.5])
    
    % See how the pitch changes in time
%     figure(3);
%     spectrogram( xk(1:N_short),[],[],[],FS );
%     xlim([200,2000]);

% line( .2626*[1,1],[-90,-20]); % middle C
% line( .440*[1,1],[-90,-20]);  % concert A
end
%% PLAY THE RECONSTRUCTIONS:
PLAY = input('Play the music clips? ','s');
if strcmpi(PLAY,'y') || strcmpi(PLAY,'yes')
    fprintf('The original recording:\n');
    sound( x_exact(1:N_short) );
    fprintf('As recovered via analysis (took %.1f seconds):\n',time.analysis);
    sound( xk(1:N_short) );
end
