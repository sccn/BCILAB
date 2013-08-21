%  DemoAnalysis.m
%
%  This is a short script that shows how to use NESTA
%   along with an analysis/synthesis operator (e.g. frame operator)
%
%  min ||W^*x||_l1 s.t. ||b - A x||_2 < delta (ANALYSIS FORMULATION
%
%  min || a ||_l1 s.t. || b - A( W a ) ||_2 < delta (SYNTHESIS FORMULATION)
%
%  Here A*A is assumed to be an orthogonal projector
%
% Written by: Stephen Becker, Caltech
% Email: srbecker@acm.caltech.edu
% Created: May 2009
%
% NESTA Version 1.1
%   See also NESTA and Core_Nesterov

clear all;clc;
Setup_Nesta         %-- setup the path for the solvers
addpath Analysis    % -- setup path for the analysis stuff

fprintf('###############################################\n\n');
fprintf('NESTA: l1 analysis \n\n');
fprintf('###############################################\n\n');
%% LOAD A SIGNAL
HANDEL = load('handel');        % this comes with Matlab
x_exact = HANDEL.y;             % It's a clip from Handel's "Hallelujah"
FS = HANDEL.Fs;                 % The sampling frequency

PLAY = input('Play the original music clip? ','s');
if strcmpi(PLAY,'y') || strcmpi(PLAY,'yes')
    sound( x_exact );
end
N = length(x_exact);
% The Hadamard code used below requires inputs of 2^k
% So, padd the signal with zeros at the end
k = nextpow2(N);
x_exact = [x_exact; zeros(2^k-N,1) ];
N_short = N;
N = 2^k;

%% Take some measurements
% This example is sufficiently large that we'll need to take
% measurements via an efficient transform, e.g. a subsampled DCT
%   or a subsampled Hadamard.  Here, we'll use a subsampled
%   Hadamard, with randomly permuted columns.

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
S.type = '()'; S.subs{1} = ROWS; S.subs{2} = ':';
upsample = @(x) subsasgn( zeros(N,size(x,2)),S,x);

% -- And make the operator
sqrtN = 1/sqrt(N);
A = @(x) sqrtN*downsample( hadamard( permute_cols(x) ) );
At = @(x) sqrtN*i_permute_cols( hadamard( upsample(x) ) );

% -- Hadamard is a mex file, written by Peter Stobbe
% If there's no pre-built executable for your system, then
% you'll need to compile it
try
    y = hadamard(x_exact);
catch
    fprintf('Can''t find hadamard.%s; please compile yourself\n',mexext);
    fprintf('To do so, run:    mex hadamard.c\n');
    fprintf('If you have a recent version of the Signal Processing toolbox\n');
    fprintf(' you could try fwht instead (use the ''hadamard'' ordering)\n');
end

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

% semilogy( sort(abs(y),'descend' ) ); title('Coefficients for signal');

%% RECONSTRUCT via NESTA, using "ANALYSIS"


opts = [];
opts.U = psiT;
opts.Ut = psi;
opts.Verbose = 10;

b = A(x_exact);     % the data
delta = 1e-2;
muf = 1e-2;
disp('------- Reconstruction via Analysis -----------');
tic
[xk,niter,resid,outData] = NESTA(A,At,b,muf,delta,opts);
time.analysis = toc;

%% RECONSTRUCT via NESTA, using "SYNTHESIS"

opts = [];
opts.Verbose = 10;
AA = @(x) A( psi( x ) );
AAt =@(y) psiT( At( y ) );

delta = 1e-2;
mu = 1e-2;
disp('------- Reconstruction via Synthesis -----------');
tic
[coeff_k,niter,resid,outData] = NESTA(AA,AAt,b,muf,delta,opts);
xk2 = psi(coeff_k);
time.synthesis = toc;
%% PLOT IN TIME
figure(1); clf;
subplot(2,2,1);
coeff = psiT( xk );
semilogy( sort(abs(coeff),'descend') )
hold all
semilogy( sort(abs(coeff_k),'descend') );
semilogy( sort(abs(psiT(x_exact)),'descend') );
title('Dictionary coefficients of recovered signal');
xlabel('Sorted coefficients');
ylabel('Magnitude');
legend('recovered via analysis','recovered via synthesis',...
    'original signal');

% cutoff = 1e-2;
% line( [0,length(coeff)], cutoff*[1,1] );
% xk2 = psi( coeff.*( abs(coeff)>cutoff ) );

subplot(2,2,2);
plot( (1:N)/FS, xk );
xlabel('time (in seconds)');
ylim([-.8,.8]);
title('As recovered via analysis');

subplot(2,2,3);
plot( (1:N)/FS, xk2 );
xlabel('time (in seconds)');
ylim([-.8,.8]);
title('As recovered via synthesis');

subplot(2,2,4);
plot( (1:N)/FS, x_exact,'r' );
xlabel('time (in seconds)');
ylim([-.8,.8]);
title('Original signal');
%% PLOT IN FREQUENCY
figure(2); clf;
subplot(3,1,1);
periodogram( xk(1:N_short) ,[],[],FS);
title('As recovered via analysis');
xlim([0,1.5])

subplot(3,1,2);
periodogram( xk2(1:N_short) ,[],[],FS);
title('As recovered via synthesis');
xlim([0,1.5])

subplot(3,1,3);
periodogram( x_exact(1:N_short) ,[],[],FS);
title('Original signal');
xlim([0,1.5])
%% PLAY THE RECONSTRUCTIONS:
PLAY = input('Play the music clips? ','s');
if strcmpi(PLAY,'y') || strcmpi(PLAY,'yes')
    fprintf('The original recording:\n');
    sound( x_exact(1:N_short) );
    fprintf('As recovered via analysis (took %.1f seconds):\n',time.analysis);
    sound( xk(1:N_short) );
    fprintf('As recovered via synthesis (took %.1f seconds):\n',time.synthesis);
    sound( xk2(1:N_short) )
end
