function run_tests( LOGFILE )
if nargin < 1, LOGFILE = []; end

%{
    Code for section 6.2 from the NESTA paper

    This generates the data and save to a .mat file
    You can choose either analysis or synthesis,
    to generate the plots in fig. 6.1
    (to actually make the plots, see analyze_tests.m)

    Parameters: N = 2^20, 60 dB
    Either synthesis ("id = 2") or analysis ("id = 3")
%}

ANALYSIS = false;


snr = 60;

addpath ../NESTA

PLOT = false;   % make this "true" to see some plots
    
if ~isempty(LOGFILE)
    fp = fopen(LOGFILE,'w');
    LOG = true;
else
    fp = 1;
    LOG = false;
end

SAVE = true;    % make this "false" if you don't want to save to disk
if LOGFILE, SAVE = true; end

if ANALYSIS, ID = 3; else ID = 2; end
date = datestr(now,29);
FILE = sprintf('radar_id%02d_%s',ID,date);

N = 2^20;
%N = 2^12;
M = round( .30 * N);

% --- Setup sampling operator
seed = 2009 + 2;
randn('state',seed); rand('state',seed);

perm = randperm(N);
ROWS = sort(randsample(N,M));

param = [];
param.N = N;
param.M = M;
param.seed = seed;
param.perm = perm;
param.rows = ROWS;

% --- Randomly permute the columns ---
permute_cols = @(x) x(perm,:);
Sperm.type = '()'; Sperm.subs{1} = perm; Sperm.subs{2} = ':';
i_permute_cols = @(x) subsasgn( zeros(size(x)),Sperm,x);

% --- Randomly subsample the rows
downsample = @(x) x(ROWS,:);
S.type = '()'; S.subs{1} = ROWS; S.subs{2} = ':';
upsample = @(x) subsasgn( zeros(N,size(x,2)),S,x);

% -- And make the operator
sqrtN = 1/sqrt(N);
A = @(x) sqrtN*downsample( my_hadamard( permute_cols(x) ) );
At = @(x) sqrtN*i_permute_cols( my_hadamard( upsample(x) ) );

% make versions of these that do NOT increment the counter:
A_noCnt = @(x) sqrtN*downsample( hadamard( permute_cols(x) ) );
At_noCnt = @(x) sqrtN*i_permute_cols( hadamard( upsample(x) ) );

%% Setup a dictionary
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
%% test with a realistic signal
% N = 2^18;
[x,window2,window3,f2,f3,t,FS] = radar_signal(N); x= x.';
b_exact = A_noCnt(x);
param.x = x;
param.window2 = window2;
param.window3 = window3;
param.f2 = f2;
param.f3 = f3;
param.FS = FS;
if PLOT
    % plot( t, abs([window2;window3]) ); xlim([0,t(end)]);
    %
    % pwelch( x, [], [],2*N, FS );
    % spectrogram(x,N/16,N/32,N/32/4,FS);
    [Pxx,FreqList] = periodogram(x,[], N, FS );
    plot( FreqList, db(Pxx,'power') ); ylabel('dB');xlabel('Frequency (Hz)');
    y = [-160,-120];
    f1 = .8345*1e9;
    line( f1*[1,1], y , 'color','k');
    line( f2*[1,1], y , 'color','r');
    for ff = f3
        line( ff*[1,1], y , 'color','g');
    end
end
%%

sigma = makeSnr( b_exact, snr );
param.snr = snr;
param.sigma = sigma;
noise = sigma * randn(size(b_exact));
b = b_exact + noise;


% calculating the residual costs one iteration, so don't actually do
%   it for the final tests:

% residFcn = @(y) norm( A_noCnt(y) - b );   % costly
residFcn = @(y) 0;                          % free
errFcn_l2 = @(y) norm(y-x)/norm(x);         % cheap
errFcn_lInf = @(y) norm( x(:) - y(:),'inf');% cheap

%% Run Nesterov
Verbose = false;
if LOG, Verbose = false; end
MDec = false; TypeMin = 'L1';
delta = sqrt(M+2*sqrt(2*M)) * sigma;
THRESH = 5 * sigma* sqrt(M/N);
mu = sigma/5;
maxiter = 5000;
TolVar = 1e-5;
MaxIntIter = 5;

% record these choices:
param.ANALYSIS = ANALYSIS;
param.mu = mu;
param.maxiter = maxiter;
param.thresh = THRESH;
param.delta = delta;
param.TolVar = TolVar;
param.MaxIntIter = MaxIntIter;


errFcn = errFcn_l2;


fprintf(fp,'\n********* Starting Nesterov ***********\n');
if ANALYSIS, type='analysis'; else type='synthesis'; end
fprintf(fp,'Using %s, mu = %.2f*sigma, MaxIntIter %d\n',type,mu/sigma,MaxIntIter);
logData = []; tm = cputime;  my_hadamard(); counter();% zero-out counters
logFcn = @(y) [residFcn(y), errFcn_l2(y), ...
    errFcn_lInf(y) ];
titles = {'Residual','rel. l_2 error','l_{\infty}'};

fprintf(fp,'Beginning at %s\n', datestr(now));

if ANALYSIS
    % Make change of variable.  Do this before the row subsampling
    % let Hx = H*x, so x = H^-1* Hx, where H is the hadamard
    H = @(x) sqrtN * my_hadamard( permute_cols(x) );
    Hinv = @(x) sqrtN * i_permute_cols( my_hadamard( x ) );
    AA = @(x) downsample( x );
    AAt= @(x) upsample( x );

    U = @(x) psiT( Hinv( x ));
    Ut= @(x) H( psi( x ));


    [xk,niter,residuals,outputData] = Continuation_Nesterov(b,N,delta,mu,...
            AA,AAt,U,Ut,MaxIntIter,maxiter,MDec,TolVar,TypeMin,Verbose);
    xk = Hinv(xk);

else
   % Synthesis 

   AA = @(x) A( psi(x) );
   AAt= @(x) psiT( At( x ) );

   U = @(x) x;
   Ut= U;
   
   overcomplete = N_Gabor/N;
   mu = mu*overcomplete;
   TolVar = TolVar*overcomplete;
   param.TolVar = TolVar;
   param.mu = mu;

   % need to change N also
   [xk,niter,residuals,outputData] = Continuation_Nesterov(b,N_Gabor,delta,mu,...
       AA,AAt,U,Ut,MaxIntIter,maxiter,MDec,TolVar,TypeMin,Verbose);
   % xk is the dictionary coefficients
   xk = psi(xk); 
end
    
if niter >= maxiter -1
    fprintf(fp,'Warning!! Reached maxed # iterations');
end


param.err_l2 = errFcn_l2( xk );
param.err_linf = errFcn_lInf( xk );
param.niter = niter;

fprintf(fp,'l2 error: %.2e\tl_inf error: %.2e\n',param.err_l2,...
    param.err_linf );

tm = cputime - tm; fprintf(fp,'Elapsed time in Nesterov: %.2f sec\n',tm);
nCalls_gabor = counter();
nCalls_hadamard = my_hadamard();
param.nCalls_gabor = nCalls_gabor;
param.nCalls_hadamard = nCalls_hadamard;
param.tm = tm;
fprintf(fp,'# calls to hadamard measurement matrix: %d\n', nCalls_hadamard );
fprintf(fp,'# calls to gabor analysis dictionary: %d\n', nCalls_gabor );
fprintf(fp,'Finished at %s\n', datestr(now));
if SAVE
    save(FILE,'x','param','xk','logData');
    fprintf(fp,'Saved results to %s.mat\n', FILE);
end

if LOG
    fclose(fp);
elseif PLOT
    figure(1); clf;
    % plot the exact answer
    pgram = @(x) periodogram(x,[],N,FS);
    [Pxx,FreqList] = pgram(x);
    Pxx = db(Pxx,'power');
    plot( FreqList, Pxx,'linewidth',5 ); 
    ylabel('dB');xlabel('Frequency (Hz)');
    % plot the recovered answer
    [Pxx,FreqList] = pgram(xk);
    Pxx = db(Pxx,'power');
    hold on
    plot( FreqList, Pxx, '-r','linewidth',3 );
    
    % add some lines
%     y = [-160,-120];
    y = [ min( Pxx ), max( Pxx ) ];
    f1 = .8345*1e9;
    line( f1*[1,1], y , 'color','k');
    line( f2*[1,1], y , 'color','m');
    for ff = f3
        line( ff*[1,1], y , 'color','g');
    end
end
