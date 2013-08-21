%  DemoComplex.m
%
%  This is a short script that tests NESTA with complex data
%
%  Both the signal and the operators may be complex
%
%  Runs several examples:
%   (1) A is partial unitary (A*A' = I)
%   (2) A is not partial unitary
%       (a) noiseless data, so provide Cholesky factors of A*A'
%       (b) noisy data, so provide svd(A)
%       (c) noisy data, so call NESTA_UP.m, which doesn't require A*A'=I
%
%  min ||x||_l1 s.t. ||b - A x||_2 < delta
%
%  Parameters to set : N : signal size (default: 64^2)
%                      M : number of measurements (default: N/8)
%                      K : number of nonzeros entries of x (default: M/5)
%                      Dyna : the dynamic range of the original image (in dB) (default : 60)
%                      Sigma : noise level (default : 0.1)
%                      Chg : if true, performs the change of variable x <- U(x) (default : true)
%
% Written by: Stephen Becker, Caltech
% Email: srbecker@acm.caltech.edu
% Created: Nov 2009
%
% NESTA Version 1.1
%   See also NESTA and Core_Nesterov

clear all;clc; close all;
Setup_Nesta     %-- setup the path
global COMPLEX
COMPLEX = true;  % this tells msp_signal to generate complex data

N = 64^2;       %-- signal size 
% N = 1024;
M = floor(N/8); %-- number of measurements
K = floor(M/5); %-- number of nonzero elements in x
Dyna = 60; 
Sigma = 0.1;    %-- noise level

ExplicitA = true; %-- whether to use explicit matrix for A, or implicit function handles

ExplicitU = false; %-- whether to use explicit matrix for U, or implicit function handles
                  %   (has no effect if Chg = false)
                  
Chg = true;     %-- if true, performs the change of variable x <- U(x)

fprintf('###############################################\n\n');
fprintf('NESTA: demo with Complex Data.  N = %d\n\n',N);
fprintf('###############################################\n\n');

n = sqrt(N);
delta = sqrt(M + 2*sqrt(2*M))*Sigma; 

%% GENERATE SIGNAL

disp('Creating Data');
U = @(z) ifft(z)*n;  % need to normalize to make it unitary
Ut = @(z) fft(z)/n;

[x0,b,A,At,supind,Omega]=msp_signal(N,M,K,U,Ut,Dyna,Sigma,0,1);
% x0 is the original signal; supind is the support of x0

U = @(z) counter(U,z);
Ut = @(z) counter(Ut,z);

disp(' ');
if Chg fprintf('With change of variable trick');end
disp(' ');


%% or, make A implicit but U explicit
if ~Chg
    
    S.type = '()';
    S.subs{1} = Omega;
    S.subs{2} = ':';
    upsample = @(x) subsasgn( zeros(N,size(x,2)),S,x);
    downsample =@(x) x(Omega,:);

    A_0 = @(x) downsample(fft(x))/n;   % use 1/n to make inv(A)=A'
    At_0 = @(x) ifft(upsample( x ))*n;
    A = @(z) counter(A_0,z);
    At = @(z) counter(At_0,z);
    U = []; Ut = [];
end
if ExplicitA
    Amatrix = A(eye(N));
    A = Amatrix; At=[];
end
if ExplicitU
    if Chg, disp('Change-of-variables trick not recommended in this case'); end
    Umatrix = U(eye(N));
    U = Umatrix; Ut=[];
end

%% APPLY NESTA v1 with continuation

disp('Applying NESTA with continuation');
mu = 0.1*Sigma; %--- can be chosen to be small

opts = [];
opts.maxintiter = 5;
opts.TOlVar = 1e-6;
opts.verbose = 50;
opts.maxiter = 3000;
opts.U = U;
opts.Ut = Ut;
opts.stoptest = 1;  
counter();

[x_2,niter_2,resid_2,err_2] = NESTA(A,At,b,mu,delta,opts);
if Chg
    % Undo change-of-variables trick if necessary
    if ExplicitU
        x_2 = U*x_2;
    else
        x_2 = U(x_2);
    end
end

N2 = counter();fprintf('Took %d calls\n',N2);
fprintf('relative l2 norm difference between original signal and l_1 reconstruction: %.2e\n',...
    norm( x0-x_2)/norm(x0) );
%% Plot results
figure(2); clf;
if COMPLEX
    subplot(2,1,1);
    stem( real(x0),'o','markersize',7); hold all
    stem( real(x_2),'d','markerfacecolor',[0,.5,0],'markersize',5 );
    title('Real part of data');
    subplot(2,1,2);
    stem( imag(x0),'o','markersize',7); hold all
    stem( imag(x_2),'d','markerfacecolor',[0,.5,0],'markersize',5 );
    title('Imaginary part of data');
else
    stem( x0,'o','markersize',7); hold all
    stem( x_2,'d','markerfacecolor',[0,.5,0],'markersize',5 );
end
legend('Original signal','l_1 reconstruction');
% We don't necessarily expect the error to be zero, unless k is very small
%   and Sigma is zero.
fprintf('relative l2 norm difference between original signal and l_1 reconstruction: %.2e\n',...
    norm( x0-x_2)/norm(x0) );


%% --- now, try an example with a matrix that isn't a projector
COMPLEX = true;  % this tells msp_signal to generate complex data
N = 1024;       %-- signal size 
M = floor(N/8); %-- number of measurements
K = floor(M/5); %-- number of nonzero elements in x
fprintf('###############################################\n\n');
fprintf('NESTA: demo with Complex Data, non-projector.  N = %d\n\n',N);
fprintf('###############################################\n\n');
n = sqrt(N);
% -- GENERATE SIGNAL
disp('Creating Data');
randn('state',2009); rand('state',2009);
A = randn(M,N) + COMPLEX*1i*randn(M,N);
% A = orth(A')';  % for testing, we can make A into a projector
x0 = zeros(N,1);
Omega = randperm(N); Omega = sort(Omega(1:K));
x0(Omega) = 20*(randn(K,1) + COMPLEX*1i*randn(K,1));
b=A*x0;
Sigma = .1;
bNoisy = b + Sigma/sqrt(1+COMPLEX)*( randn(size(b)) + COMPLEX*1i*randn(size(b)) );

% -- option 1, for noiseless data: since A is small, we can find inv(A*A') as follows:
opts = [];
mu = 1e-4;
delta = 0;
opts.verbose = 150;
opts.tolvar = 1e-7;
R = chol( A*A' );
opts.AAtinv = @(x) R\( R'\x );
disp('----- Using noiseless data ----- ');        
[x_3,niter_3] = NESTA(A,[],b,mu,delta,opts);
fprintf('relative l2 norm difference between original signal and l_1 reconstruction: %.2e\n',...
    norm( x0-x_3)/norm(x0) );

% -- option 2, for noisy data: since A is small, we can find svd(A)
opts = [];
mu = 1e-5;
% delta = norm( bNoisy - b)*.9;
delta = Sigma/sqrt(1+COMPLEX)*sqrt( M +sqrt(8*M) );
opts.verbose = 150;
opts.tolvar = 1e-6;
[U,S,V] = svd(A,'econ');
opts.USV.U=U;
opts.USV.S=S;
opts.USV.V=V;
disp('----- Using noisy data, providing SVD(A) ----- ');  
[x_4,niter_4,resid,outData,opts] = NESTA(A,[],bNoisy,mu,delta,opts);
fprintf('relative l2 norm difference between original signal and l_1 reconstruction: %.2e\n',...
    norm( x0-x_4)/norm(x0) );

fprintf('delta is %f, norm(Ax-b) is %f\n',delta,norm(A*x_4-bNoisy));
% -- option 3, for noisy data: use NESTA_UP.  No need to find svd(A) 
opts = [];
opts.verbose = 450;
La = norm(A*A');
opts.tolvar = 1e-6;
mu = 1e-5;
Lambda = Sigma/sqrt(1+COMPLEX)*sqrt( 2*log(N) )*7;
fprintf('----- Using noisy data, and NESTA_UP with Lambda = %f -----\n',Lambda); 
[x_5,niter_5] = NESTA_UP(A,[],bNoisy,Lambda,La,mu,opts);
fprintf('relative l2 norm difference between original signal and l_1 reconstruction: %.2e\n',...
    norm( x0-x_5)/norm(x0) );

fprintf('estimate of delta(lambda) is %f (delta is %f)\n',norm(A*x_5-bNoisy),delta);