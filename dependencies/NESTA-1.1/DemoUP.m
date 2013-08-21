%{
    DemoUP.m

    A short example of using the unconstrained (aka "UP") version of NESTA

    Instead of solving 

        min || U x ||_1 s.t. || b - Ax||_2 <= delta

    the UP versions solve

        min   lambda*|| U x ||_1 + 1/2|| b - Ax ||_2^2

    The advantage is that we don't need to do any projections, so there's
    no restriction on A as in the regular version of NESTA.

    We do need to calculate ||A||^2


 Written by: Stephen Becker, Caltech
 Email: srbecker@acm.caltech.edu
 Ceated: Nov 2009

 NESTA Version 1.1

%}
    
clear all;clc; close all;
Setup_Nesta     %-- setup the path

N = 64^2;       %-- signal size 
M = floor(N/8); %-- number of measurements
K = floor(M/5); %-- number of nonzero elements in x
Sigma = 0.1;    %-- noise level

ExplicitA = false; %-- whether to use explicit matrix for A, or implicit function handles

fprintf('###############################################\n\n');
fprintf('NESTA_UP: Unconstrained version, N = %d\n\n',N);
fprintf('###############################################\n\n');

lambda = Sigma*sqrt(2*log(N));

%% GENERATE SIGNAL

disp('Creating Data');
randn('state',2009); rand('state',2009);
Amatrix = randn(M,N);    % no longer a projector
Omega = randperm(N); Omega = Omega(1:K);
x0 = zeros(N,1);
x0(Omega) = randn(K,1);
b = Amatrix*x0 + Sigma*randn(M,1);


if ExplicitA
    A = Amatrix; At=[];
    La = norm( A*A' );
else
    A = @(z) counter( @(x)Amatrix*x, z);
    At= @(z) counter( @(x)Amatrix'*x, z);
    La = my_normest( A, At, N )^2;
    % reset the counter
    nCalls = counter();
end

U = [];
Ut = [];


%% APPLY NESTA_UP v1 with continuation

disp('Applying NESTA_UP with continuation');
mu = 0.1*Sigma/La; %--- can be chosen to be small

opts = [];
opts.maxintiter = 6;
opts.TOlVar = 1e-6;
opts.verbose = 50;
opts.maxiter = 3000;
opts.U = U;
opts.Ut = Ut;
opts.stoptest = 1;  
counter();

[x_2,niter_2,resid_2,err_2,optsOut] = NESTA_UP(A,At,b,lambda,La,mu,opts);

N2 = counter();fprintf('Took %d calls\n',N2);

figure(2); clf;
stem( x0,'o','markersize',7); hold all
stem( x_2,'d','markerfacecolor',[0,.5,0],'markersize',5 );
legend('Original signal','l_1 reconstruction');
% We don't necessarily expect the error to be zero, unless k is very small
fprintf('relative l2 norm difference between original signal and l_1 reconstruction: %.2e\n',...
    norm( x0-x_2)/norm(x0) );
x_ref = x_2;

%% Using reweighting with U and Ut
for rwNo = 1:3
    fprintf('\n====== Reweighting %d =========\n\n',rwNo);
    xAbs = sort(abs(x_ref),'descend');
    cutoff = xAbs( round( .02*N) );
    u = 1./( abs(x_ref) + cutoff );
    U = spdiags( u, 0, N, N );
    Ut = U;
    normU = max(u);
    
    % give NESTA this extra information
    opts.U = U;
    opts.Ut = Ut;
    opts.normU = normU;
    opts.xplug = x_ref;
    opts.verbose = 150;
    
    [x_rw,niter_rw,resid_rw,err_rw,optsOut] = NESTA_UP(A,At,b,lambda,La,mu,opts);
    
    figure(3); clf;
    stem( x0,'o','markersize',7); hold all
    stem( x_rw,'d','markerfacecolor',[0,.5,0],'markersize',5 );
    legend('Original signal','l_1 reconstruction with reweighting');
    % We don't necessarily expect the error to be zero, unless k is very small
    fprintf('relative l2 norm difference between original signal and l_1 reconstruction: %.2e\n',...
        norm( x0-x_rw)/norm(x0) );
    
    x_ref = x_rw;
end

%% to get an even better answer, reproject
T = find( abs(x_ref) > 1e-2 );
At = Amatrix(:,T);
fprintf('\n====== Reprojectiong =========\n\n');
x_reproject = zeros(N,1);
x_reproject(T) = At\b;
    fprintf('relative l2 norm difference between original signal and l_1 reconstruction: %.2e\n',...
        norm( x0-x_reproject)/norm(x0) );