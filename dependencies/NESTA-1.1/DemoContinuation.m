%  DemoContinuation.m
%
%  This is a short script that performs comparisons with or without
%  continuation when solving
%
%  min ||x||_l1 s.t. ||b - A x||_2 <= delta
%
%  Here A*A is assumed to be an orthogonal projector
%
%  Parameters to set : N : signal size (default: 64^2)
%                      M : number of measurements (default: N/8)
%                      K : number of nonzeros entries of x (default: M/5)
%                      Dyna : the dynamic range of the original image (in dB) (default : 60)
%                      Sigma : noise level (default : 0.1)
%                      Chg : if true, performs the change of variable x <- U(x) (default : true)
%
% Note: for continuation, call NESTA.m
%   without continuation, call Core_Nesterov.m
%   Not only do we recommend using NESTA with continuation (NESTA.m),
%   but this version also provides more error checking than Core_Nesterov.m
%
% Written by: Stephen Becker, Caltech
% Email: srbecker@acm.caltech.edu
% Created: May 2009
% Modified: May 2009 - Jerome Bobin
% Modified: Nov 2009 - Stephen Becker - showing how to "undo" after the
%   change-of-variables trick, and how to do reprojection
%
% NESTA Version 1.1
%   See also NESTA and Core_Nesterov

clear all;clc; close all;
Setup_Nesta     %-- setup the path

N = 64^2;       %-- signal size 
M = floor(N/8); %-- number of measurements
K = floor(M/5); %-- number of nonzero elements in x
Dyna = 60; 
Sigma = 0.1;    %-- noise level

ExplicitA = false; %-- whether to use explicit matrix for A, or implicit function handles

ExplicitU = false; %-- whether to use explicit matrix for U, or implicit function handles
                  %   (has no effect if Chg = false)

Chg = true;     %-- if true, performs the change of variable x <- U(x)

fprintf('###############################################\n\n');
fprintf('NESTA: experiments with or without continuation, N = %d\n\n',N);
fprintf('###############################################\n\n');

n = sqrt(N);
delta = sqrt(M + 2*sqrt(2*M))*Sigma; 

%% GENERATE SIGNAL

disp('Creating Data');
U = @(z) idct(z);
Ut = @(z) dct(z);
[x0,b,A,At,supind,Omega]=msp_signal(N,M,K,U,Ut,Dyna,Sigma,0,1);
% x0 is the original signal; supind is the support of x0

U = @(z) counter(U,z);
Ut = @(z) counter(Ut,z);

disp(' ');
if Chg fprintf('With change of variable trick');end
disp(' ');

% -- no change-of-variables trick:
if ~Chg
    
    S.type = '()';
    S.subs{1} = Omega;
    S.subs{2} = ':';
    upsample = @(x) subsasgn( zeros(N,size(x,2)),S,x);

%     A = @(x) A_Subset(dct(x),Omega);
    downsample = @(x) x(Omega,:);
    A = @(x) downsample(dct(x));
    At = @(x) idct(upsample( x ));
    A = @(z) counter(A,z);
    At = @(z) counter(At,z);
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


%% APPLY NESTA v1

disp('Applying Nesterov');
mu = 0.1*Sigma; %--- can be chosen to be small

opts = [];
opts.TOlVar = 1e-6;  % capitalization not important
opts.verbose = 50;
opts.maxiter = 3000;
opts.U = U;
opts.Ut = Ut;

[x_1,niter_1,resid_1,err_1] = Core_Nesterov(A,At,b,mu,delta,opts);

N1 = counter(); fprintf('Took %d calls\n',N1);
%%
% plot( resid_1(:,2) )  % see fmu as a function of iteration count

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

[x_2,niter_2,resid_2,err_2,optsOut] = NESTA(A,At,b,mu,delta,opts);

N2 = counter();fprintf('Took %d calls\n',N2);
%% Compare decay of residuals
figure(1);clf;plot( resid_1(:,2) ); hold all; plot( resid_2(:,2))  % see fmu as a function of iteration count
% plot( resid_2(:,1) ); line( [1,niter_nepa2],delta*[1,1]); % residuals

legend('NESTA','NESTA with continuation');
title('Value of f_\mu'); xlabel('iteration');

%% Display discrepancy between continuation and no-continuation versions
fprintf('w/o continuation: iter: %d,\tcalls: %g\t-- Per iteration : %g\n',niter_1,N1,round(N1/niter_1));
fprintf('w/  continuation: iter: %d,\tcalls: %g\t-- Per iteration : %g\n',niter_2,N2,round(N2/niter_2));
fprintf('Relative error between the two solutions, in dB:\t%g dB\n',20*log10(norm(x_1-x_2)/norm(x_2)));
fprintf('Relative error between the two solutions:\t\t %g\n',norm(x_1-x_2)/norm(x_2));

%% Undo change-of-variables trick if necessary
if Chg
    if isa(U,'function_handle')
        x_1 = U(x_1);
        x_2 = U(x_2);
    else
        x_1 = U*x_1;
        x_2 = U*x_2;
    end
end
figure(2); clf;
stem( x0,'o','markersize',7); hold all
stem( x_2,'d','markerfacecolor',[0,.5,0],'markersize',5 );
legend('Original signal','l_1 reconstruction');
% We don't necessarily expect the error to be zero, unless k is very small
fprintf('relative l2 norm difference between original signal and l_1 reconstruction:      %.2e\n',...
    norm( x0-x_2)/norm(x0) );
%% Optional: debias.
% Idea: l1 recovery is very good at finding the support of a signal,
% but the actual values are not always correct and can benefit from
% de-biasing, aka reprojection.  The idea is that we fix the support, T,
% and then solve a least-square problem.

xSorted = sort(abs(x_2),'descend');
cutoff = xSorted(M);
T = find( abs(x_2) >= max(cutoff,1) );

% If the measurement matrix Phi is explicit, this is as simple.
% Our measurement matrix wasn't explicit, but we can make it explicit:
if N < 10000   % don't try this with a huge matrix!
    Phi = zeros(M,length(T));
    e = zeros(N,1);
    if ~Chg && ExplicitA
        Phi = A;
    else
        for i = 1:length(T)
            e(T(i)) = 1;
            if Chg
                if ExplicitA
                    Phi(:,i) = A*Ut(e);
                else
                    Phi(:,i) = A(Ut(e));
                end
            else
                Phi(:,i) = A(e);
            end
            e(T(i)) = 0;
        end
    end
    x_debias = zeros(N,1);
    x_debias(T) = Phi\b;
    fprintf('relative l2 norm difference between original signal and l_1 reconstruction:      %.2e\n',...
        norm( x0-x_2)/norm(x0) );
    fprintf('relative l2 norm difference between original signal and debiased reconstruction: %.2e\n',...
        norm( x0-x_debias)/norm(x0) );
    figure(3); clf;
    stem( x0,'o','markersize',7 ); hold all
    stem( x_debias,'d','markerfacecolor',[0,.5,0],'markersize',5);
    title('After debiasing');
end
%%
% Another way to debias, if Phi is only given as a function handle
% Note: MATLAB's LSQR routine in versions prior to R2009a was faulty
%   If you have an older version, I suggest downloading Saunder's
%   implementation of LSQR from his Stanford website.
if verLessThan('matlab','7.8')
    disp('Sorry, LSQR implementation may be faulty; try upgrading Matlab');
else
    downsampleT = @(x) x(T);
    ST.type = '()';
    ST.subs{1} = T;
    ST.subs{2} = ':';
    upsampleT = @(x) subsasgn( zeros(N,size(x,2)),ST,x);
    
    if ExplicitA
        AA = @(x) A*x;
        AAt =@(x) A'*x;
    else
        AA = A; AAt = At;
    end
    
    if Chg
        Phi =  @(x) AA(Ut(upsampleT(x)));
        PhiT = @(x) downsampleT(U(AAt(x)));
    else
        Phi =  @(x) AA(upsampleT(x));
        PhiT = @(x) downsampleT(AAt(x));
    end
    F = @(x,transp) LSQRwrapper(x,transp,Phi,PhiT);
    TOL = 1e-6; MAXIT = 50; X0 = x_2(T);
    
    [x_debiasT,flag] = lsqr( F, b, TOL, MAXIT, [],[], X0 );
    
    x_debias = zeros(N,1);
    x_debias(T) = x_debiasT;
        fprintf('relative l2 norm difference between original signal and l_1 reconstruction:      %.2e\n',...
        norm( x0-x_2)/norm(x0) );
    fprintf('relative l2 norm difference between original signal and debiased reconstruction: %.2e\n',...
        norm( x0-x_debias)/norm(x0) );
        
end