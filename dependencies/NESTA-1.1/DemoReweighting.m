%  DemoReweighting.m
%
%  This is a short script that shows how to use NESTA
%   for reweighted l1
%
%  min ||Wx||_l1 s.t. ||b - A x||_2 < delta 
%
%   Where W is a diagonal operator that uses the previous estimate of x
%   to make ||x||_l1 more closely resemble ||x||_l0
%
%  Here A*A is assumed to be an orthogonal projector
%
% Written by: Stephen Becker, Caltech
% Email: srbecker@acm.caltech.edu
% Created: June 2009
%
% NESTA Version 1.1
%   See also NESTA and Core_Nesterov

clear all;clc;
Setup_Nesta         %-- setup the path for the solvers

fprintf('###############################################\n\n');
fprintf('NESTA: reweighting demo \n\n');
fprintf('###############################################\n\n');
%% LOAD A SIGNAL

randn('state',123); rand('state',123);
N = 400;
k = round(N/12);
omega = randsample(N,k);
x = zeros(N,1);
x(omega) = 1+randn(k,1);
x_exact = x;

M = round(N/4);  % number of measurements
seed = 1234;     % make it reproducible
randn('state',seed); rand('state',seed);

A = randn(M,N);  % Gaussian measurements
A = orth(A')';   % with orthogonal rows

%% RECONSTRUCT via NESTA

opts = [];
opts.Verbose = 0;
opts.tolvar = 1e-8;

b = A*x_exact;     % the data
At = [];
delta = 0;
muf = 1e-8;
[xk,niter,resid,outData] = NESTA(A,At,b,muf,delta,opts);

fprintf('No reweighting\n\tError is %.2e\n', norm(xk-x_exact)/norm(x_exact) );

figure(1); clf;
stem(x_exact); hold all
stem( xk )
title('no reweighting');

%% REWEIGHTING
for rw = 1:5
    fprintf('Reweighting %d time\n',rw );
%     opts.U = @(y) y./( abs(xk) + .1);
    % or...
    opts.U = spdiags( 1./(abs(xk)+.1),0,N,N);
    opts.Ut = opts.U;
    opts.normU = max( 1./(abs(xk)+.1) );
    opts.xplug = xk;  % use old solution as starting value
    [xk,niter,resid,outData,optsOut] = NESTA(A,At,b,muf,delta,opts);
    fprintf('\tError is %.2e\n', norm(xk-x_exact)/norm(x_exact) );
end
figure(2); clf;
stem(x_exact); hold all
stem( xk )
title('with reweighting');
