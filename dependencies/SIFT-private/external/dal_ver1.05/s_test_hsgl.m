% s_test_hsgl - a sample program for sparse connectivity inference.
%
% This program infers a sparsely connected multivariate AR
% model assuming that the sources are independent and super
% Gaussian. We use the hyperbolic secant likelihood.
%
% In our IEEE TBME paper, we also estimate mixing matrix
% through an EM algorithm. This is omitted here for the sake
% of simplicity. Thus, all state variables are directly measured.
%
% Reference:
% Modeling sparse connectivity between underlying brain sources
% for EEG/MEG. Stefan Haufe, Ryota Tomioka, Guido Nolte,
% Klaus-Robert Mueller, and Motoaki Kawanabe, IEEE
% Trans. Biomed. Eng. 57(8), pp. 1954-1963, 2010.
% 
% Copyright(c) 2009-2011 Ryota Tomioka
%              2009      Stefan Haufe
% This software is distributed under the MIT license. See license.txt


M=20;
N=5000;
k=10;
P=3;

H0=randmvar(M,P,k);
Z=mvarfilter(H0, 2/pi*log(tan(pi*rand(N,M)/2)))';
lmd=150;

X = [];
for jj = 1:P
  X = cat(3, X, Z(:, P+1-jj:end-jj)'); 
end
X = reshape(permute(X, [1 3 2]), N-P, M*P);
A = kron(speye(M), X);

%% Indices for the diagonal elements
indsB = vec(repmat((1:P*(M+1):P*M^2), P, 1) + repmat((0:(P-1))', 1, M))';

%% Indices for the off-diagonal elements
indsA = setdiff_bc(1:P*M^2, indsB);

%% Design corresponding to the diagonal elements (self connection)
B = A(:, indsB);

%% Design corresponding to the off-diagonal elements (connection to others)
A = A(:, indsA);

%% Target
Y = vec(Z(:, P+1:end)');

% Group penalize AR coefficients with different time-lags
% between sources (M*(M-1) groups of size P).
% Self connections are put into one group (size P*M).
[H, U, status] = dalhsgl(zeros(P*M^2,1),[], [A B],[], ...
                         Y, lmd, 'blks', [P*ones(1,M*(M-1)),P*M]);

H2 = zeros(P,M,M);
H2(indsA) = vec(H(1:(P*M*(M-1))));       % non-diagonal elements
H2(indsB) = vec(H(((P*M*(M-1))+1):end)); % diagonal elements

H = permute(full(H2), [3 2 1]);

figure;
for jj=1:P
  subplot(2,P,jj); imagesc(H0(:,:,jj))              
  subplot(2,P,jj+P); imagesc(H(:,:,jj));
end

subplot(2,P,round(P/2));
title('True');

subplot(2,P,round(P/2)+P);
title(sprintf('Estimated (lambda=%g)',lmd));

