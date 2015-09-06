function [v]=tvarsim(w,A,C,n_ntr,ndisc,beta,alpha,distribution)
%  v=tvarsim(w,A,C,n) simulates n time steps of the AR(p) process
%
%     v(k,:)' = w' + A1*v(k-1,:)' +...+ Ap*v(k-p,:)' + eta(k,:)',
%
%  where A=[A1 ... Ap] is the coefficient matrix, and w is a vector of
%  intercept terms that is included to allow for a nonzero mean of the
%  process. The vectors eta(k,:) are independent Gaussian noise
%  vectors with mean zero and covariance matrix C. The noise has a generalized Gaussian
%  distribution with shape parameter beta (def=2) and scale parameter alpha (def=1).
%  If beta=2 (default), the noise is normally distributed.
%  If beta < 2, the noise is super-gaussian (beta=1 is the laplacian).
%  If beta > 2, the noise is super-gaussian
%  (beta=Inf is the uniform distribution).
%
%  The matrix A may also be a cell array equal to n (the length of the
%  process). In this case, A{t} is interpreted as the time-varying VAR
%  coefficient matrix (same dimension as above) at time t.
%
%  The p vectors of initial values for the simulation are taken to
%  be equal to the mean value of the process. (The process mean is
%  calculated from the parameters A and w.) To avoid spin-up effects,
%  the first 10^3 time steps are discarded. Alternatively,
%  tvarsim(w,A,C,n,ndisc) discards the first ndisc time steps.
%
%  tvarsim(w,A,C,[n, ntr]) generates ntr realizations (trials) of
%  length n of the AR(p) process, which are output as the matrices
%  v(:,:,itr) with itr=1,...,ntr.

% -------------------------------------------------------------------------
% IMPORTANT NOTE:
% This function is a modified version of the function arsim.m included in
% ARfit. For the latest release of ARFit please visit ARfit's official site:
% http://www.gps.caltech.edu/~tapio/arfit/
%
% The modifications were performed by:
%
% Tim Mullen
% tim@sccn.ucsd.edu
% http://antillipsi.net
%
% Germán Gómez-Herrero
% german.gomezherrero@tut.fi
% http://www.cs.tut.fi/~gomezher/index.htm
%
%
% This modified version of arsim.m allows the user to specify the
% distribution of the residuals via the parameter beta. It requires
% function ggrnd.m (modification by Gomez-Herrero)
% Additionally, this version allows the user to input a cell array of MVAR
% coefficients (one for each time-point) and thereby generate a time-varying
% VAR[p] process (modification by Tim Mullen)
% -------------------------------------------------------------------------

%  Last modification on 12-May-2011 by Tim Mullen
%  tim@sccn.ucsd.edu
%  Previous modification on 19-March-2008 by German Gomez-Herrero,
%  german.gomezherrero@tut.fi

%  Modified 13-Oct-00
%  Author: Tapio Schneider
%          tapio@gps.caltech.edu

if nargin < 8
    distribution = 'gengauss';
end
if nargin < 7
    % gen. gaussian scale
    alpha = 1;
end
if nargin < 6
    % gen. gaussian shape
    beta = 2;
end
if nargin < 5
    % Discard the first ndisc time steps; if ndisc is not given as input
    % argument, use default
    ndisc = 10^3;
end

if ~iscell(A)
    A = {A};
    istv = false;     % VAR not time-varying
else
    istv = true;      % VAR time-varying
end

m       = size(C,1);                  % dimension of state vectors
p       = size(A{1},2)/m;             % order of process
n       = n_ntr(1);                   % number of time steps
if size(n_ntr) == 1
    ntr   = 1;
else
    ntr   = n_ntr(2);
end

if (p ~= round(p))
    error('Bad arguments.');
end

if (length(w) ~= m | min(size(w)) ~= 1)
    error('Dimensions of arguments are mutually incompatible.')
end
w       = w(:)';                      % force w to be row vector


% Check whether specified model is stable
for k=1:length(A)
    A1 	  = [A{k}; eye((p-1)*m) zeros((p-1)*m,m)];
    lambda  = eig(A1);
    if any(abs(lambda) > 1)
        warning('The specified AR model is unstable at timeindex %d.',k)
    end
end


% Compute Cholesky factor of covariance matrix C
[R, err]= chol(C);                    % R is upper triangular
if err ~= 0
    error('Covariance matrix not positive definite.')
end

% Get ntr realizations of ndisc+n independent Gaussian
% pseudo-random vectors with covariance matrix C=R'*R
randvec = zeros([ndisc+n,m ntr]);

if iscell(distribution)
    % distribution param is of form
    % {@func, parameter_array), and where @func is of form
    % f(x,param1,param2,...) for uniformly random x
    fun = distribution{1};
    randvec = fun(rand(ndisc+n,m,ntr),distribution{2:end});
elseif ischar(distribtion)
    if strcmpi(distribution,'hsec')
        % draw from hyperbolic secant distribution
        randvec = alpha + beta*log(tan(pi*rand(ndisc+n,m,ntr)/2));  % beta = 2/pi
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % These are the modifications performed by Germán Gómez-Herrero
        %   for itr=1:ntr
        %     randvec(:, :, itr) = randn([ndisc+n,m])*R;
        %   end
        for itr = 1:ntr
            randvec(:,:,itr) = reshape(ggrnd((ndisc+n)*m,beta,alpha),ndisc+n,m)*R;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
else
    error('tvarsim:bad_distribution','Unknown form for noise distribution');
end

% Add intercept vector to random vectors
randvec = randvec + repmat(w, [ndisc+n, 1, ntr]);

% Get transpose of system matrix A (use transpose in simulation because
% we want to obtain the states as row vectors)
for k=1:length(A)
    AT{k}      = A{k}';
end

% Take the p initial values of the simulation to equal the process mean,
% which is calculated from the parameters A and w
if any(w)
    %  Process has nonzero mean    mval = inv(B)*w'    where
    %             B = eye(m) - A1 -... - Ap;
    %  Assemble B
    B 	 = eye(m);
    for j=1:p
        B = B - A{1}(:, (j-1)*m+1:j*m);
    end
    %  Get mean value of process
    mval = w / B';
    
    %  The optimal forecast of the next state given the p previous
    %  states is stored in the vector x. The vector x is initialized
    %  with the process mean.
    x    = repmat(mval, [p, 1]);
else
    %  Process has zero mean
    x    = zeros(p,m);
end

% Initialize state vectors
u      = repmat([x; zeros(ndisc+n,m)], [1, 1, ntr]);

% Simulate ntr realizations of n+ndisc time steps. In order to be
% able to make use of Matlab's vectorization capabilities, the
% cases p=1 and p>1 must be treated separately.
if p==1
    for itr=1:ntr
        for k=2:ndisc+n+1;
            x(1,:) = u(k-1,:,itr)*AT{(k-1)*istv+1};
            u(k,:,itr) = x + randvec(k-1,:,itr);
        end
    end
else
    for itr=1:ntr
        for k=p+1:ndisc+n+p;
            for j=1:p;
                x(j,:) = u(k-j,:,itr)*AT{(k-1)*istv+1}((j-1)*m+1:j*m,:);
            end
            u(k,:,itr) = sum(x)+randvec(k-p,:,itr);
        end
    end
end

% return only the last n simulated state vectors
v = u(ndisc+p+1:ndisc+n+p,:,:);





