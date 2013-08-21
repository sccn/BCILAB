function [xk,niter,residuals,outputData,opts] =NESTA(A,At,b,muf,delta,opts)
% [xk,niter,residuals,outputData] =NESTA(A,At,b,muf,delta,opts)
%
% Solves a L1 minimization problem under a quadratic constraint using the
% Nesterov algorithm, with continuation:
%
%     min_x || U x ||_1 s.t. ||y - Ax||_2 <= delta
% 
% Continuation is performed by sequentially applying Nesterov's algorithm
% with a decreasing sequence of values of  mu0 >= mu >= muf
%
% The primal prox-function is also adapted by accounting for a first guess
% xplug that also tends towards x_muf 
%
% The observation matrix A is a projector
%
% Inputs:   A and At - measurement matrix and adjoint (either a matrix, in which
%               case At is unused, or function handles).  m x n dimensions.
%           b   - Observed data, a m x 1 array
%           muf - The desired value of mu at the last continuation step.
%               A smaller mu leads to higher accuracy.
%           delta - l2 error bound.  This enforces how close the variable
%               must fit the observations b, i.e. || y - Ax ||_2 <= delta
%               If delta = 0, enforces y = Ax
%               Common heuristic: delta = sqrt(m + 2*sqrt(2*m))*sigma;
%               where sigma=std(noise).
%           opts -
%               This is a structure that contains additional options,
%               some of which are optional.
%               The fieldnames are case insensitive.  Below
%               are the possible fieldnames:
%               
%               opts.xplug - the first guess for the primal prox-function, and
%                 also the initial point for xk.  By default, xplug = At(b)
%               opts.U and opts.Ut - Analysis/Synthesis operators
%                 (either matrices of function handles).
%               opts.normU - if opts.U is provided, this should be norm(U)
%                   otherwise it will have to be calculated (potentially
%                   expensive)
%               opts.MaxIntIter - number of continuation steps.
%                 default is 5
%               opts.maxiter - max number of iterations in an inner loop.
%                 default is 10,000
%               opts.TolVar - tolerance for the stopping criteria
%               opts.stopTest - which stopping criteria to apply
%                   opts.stopTest == 1 : stop when the relative
%                       change in the objective function is less than
%                       TolVar
%                   opts.stopTest == 2 : stop with the l_infinity norm
%                       of difference in the xk variable is less
%                       than TolVar
%               opts.TypeMin - if this is 'L1' (default), then
%                   minimizes a smoothed version of the l_1 norm.
%                   If this is 'tv', then minimizes a smoothed
%                   version of the total-variation norm.
%                   The string is case insensitive.
%               opts.Verbose - if this is 0 or false, then very
%                   little output is displayed.  If this is 1 or true,
%                   then output every iteration is displayed.
%                   If this is a number p greater than 1, then
%                   output is displayed every pth iteration.
%               opts.fid - if this is 1 (default), the display is
%                   the usual Matlab screen.  If this is the file-id
%                   of a file opened with fopen, then the display
%                   will be redirected to this file.
%               opts.errFcn - if this is a function handle,
%                   then the program will evaluate opts.errFcn(xk)
%                   at every iteration and display the result.
%                   ex.  opts.errFcn = @(x) norm( x - x_true )
%               opts.outFcn - if this is a function handle, 
%                   then then program will evaluate opts.outFcn(xk)
%                   at every iteration and save the results in outputData.
%                   If the result is a vector (as opposed to a scalar),
%                   it should be a row vector and not a column vector.
%                   ex. opts.outFcn = @(x) [norm( x - xtrue, 'inf' ),...
%                                           norm( x - xtrue) / norm(xtrue)]
%               opts.AAtinv - this is an experimental new option.  AAtinv
%                   is the inverse of AA^*.  This allows the use of a 
%                   matrix A which is not a projection, but only
%                   for the noiseless (i.e. delta = 0) case.
%               opts.USV - another experimental option.  This supercedes
%                   the AAtinv option, so it is recommended that you
%                   do not define AAtinv.  This allows the use of a matrix
%                   A which is not a projection, and works for the
%                   noisy ( i.e. delta > 0 ) case.
%                   opts.USV should contain three fields: 
%                   opts.USV.U  is the U from [U,S,V] = svd(A)
%                   likewise, opts.USV.S and opts.USV.V are S and V
%                   from svd(A).  S may be a matrix or a vector.
%
%  Outputs:
%           xk  - estimate of the solution x
%           niter - number of iterations
%           residuals - first column is the residual at every step,
%               second column is the value of f_mu at every step
%           outputData - a matrix, where each row r is the output
%               from opts.outFcn, if supplied.
%           opts - the structure containing the options that were used    
%
% Written by: Jerome Bobin, Caltech
% Email: bobin@acm.caltech.edu
% Created: February 2009
% Modified (version 1.0): May 2009, Jerome Bobin and Stephen Becker, Caltech
% Modified (version 1.1): Nov 2009, Stephen Becker, Caltech
%
% NESTA Version 1.1
%   See also Core_Nesterov


if nargin < 6, opts = []; end
if isempty(opts) && isnumeric(opts), opts = struct; end

%---- Set defaults
fid = setOpts('fid',1);
Verbose = setOpts('Verbose',true);
function printf(varargin), fprintf(fid,varargin{:}); end
MaxIntIter = setOpts('MaxIntIter',5,1);
TypeMin = setOpts('TypeMin','L1');
TolVar = setOpts('tolvar',1e-5);
[U,U_userSet] = setOpts('U', @(x) x );
if ~isa(U,'function_handle')
    Ut = setOpts('Ut',[]);
else
    Ut = setOpts('Ut', @(x) x );
end
xplug = setOpts('xplug',[]);
normU = setOpts('normU',[]);  % so we can tell if it's been set

residuals = []; outputData = [];
AAtinv = setOpts('AAtinv',[]);
USV = setOpts('USV',[]);
if ~isempty(USV)
    if isstruct(USV)
        Q = USV.U;  % we can't use "U" as the variable name
                    % since "U" already refers to the analysis operator
        S = USV.S;
        if isvector(S), s = S; %S = diag(s);
        else s = diag(S); end
        %V = USV.V;
    else
        error('opts.USV must be a structure');
    end
end

% -- We can handle non-projections IF a (fast) routine for computing
%    the psuedo-inverse is available.
%    We can handle a nonzero delta, but we need the full SVD
if isempty(AAtinv) && isempty(USV)
    % Check if A is a partial isometry, i.e. if AA' = I
    z = randn(size(b));
    if isa(A,'function_handle'), AAtz = A(At(z));
    else AAtz = A*(A'*z); end
    if norm( AAtz - z )/norm(z) > 1e-8
        error('Measurement matrix A must be a partial isometry: AA''=I');
    end
end

% -- Find a initial guess if not already provided.
%   Use least-squares solution: x_ref = A'*inv(A*A')*b
% If A is a projection, the least squares solution is trivial
if isempty(xplug) || norm(xplug) < 1e-12
    if ~isempty(USV) && isempty(AAtinv)
        AAtinv = Q*diag( s.^(-2) )*Q';
    end
    if ~isempty(AAtinv)
        if delta > 0 && isempty(USV)
            error('delta must be zero for non-projections');
        end
        if isa(AAtinv,'function_handle')
            x_ref = AAtinv(b);
        else
            x_ref = AAtinv * b;
        end
    else
        x_ref = b;
    end
    
    if isa(A,'function_handle')
        x_ref=At(x_ref);
    else
        x_ref = A'*x_ref;
    end

    if isempty(xplug)
        xplug = x_ref;
    end
    % x_ref itself is used to calculate mu_0
    %   in the case that xplug has very small norm
else
    x_ref = xplug;
end

% use x_ref, not xplug, to find mu_0
if isa(U,'function_handle')
    Ux_ref = U(x_ref);
else
    Ux_ref = U*x_ref;
end
switch lower(TypeMin)
    case 'l1'
        mu0 = 0.9*max(abs(Ux_ref));
    case 'tv'
        mu0 = ValMUTv(Ux_ref);
end

% -- If U was set by the user and normU not supplied, then calcuate norm(U)
if U_userSet && isempty(normU)
    % simple case: U*U' = I or U'*U = I, in which case norm(U) = 1
    z = randn(size(xplug));
    if isa(U,'function_handle'), UtUz = Ut(U(z)); else UtUz = U'*(U*z); end
    if norm( UtUz - z )/norm(z) < 1e-8
        normU = 1;
    else
        z = randn(size(Ux_ref));
        if isa(U,'function_handle')
            UUtz = U(Ut(z)); 
        else
            UUtz = U*(U'*z);
        end
        if norm( UUtz - z )/norm(z) < 1e-8
            normU = 1;
        end
    end
    
    if isempty(normU)
        % have to actually calculate the norm
        if isa(U,'function_handle')
            [normU,cnt] = my_normest(U,Ut,length(xplug),1e-3,30);
            if cnt == 30, printf('Warning: norm(U) may be inaccurate\n'); end
        else
            [mU,nU] = size(U);
            if mU < nU, UU = U*U'; else UU = U'*U; end 
            % last resort is to call MATLAB's "norm", which is slow
            if norm( UU - diag(diag(UU)),'fro') < 100*eps
                % this means the matrix is diagonal, so norm is easy:
                normU = sqrt( max(abs(diag(UU))) );
            elseif issparse(UU)
                normU = sqrt( normest(UU) );
            else
                if min(size(U)) > 2000
                    % norm(randn(2000)) takes about 5 seconds on my PC
                    printf('Warning: calculation of norm(U) may be slow\n');
                end
                normU = sqrt( norm(UU) );
            end
        end
    end
    opts.normU = normU;
end
        

niter = 0;
Gamma = (muf/mu0)^(1/MaxIntIter);
mu = mu0;
Gammat= (TolVar/0.1)^(1/MaxIntIter);
TolVar = 0.1;
 
for nl=1:MaxIntIter
    
    mu = mu*Gamma;
    TolVar=TolVar*Gammat;    opts.TolVar = TolVar;
    opts.xplug = xplug;
    if Verbose, printf('\tBeginning %s Minimization; mu = %g\n',opts.TypeMin,mu); end
    [xk,niter_int,res,out,optsOut] = Core_Nesterov(...
        A,At,b,mu,delta,opts);
    
    xplug = xk;
    niter = niter_int + niter;
    
    residuals = [residuals; res];
    outputData = [outputData; out];

end
opts = optsOut;


%---- internal routine for setting defaults
function [var,userSet] = setOpts(field,default,mn,mx)
    var = default;
    % has the option already been set?
    if ~isfield(opts,field) 
        % see if there is a capitalization problem:
        names = fieldnames(opts);
        for i = 1:length(names)
            if strcmpi(names{i},field)
                opts.(field) = opts.(names{i});
                opts = rmfield(opts,names{i});
                break;
            end
        end
    end
    if isfield(opts,field) && ~isempty(opts.(field))
        var = opts.(field);  % override the default
        userSet = true;
    else
        userSet = false;
    end
    % perform error checking, if desired
    if nargin >= 3 && ~isempty(mn)
        if var < mn
            printf('Variable %s is %f, should be at least %f\n',...
                field,var,mn); error('variable out-of-bounds');
        end
    end
    if nargin >= 4 && ~isempty(mx)
        if var > mx
            printf('Variable %s is %f, should be at least %f\n',...
                field,var,mn); error('variable out-of-bounds');
        end
    end
    opts.(field) = var;
end




%---- internal routine for setting mu0 in the tv minimization case
function th=ValMUTv(x)

    N = length(x);n = floor(sqrt(N));
    Dv = spdiags([reshape([-ones(n-1,n); zeros(1,n)],N,1) ...
            reshape([zeros(1,n); ones(n-1,n)],N,1)], [0 1], N, N);
        Dh = spdiags([reshape([-ones(n,n-1) zeros(n,1)],N,1) ...
            reshape([zeros(n,1) ones(n,n-1)],N,1)], [0 n], N, N);
        D = sparse([Dh;Dv]);


    Dhx = Dh*x;
    Dvx = Dv*x;
    
    sk = sqrt(abs(Dhx).^2 + abs(Dvx).^2);
    th = max(sk);

end

end %-- end of NESTA function

%%%%%%%%%%%% POWER METHOD TO ESTIMATE NORM %%%%%%%%%%%%%%%
% Copied from MATLAB's "normest" function, but allows function handles, not just sparse matrices
function [e,cnt] = my_normest(S,St,n,tol, maxiter)
%MY_NORMEST Estimate the matrix 2-norm via power method.
    if nargin < 4, tol = 1.e-6; end
    if nargin < 5, maxiter = 20; end
    if isempty(St)
        St = S;  % we assume the matrix is symmetric;
    end
    x = ones(n,1);
    cnt = 0;
    e = norm(x);
    if e == 0, return, end
    x = x/e;
    e0 = 0;
    while abs(e-e0) > tol*e && cnt < maxiter
       e0 = e;
       Sx = S(x);
       if nnz(Sx) == 0
          Sx = rand(size(Sx));
       end
       e = norm(Sx);
       x = St(Sx);
       x = x/norm(x);
       cnt = cnt+1;
    end
end
