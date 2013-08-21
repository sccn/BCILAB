% [xk,niter,residuals,outputData,opts] =NESTA_UP(A,At,b,Lambda,La,muf,opts)
%
% Solves a L1 or TV + quadratic term minimization problem, with continuation:
%
%     min_x lambda || U x ||_1 + 1/2||b - Ax||_2^2 
% 
% Continuation is performed by sequentially applying Nesterov's algorithm
% with a decreasing sequence of values of  mu0 >= mu >= muf
%
% The primal prox-function is also adapted by accounting for a first guess
% xplug that also tends towards x_muf 
%
% The observation matrix A need not be a projector
%
% Inputs:   A and At - measurement matrix and adjoint (either a matrix, in which
%               case At is unused, or function handles), dimensions m x n.
%           b   - Observed data, a m x 1 array
%           Lambda - Lagrange multiplier
%               Common heuristic: Lambda = sigma*sqrt(2*log(n)),
%               where sigma=std(noise).
%           La - Lipschitz constant of the quadratic term; La = ||A||^2
%           muf - The desired value of mu at the last continuation step.
%               A smaller mu leads to higher accuracy.
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
%  Outputs:
%           xk  - estimate of the solution x
%           niter - number of iterations
%           residuals - first column is the residual at every step,
%               second column is the value of f_mu at every step
%            outputData - a matrix, where each row r is the output
%               from opts.outFcn, if supplied.
%           out - structure with the parameters used
%       
% Written by: Jerome Bobin, Caltech
% Email: bobin@acm.caltech.edu
% Created: May 2009
% Modified (version 1.1): Nov 2009, Stephen Becker, Caltech
%
% NESTA Version 1.1
%   See also Core_Nesterov_UP


function [xk,niter,residuals,outputData,opts] = NESTA_UP(A,At,b,Lambda,La,muf,opts)

%---- Set defaults
if nargin < 7, opts = []; end
if isempty(opts) && isnumeric(opts), opts = struct; end
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

if isempty(xplug) || norm(xplug) < 10*eps
    if isa(A,'function_handle')
        x_ref=At(b);
    else
        x_ref = A'*b;
    end
    if isempty(xplug)
        xplug = x_ref;
    end
else
    x_ref = xplug;
end

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
        if isa(U,'function_handle'), UUtz = U(Ut(z)); else UUtz = U*(U'*z); end
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


muL = Lambda/La;
mu0 = max(mu0,muL);
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
    [xk,niter_int,res,out,optsOut] = Core_Nesterov_UP(...
        A,At,b,Lambda,La,mu,opts);
    
    xplug = xk;
    niter = niter_int + niter;
    
    residuals = [residuals; res];
    outputData = [outputData; out];

end
opts=optsOut;


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
                opts=rmfield(opts,names{i});
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

end %-- end of NESTA_UP function
