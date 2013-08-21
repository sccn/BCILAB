% [xk,niter,residuals,outputData,opts] =Core_Nesterov_UP(A,At,b,Lambda,La,mu,opts)
%
% Solves a L1 or TV minimization + quadratic termproblem
% Nesterov algorithm, without continuation:
%
%     min_x lambda || U x ||_1 + 1/2 ||b - Ax||_2^2 
% 
% If continuation is desired, see the function NESTA.m
%
% The primal prox-function is also adapted by accounting for a first guess
% xplug that also tends towards x_muf 
%
% The observation matrix A need not be a projector
%
% Inputs:   A and At - measurement matrix and adjoint (either a matrix, in which
%               case At is unused, or function handles)
%           b   - Observed data, a m x 1 array
%           Lambda - Lagrange multiplier
%               Common heuristic: Lambda = sigma*sqrt(2*log(n)),
%               where sigma=std(noise).
%           La - Lipschitz constant of the quadratic term; La =  ||A||^2
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
%           outputData - a matrix, where each row r is the output
%               from opts.outFcn, if supplied.
%           opts - the options that were used
%
% Written by: Jerome Bobin, Caltech
% Email: bobin@acm.caltech.edu
% Created: May 2009
% Modified: Nov 2009 Stephen Becker
%
% NESTA Version 1.1
%   See also NESTA_UP


function [xk,niter,residuals,outputData,opts] = Core_Nesterov_UP(A,At,b,Lambda,La,mu,opts)

%---- Set defaults
%opts = [];
if nargin < 7, opts = []; end
if isempty(opts) && isnumeric(opts), opts = struct; end
fid = setOpts('fid',1);
function printf(varargin), fprintf(fid,varargin{:}); end
maxiter = setOpts('maxiter',10000,0);
TolVar = setOpts('TolVar',1e-5);
TypeMin = setOpts('TypeMin','L1');
Verbose = setOpts('Verbose',true);
errFcn = setOpts('errFcn',[]);
outFcn = setOpts('outFcn',[]);
stopTest = setOpts('stopTest',1,1,2);
U = setOpts('U', @(x) x );
if ~isa(U,'function_handle')
    Ut = setOpts('Ut',[]);
else
    Ut = setOpts('Ut', @(x) x );
end
normU = setOpts('normU',1);
xplug = setOpts('xplug',[]);
if isempty(xplug)
    if isa(A,'function_handle')
        xplug=At(b);
    else
        xplug = A'*b;
    end
end

%---- Initialization
N = length(xplug);
wk = zeros(N,1); 
xk = xplug;

%---- Init Variables
Ak= 0;
Lmu = normU/mu;
yk = xk;
zk = xk;
fmean = realmin/10;
OK = 0;
n = floor(sqrt(N));

%---- TV Minimization
if strcmp(TypeMin,'tv')
    Lmu = 8*Lmu;
    Dv = spdiags([reshape([-ones(n-1,n); zeros(1,n)],N,1) ...
        reshape([zeros(1,n); ones(n-1,n)],N,1)], [0 1], N, N);
    Dh = spdiags([reshape([-ones(n,n-1) zeros(n,1)],N,1) ...
        reshape([zeros(n,1) ones(n,n-1)],N,1)], [0 n], N, N);
    D = sparse([Dh;Dv]);
end

Lmu = Lambda*Lmu + La;

Lmu1 = 1/Lmu;
SLmu = sqrt(Lmu);
SLmu1 = 1/sqrt(Lmu);

%---- setup data storage variables
[DISPLAY_ERROR, RECORD_DATA] = deal(false);
outputData = deal([]);
residuals = zeros(maxiter,2);
if ~isempty(errFcn), DISPLAY_ERROR = true; end
if ~isempty(outFcn) && nargout >= 4
    RECORD_DATA = true;
    outputData = zeros(maxiter, size(outFcn(xplug),2) );
end

for k = 0:maxiter-1,
    
   %---- Dual problem
   
   if strcmp(TypeMin,'L1')  [df,fx,val,uk] = Perform_L1_Constraint(xk,mu,U,Ut);end
   
   if strcmp(TypeMin,'tv')  [df,fx] = Perform_TV_Constraint(xk,mu,Dv,Dh,D,U,Ut);end
   
   if isa(A,'function_handle'),
        res = A(xk) - b;
        Ares = At(res);
   else
       res = A*xk - b;
       Ares = A'*(res);
   end
   
   df = Lambda*df + Ares;
   fx= Lambda*fx + 1/2*norm(res)^2;
   
   %---- Updating yk 
    
    yk = xk - Lmu1*df;
    
    residuals(k+1,1) = norm(res);    % the residual
    residuals(k+1,2) = fx;              % the value of the objective
    %--- if user has supplied a function, apply it to the iterate
    if RECORD_DATA
        outputData(k+1,:) = outFcn(xk);
    end
    
    %--- Stopping criterion
    qp = abs(fx - mean(fmean))/mean(fmean);
    
    switch stopTest
        case 1
            % look at the relative change in function value
            if qp <= TolVar && OK; break;end
            if qp <= TolVar && ~OK; OK=1; end
        case 2
            % look at the l_inf change from previous iterate
            if k >= 1 && norm( xk - xold, 'inf' ) <= TolVar
                break
            end
    end
    fmean = [fx,fmean];
    if (length(fmean) > 10) fmean = fmean(1:10);end
    
    %--- Updating zk
  
    apk =0.5*(k+1);
    Ak = Ak + apk; 
    tauk = 2/(k+3); 
    
    wk =  apk*df + wk;
    
    zk = xplug - Lmu1*wk;

    %--- Updating xk
    
    xkp = tauk*zk + (1-tauk)*yk;
    xk=xkp; 
    
    %--- display progress if desired
    if ~mod(k+1,Verbose )
        printf('Iter: %3d  ~ fmu: %.3e ~ Rel. Variation of fmu: %.2e ~ Residual: %.2e',...
            k+1,fx,qp,residuals(k+1,1) ); 
        %--- if user has supplied a function to calculate the error,
        % apply it to the current iterate and dislay the output:
        if DISPLAY_ERROR, printf(' ~ Error: %.2e',errFcn(xk)); end
        printf('\n');
    end

end

niter = k+1; 

residuals = residuals(1:niter,:);
if RECORD_DATA
    outputData = outputData(1:niter,:);
end
%end

%---- internal routine for setting defaults
    function var = setOpts(field,default,mn,mx)
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

end %% end of main Core_Nesterov routine

%%%%%%%%%%%% PERFORM THE L1 CONSTRAINT %%%%%%%%%%%%%%%%%%

function [df,fx,val,uk] = Perform_L1_Constraint(xk,mu,U,Ut)

    if isa(U,'function_handle')
        uk = U(xk);
    else
        uk = U*xk;
    end

    fx = uk;
    
    uk = uk./max(mu,abs(uk));
    val = real(uk'*fx);
    fx = real(uk'*fx - mu/2*norm(uk)^2);
    
    if isa(Ut,'function_handle')
        df = Ut(uk);
    else
        df = U'*uk;
    end

end

%%%%%%%%%%%% PERFORM THE TV CONSTRAINT %%%%%%%%%%%%%%%%%%

function [df,fx] = Perform_TV_Constraint(xk,mu,Dv,Dh,D,U,Ut)

    if isa(U,'function_handle')
        x = U(xk);
    else
        x = U*xk;
    end
    df = zeros(size(x));
    
    Dhx = Dh*x;
    Dvx = Dv*x;
    
    tvx = sum(sqrt(abs(Dhx).^2+abs(Dvx).^2));
    w = max(mu,sqrt(abs(Dhx).^2 + abs(Dvx).^2));
    uh = Dhx ./ w;
    uv = Dvx ./ w;
    u = [uh;uv];
    fx = u'*D*x - mu/2 * norm(u)^2;
    if isa(Ut,'function_handle')
        df = Ut(D'*u);
    else
        df = U'*(D'*u);
    end

end