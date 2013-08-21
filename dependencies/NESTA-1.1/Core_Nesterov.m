function [xk,niter,residuals,outputData,opts] = Core_Nesterov(...
    A,At,b,mu,delta,opts)
% [xk,niter,residuals,outputData,opts] =Core_Nesterov(A,At,b,mu,delta,opts)
%
% Solves a L1 minimization problem under a quadratic constraint using the
% Nesterov algorithm, without continuation:
%
%     min_x || U x ||_1 s.t. ||y - Ax||_2 <= delta
% 
% If continuation is desired, see the function NESTA.m
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
%                   If the SVD of A is U*S*V', then AAtinv = U*(S^{-2})*U'.
%               opts.USV - another experimental option.  This supercedes
%                   the AAtinv option, so it is recommended that you
%                   do not define AAtinv.  This allows the use of a matrix
%                   A which is not a projection, and works for the
%                   noisy ( i.e. delta > 0 ) case.
%                   opts.USV should contain three fields: 
%                   opts.USV.U  is the U from [U,S,V] = svd(A)
%                   likewise, opts.USV.S and opts.USV.V are S and V
%                   from svd(A).  S may be a matrix or a vector.
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
% Modified: May 2009, Jerome Bobin and Stephen Becker, Caltech
% Modified: Nov 2009, Stephen Becker
%
% NESTA Version 1.1
%   See also NESTA

%---- Set defaults
% opts = [];
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
xplug = setOpts('xplug',[]);
normU = setOpts('normU',1);

if delta < 0, error('delta must be greater or equal to zero'); end

if isa(A,'function_handle')
    Atfun = At;
    Afun = A;
else
    Atfun = @(x) A'*x;
    Afun = @(x) A*x;
end
Atb = Atfun(b);

AAtinv = setOpts('AAtinv',[]);
USV = setOpts('USV',[]);
if ~isempty(USV)
    if isstruct(USV)
        Q = USV.U;  % we can't use "U" as the variable name
                    % since "U" already refers to the analysis operator
        S = USV.S;
        if isvector(S), s = S; S = diag(s);
        else s = diag(S); end
        V = USV.V;
    else
        error('opts.USV must be a structure');
    end
    if isempty(AAtinv)
        AAtinv = Q*diag( s.^(-2) )*Q';
    end
end
% --- for A not a projection (experimental)
if ~isempty(AAtinv)
    if isa(AAtinv,'function_handle')
        AAtinv_fun = AAtinv;
    else
        AAtinv_fun = @(x) AAtinv * x;
    end
    
    AtAAtb = Atfun( AAtinv_fun(b) );

else
    % We assume it's a projection
    AtAAtb = Atb;
    AAtinv_fun = @(x) x;
end

if isempty(xplug)
    xplug = AtAAtb; 
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

%---- Computing Atb
Atb = Atfun(b);
Axk = Afun(xk);% only needed if you want to see the residuals
% Axplug = Axk;


%---- TV Minimization
if strcmpi(TypeMin,'TV')
    Lmu = 8*Lmu;
    Dv = spdiags([reshape([-ones(n-1,n); zeros(1,n)],N,1) ...
        reshape([zeros(1,n); ones(n-1,n)],N,1)], [0 1], N, N);
    Dh = spdiags([reshape([-ones(n,n-1) zeros(n,1)],N,1) ...
        reshape([zeros(n,1) ones(n,n-1)],N,1)], [0 n], N, N);
    D = sparse([Dh;Dv]);
end


Lmu1 = 1/Lmu;
% SLmu = sqrt(Lmu);
% SLmu1 = 1/sqrt(Lmu);
lambdaY = 0;
lambdaZ = 0;

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
   
   if strcmpi(TypeMin,'L1')  [df,fx,val,uk] = Perform_L1_Constraint(xk,mu,U,Ut);end
   
   if strcmpi(TypeMin,'TV')  [df,fx] = Perform_TV_Constraint(xk,mu,Dv,Dh,D,U,Ut);end
   
   %---- Primal Problem
   
   %---- Updating yk 
    
    %
    % yk = Argmin_x Lmu/2 ||x - xk||_l2^2 + <df,x-xk> s.t. ||b-Ax||_l2 < delta
    % Let xp be sqrt(Lmu) (x-xk), dfp be df/sqrt(Lmu), bp be sqrt(Lmu)(b- Axk) and deltap be sqrt(Lmu)delta
    % yk =  xk + 1/sqrt(Lmu) Argmin_xp 1/2 || xp ||_2^2 + <dfp,xp> s.t. || bp - Axp ||_2 < deltap
    %
    
    
    cp = xk - 1/Lmu*df;  % this is "q" in eq. (3.7) in the paper
    
    Acp = Afun( cp );
    if ~isempty(AAtinv) && isempty(USV)
        AtAcp = Atfun( AAtinv_fun( Acp ) );
    else
        AtAcp = Atfun( Acp );
    end
    
    residuals(k+1,1) = norm( b-Axk);    % the residual
    residuals(k+1,2) = fx;              % the value of the objective
    %--- if user has supplied a function, apply it to the iterate
    if RECORD_DATA
        outputData(k+1,:) = outFcn(xk);
    end
    
    if delta > 0
        if ~isempty(USV)
            % there are more efficient methods, but we're assuming
            % that A is negligible compared to U and Ut.
            % Here we make the change of variables x <-- x - xk
            %       and                            df <-- df/L
            dfp = -Lmu1*df;  Adfp = -(Axk - Acp);
            bp = b - Axk;
            deltap = delta;
            % Check if we even need to project:
            if norm( Adfp - bp ) < deltap
                lambdaY = 0;  projIter = 0;
                % i.e. projection = dfp;
                yk = xk + dfp;
                Ayk = Axk + Adfp;
            else
                lambdaY_old = lambdaY;
                [projection,projIter,lambdaY] = fastProjection(Q,S,V,dfp,bp,...
                    deltap, .999*lambdaY_old );
                if lambdaY > 0, disp('lambda is positive!'); keyboard; end
                yk = xk + projection;
                Ayk = Afun(yk);
                % DEBUGGING
%                 if projIter == 50
%                     fprintf('\n Maxed out iterations at y\n');
%                     keyboard
%                 end
            end
        else
            lambda = max(0,Lmu*(norm(b-Acp)/delta - 1));gamma = lambda/(lambda + Lmu);
            yk = lambda/Lmu*(1-gamma)*Atb + cp - gamma*AtAcp;
            % for calculating the residual, we'll avoid calling A()
            % by storing A(yk) here (using A'*A = I):
            Ayk = lambda/Lmu*(1-gamma)*b + Acp - gamma*Acp;
        end
    else
        % if delta is 0, the projection is simplified:
        yk = AtAAtb + cp - AtAcp;
        Ayk = b;
    end

    % DEBUGGING
%     if norm( Ayk - b ) > (1.05)*delta
%         fprintf('\nAyk failed projection test\n');
%         keyboard;
%     end
    
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
    
    %
    % zk = Argmin_x Lmu/2 ||b - Ax||_l2^2 + Lmu/2||x - xplug ||_2^2+ <wk,x-xk> 
    %   s.t. ||b-Ax||_l2 < delta
    %
    
    cp = xplug - 1/Lmu*wk;
    
    Acp = Afun( cp );
    if ~isempty( AAtinv ) && isempty(USV)
        AtAcp = Atfun( AAtinv_fun( Acp ) );
    else
        AtAcp = Atfun( Acp );
    end
    
    if delta > 0
        if ~isempty(USV)
            % Make the substitution wk <-- wk/K
                 
%             dfp = (xplug - Lmu1*wk);  % = cp
%             Adfp= (Axplug - Acp);
            dfp = cp; Adfp = Acp; 
            bp = b;
            deltap = delta;            
%             dfp = SLmu*xplug - SLmu1*wk;
%             bp = SLmu*b;
%             deltap = SLmu*delta;

            % See if we even need to project:
            if norm( Adfp - bp ) < deltap
                zk = dfp;
                Azk = Adfp;
            else
                [projection,projIter,lambdaZ] = fastProjection(Q,S,V,dfp,bp,...
                    deltap, .999*lambdaZ );
                if lambdaZ > 0, disp('lambda is positive!'); keyboard; end
                zk = projection;
                %             zk = SLmu1*projection;
                Azk = Afun(zk);
            
                % DEBUGGING:
%                 if projIter == 50
%                     fprintf('\n Maxed out iterations at z\n');
%                     keyboard
%                 end
            end
        else
            lambda = max(0,Lmu*(norm(b-Acp)/delta - 1));gamma = lambda/(lambda + Lmu);
            zk = lambda/Lmu*(1-gamma)*Atb + cp - gamma*AtAcp;
            % for calculating the residual, we'll avoid calling A()
            % by storing A(zk) here (using A'*A = I):
            Azk = lambda/Lmu*(1-gamma)*b + Acp - gamma*Acp;
        end
    else
        % if delta is 0, this is simplified:
        zk = AtAAtb + cp - AtAcp;
        Azk = b;
    end
    
    % DEBUGGING
%     if norm( Ayk - b ) > (1.05)*delta
%         fprintf('\nAzk failed projection test\n');
%         keyboard;
%     end

    %--- Updating xk
    
    xkp = tauk*zk + (1-tauk)*yk;
    xold = xk;
    xk=xkp; 
    Axk = tauk*Azk + (1-tauk)*Ayk;
    
    if ~mod(k,10), Axk = Afun(xk); end   % otherwise slowly lose precision
    % DEBUG
%     if norm(Axk - Afun(xk) ) > 1e-6, disp('error with Axk'); keyboard; end
    
    %--- display progress if desired
    if ~mod(k+1,Verbose )
        printf('Iter: %3d  ~ fmu: %.3e ~ Rel. Variation of fmu: %.2e ~ Residual: %.2e',...
            k+1,fx,qp,residuals(k+1,1) ); 
        %--- if user has supplied a function to calculate the error,
        % apply it to the current iterate and dislay the output:
        if DISPLAY_ERROR, printf(' ~ Error: %.2e',errFcn(xk)); end
        printf('\n');
    end
    if abs(fx)>1e20 || abs(residuals(k+1,1)) >1e20 || isnan(fx)
        error('Nesta: possible divergence or NaN.  Bad estimate of ||A''A||?');
    end

end

niter = k+1; 

%-- truncate output vectors
residuals = residuals(1:niter,:);
if RECORD_DATA,     outputData = outputData(1:niter,:); end



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
    fx = real(u'*D*x - mu/2 * 1/numel(u)*sum(u'*u));
    if isa(Ut,'function_handle')
        df = Ut(D'*u);
    else
        df = U'*(D'*u);
    end
end

