function [z, u, history] = admm_gl(varargin)

% group_lasso  Solve group lasso problem via ADMM
%
% [x, history] = group_lasso(A, y, blks, ...);
% 
% solves the following problem via ADMM:
%
%   minimize 1/2*|| Ax - y ||_2^2 + \lambda sum(norm(x_i))
%
% The input p is a K-element vector giving the block sizes n_i, so that x_i
% is in R^{n_i}.
% 
% The solution is returned in the vector x.
%
% history is a structure that contains the objective value, the primal and 
% dual residual norms, and the tolerances for the primal and dual residual 
% norms at each iteration.
% 
% rho is the augmented Lagrangian parameter. 
%
% alpha is the over-relaxation parameter (typical values for alpha are 
% between 1.0 and 1.8).
%
%
% More information can be found in the papers found in [1-2]
%
% 
% Author: Based on example code by Boyd et al [1-2]. 
%         Tim Mullen, 2012, SCCN/INC/UCSD
%
% References:
% [1] http://www.stanford.edu/~boyd/papers/distr_opt_stat_learning_admm.html
% [2] http://www.stanford.edu/~boyd/papers/admm/group_lasso/group_lasso_example.html
g = arg_define(varargin, ...
                arg_norep({'A','DesignMatrix'},mandatory,[],'The design matrix. This is the data matrix (ie X).'), ...
                arg_norep({'y','TargetVector'},mandatory,[],'The target vector'), ...
                arg_norep({'blks','Blocks'},mandatory,[],'Block (Group) size. Array containing the number of elements in each block (group) of coefficients. ie: blks = [5 5 10] is two blocks of size 5 followed by a block of size 10.'), ...
                arg_nogui({'designMatrixBlockSize','DesignMatrixBlockSize'},[],[],'Design matrix structure. If empty, do not exploit structure for faster computation. If non-empty, the design matrix consists of identical blocks along the main diagonal. ''DesignMatrixBlockSize'' is the [numrow numcol] size of each block.'), ...
                arg({'lambda','ReguParamLambda','RegularizationParam'},[],[0 Inf],'Regularization parameter (lambda)'), ...
                arg({'rho','AugLagrangParamRho','AugmentedLagrangianParam'},1.0,[0 Inf],'Initial value of augmented Lagrangian parameter (rho)'), ...
                arg({'alpha','OverRelaxationParamAlpha','OverRelaxationParam','Alpha'},1.7,[0 Inf],'Over-relaxation parameter (alpha). Typical values are between 1.5 and 1.8'), ...
                arg_nogui({'z_init','InitialState','initAR'},[],[],'Initial state vector','shape','matrix'), ...
                arg_nogui({'u_init','InitialDual'},[],[],'Initial value of u','shape','matrix'), ...
                arg({'verb','Verbosity'},false,[],'Verbose output'), ...
                arg({'max_iter','MaxIterations'},1000,[],'Maximum number of iterations'), ...
                arg({'compute_objval','ComputeObjectiveValue'},false,[],'Compute objective value function. Slower processing. Useful mainly for diagnostics'), ...
                arg({'compute_history','ComputeHistory'},false,[],'Compute history. Produces a history of intermediate values during optimization.'), ...
                arg({'rho_update','RhoUpdate'},false,[],'Update Rho. Whether to update rho dynamically according to 3.4.1 in [1]. Note, this can sometimes cause r_norm, s_norm to "blow up"'), ...
                arg({'rho_cutoff','RhoUpdateThreshold'},10.0,[],'Rho update threshold.'), ...
                arg({'rho_incr','RhoUpdateIncr'},2.0,[],'Rho update increment factor.'), ...
                arg({'rho_decr','RhoUpdateDecr'},2.0,[],'Rho update decrement factor.'), ...
                arg({'lambda_update','LambdaUpdate'},false,[],'Update Lambda. Whether to decrease lambda (regularization) dynamically if convergence is slow.'), ...
                arg({'lambda_update_thresh','LambdaUpdateThreshold'},10^-5,[],'Update lambda if the change in r_norm is less than this','type','denserealdouble'), ...
                arg({'lambda_update_count','LambdaUpdateCount'},10,[],'Number of iterations before updating lambda. Update lambda if r_norm has not significantly changed after this many iterations.'), ...
                arg({'lambda_update_factor','LambdaUpdateFactor'},10,[],'Update factor by which to divide lambda.'), ...
                arg({'abstol','AbsoluteTolerance'},1e-4,[0 Inf],'Primal/Dual absolute tolerance'), ...
                arg({'reltol','RelativeTolerance'},1e-2,[0 Inf],'Primal/Dual relative tolerance'), ...
                arg_subswitch({'x_update','SolverMethod'},'direct',...
                    {'direct' {}, ...
                    'iterative' ...
                        { ...
                            arg({'max_iters','MaxIterations'},[],[],'Number of iterations. If unspecified, a suitable number is automatically determined by lsqr()'), ...
                            arg({'tol','Tolerance'},1e-6,[],'Stopping tolerance. Iteration terminates when residual norm is smaller than this value. Default value is 1e-6. However, this can be relaxed substantially as described in section 4.3 of [1].'), ...
                        }, ...
                    },{'Method used to update parameters x.', sprintf('\n\nIf ''direct'' then a cholesky factorization of A (data/design matrix) is pre-cached and used to update x directly. This is usually a good choice, unless A is very large.\n\nIf ''iterative'', then x is updated using an iterative method (lsqr.m)')}) ...
                );
     
arg_toworkspace(g);

if ~isempty(g.designMatrixBlockSize)
    % extract first block of A
    bA = full(A(1:g.designMatrixBlockSize(1),1:g.designMatrixBlockSize(2))); 
end

if verb
    t_start = tic;
end

%% defaults
g.compute_objval = nargout > 2 || g.compute_objval;


%% Data preprocessing

if isempty(g.designMatrixBlockSize)
    [m, n] = size(A);
    Aty = A'*y;
    if (sum(blks) ~= n)
        error('invalid partition'); 
    end
else
    n = sum(blks); 
    m = n*g.designMatrixBlockSize(1)/g.designMatrixBlockSize(2);
    Aty = reshape(bA'*reshape(y,g.designMatrixBlockSize(1),[]),[],1);
end
    
% cumulative partition
cum_part = cumsum(blks);
% check if blocks are equal size
eqSizeBlocks = all(blks==blks(1));
nblks        = length(blks);

%% select lambda using heuristic
% example based on reference [2]

if isempty(g.lambda)
    if g.verb, fprintf('Using heuristic lambda selection...'); end
        
    % compute a rough estimate of the variance for each block
    % and store the maximum of the variance estimates.
    if eqSizeBlocks
        % fast path
        lambda_max  = max( norms( reshape(Aty,n/nblks,nblks) ) );
    else
        % slow path
        lambda_max  = -Inf;
        start_ind   = 1;
        for i = 1:nblks
            sel        = start_ind:cum_part(i);
            lambda_max = max( lambda_max, norm(Aty(sel)) );
            start_ind  = cum_part(i) + 1;
        end
    end

    % regularization parameter as a fraction of max block variance
    % (this is our a-priori estimate of the noise variance)
    g.lambda = 0.1*lambda_max;
    
    if g.verb, fprintf('...g.lambda set to %0.10g\n',g.lambda); end
end


%% ADMM solver

x = zeros(n,1);
if isempty(u_init)
    u = zeros(n,1);
else
    u = u_init;
end
if isempty(z_init)
    z = zeros(n,1);  
else
    z = z_init;
end

history = struct(...
            'r_norm'    ,   nan(max_iter,1), ...
            's_norm'    ,   nan(max_iter,1), ...
            'eps_pri'   ,   nan(max_iter,1), ...
            'eps_dual'  ,   nan(max_iter,1), ...
            'objval'    ,   nan(max_iter,1), ...
            'lsqr_iters',   nan(max_iter,1)  ...
            );

if strcmp(x_update.arg_selection,'direct')
    if isempty(designMatrixBlockSize)
        % pre-factor    
        [L,U] = factor(A, g.rho, g.designMatrixBlockSize);
    else
        % pre-calculate the update matrix bM for one block
        if( m >= n ) % if skinny
            bM = inv(bA'*bA + rho*eye(g.designMatrixBlockSize(2))); % TODO: make sure that this is correct
        else            
            bM = bA'/(eye(g.designMatrixBlockSize(1)) + 1/rho*(bA*bA'))*bA;
        end
    end
end

if verb
    fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', ...
      'r norm', 'eps pri', 's norm', 'eps dual', 'objective','lsqr iters');
end

lambda_counter = 0;
is_direct = strcmp(x_update.arg_selection,'direct');
is_iterative = strcmp(x_update.arg_selection,'iterative');

for k = 1:max_iter

    if is_direct
        % x-update using direct method
        q = Aty + g.rho*(z - u);    % temporary value
        if( m >= n )    % if skinny
            if isempty(g.designMatrixBlockSize)
                x = U \ (L \ q);
            else
                x = reshape(bM*reshape(q,g.designMatrixBlockSize(2),[]),[],1); %bM * q;  TODO: check this
            end
        else            % if fat
            if isempty(g.designMatrixBlockSize)
                x = q/g.rho - (A'*(U \ ( L \ (A*q) )))/g.rho^2;
            else
                x = q/g.rho - reshape(bM*reshape(q,g.designMatrixBlockSize(2),[]),[],1)/g.rho^2;
            end
        end
    else
        % x-update using iterative method
        % this uses Matlab's lsqr(); uses previous x to warm start
        [x, flag, relres, iters] = lsqr([A; sqrt(g.rho)*speye(n)], ...
            [y; sqrt(g.rho)*(z-u)],x_update.tol,x_update.max_iters,[],[],x);
    end

    % z-update
    zold = z;
    x_hat = g.alpha*x + (1-g.alpha)*zold;
    if eqSizeBlocks
        % all blocks are same length, so reshape to [blocklen x numblks]
        % and regularize blocks all at once for speed
        z = shrinkage(  reshape(x_hat,n/nblks,nblks) ...
                      + reshape(u,n/nblks,nblks),    ...
                        g.lambda/g.rho );
        z = reshape(z,n,1);
    else
        % blocks are different size, so regularize each in a loop
        start_ind   = 1;
        for i = 1:nblks
            % regularize each block
            sel = start_ind:cum_part(i);
            z(sel) = shrinkage(x_hat(sel) + u(sel), g.lambda/g.rho);
            start_ind = cum_part(i) + 1;
        end
    end
    
    % u-update
    u = u + (x_hat - z);
    
    % diagnostics, reporting, termination checks
    if g.compute_history || g.verb || rho_update || lambda_update
        if is_iterative
            history.lsqr_iters(k) = iters;
        end
        history.r_norm(k)  = sqrt(sum((x-z).^2));    
        history.eps_pri(k) = sqrt(n)*g.abstol + g.reltol*max(sqrt(sum(x.^2)), sqrt(sum(z.^2)));
        history.s_norm(k)  = g.rho*sqrt(sum((z - zold).^2));
        history.eps_dual(k)= sqrt(n)*g.abstol + g.reltol*g.rho*sqrt(sum(u.^2));
        if g.compute_objval
            history.objval(k)  = objective(A, y, g.lambda, cum_part, x, z, eqSizeBlocks);
        end        
        if verb
            fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\t%d\n', k, ...
                history.r_norm(k), history.eps_pri(k),  ...
                history.s_norm(k), history.eps_dual(k), ...
                history.objval(k), history.lsqr_iters(k));
        end        
        if (history.r_norm(k) < history.eps_pri(k) && history.s_norm(k) < history.eps_dual(k))
             break;
        end
    else
        if g.rho*sqrt(sum((z - zold).^2)) < sqrt(n)*g.abstol + g.reltol*g.rho*sqrt(sum(u.^2)) 
            if sqrt(sum((x-z).^2)) < sqrt(n)*g.abstol + g.reltol*max(sqrt(sum(x.^2)), sqrt(sum(z.^2)))
                break;
            end
        end
    end

    % update g.rho              
    if rho_update
        refactor = false;
        if history.r_norm(k) > rho_cutoff * history.s_norm(k)
            g.rho = g.rho .* rho_incr;
            u   = u ./ rho_incr;
            refactor = true;
            if verb, fprintf('  incrementing g.rho to %5.5g.\n',g.rho); end
        elseif history.s_norm(k) > rho_cutoff * history.r_norm(k)
            g.rho = g.rho ./ rho_incr;
            u   = u .* rho_incr;
            refactor = true;
            if verb, fprintf('  decrementing g.rho to %5.5g.\n',g.rho); end
        end
        
        % update cholesky factorization using new g.rho
        if refactor && strcmp(x_update.arg_selection,'direct')
            if verb, fprintf('  updating rho...'); end
            if isempty(designMatrixBlockSize)
                % pre-factor    
                [L,U] = factor(A, g.rho, g.designMatrixBlockSize);
            else
                % pre-calculate the update matrix bM for one block
                if( m >= n ) % if skinny
                    bM = inv(bA'*bA + rho*eye(g.designMatrixBlockSize(2))); % TODO: make sure that this is correct
                else
                    bM = bA'/(eye(g.designMatrixBlockSize(1)) + 1/rho*(bA*bA'))*bA;
                end
            end
            if verb, fprintf('  done.\n'); end
        end
    end
    
    if lambda_update && k > 1
        % check if convergence has not improved in a while
        % and, if so, relax shrinkage a bit
        if  abs(history.r_norm(k)-history.r_norm(k-1)) < lambda_update_thresh ...
            && abs(history.s_norm(k)-history.s_norm(k-1)) < lambda_update_thresh
            lambda_counter = lambda_counter + 1;
            
            if lambda_counter > lambda_update_count
                g.lambda = g.lambda/lambda_update_factor;
                lambda_counter = 0;
                if verb, fprintf('  relaxing g.lambda to %0.7g.\n',g.lambda); end
            end
            
        end
    end
    
    if verb && ~mod(k,15)
    	fprintf('\n%3s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n\n',   ...
                'iter', 'r norm', 'eps pri', 's norm', 'eps dual', ...
                'objective', 'lsqr iters');
    end
end

if k>=max_iter && verb
    fprintf('glADMM:Maximum number of iterations (%d) reached\n',max_iter);
end

if verb
    toc(t_start);
end

end

%% objective function
function p = objective(A, y, lambda, cum_part, x, z, eqSizeBlocks)
    if eqSizeBlocks
        p = 1/2*sum((A*x - y).^2) + lambda*sum(sqrt(sum(reshape(z.*conj(z),cum_part(1),[]))));
    else
        obj = 0;    
        start_ind = 1;
        for i = 1:length(cum_part),
            sel = start_ind:cum_part(i);
            obj = obj + norm(z(sel));
            start_ind = cum_part(i) + 1;
        end
        p = ( 1/2*sum((A*x - y).^2) + lambda*obj );
    end
end

%% shrinkage function
function z = shrinkage(x, kappa)
    % OPTIMIZATION NOTE: the call to norms() can be replaced
    % with sqrt(mtimesx(x,'t',x)) or a call to DNorm2.m (FEX)
    % This will produce faster performance, but requires compilation
    % z = bsxfun(@times,pos(1 - kappa./norms(x)),x);
    z = bsxfun(@times,max(0,1 - kappa./sqrt(sum(x.*conj(x)))),x);    
end

%% cholesky factorization
function [L U] = factor(A, rho, blkSz)
    if nargin < 3, blkSz = []; end
    [m, n] = size(A);
    if ( m >= n )    % if skinny
       L = chol( outerprod(A,'AtA',blkSz) + rho*speye(n), 'lower' );
    else            % if fat
       L = chol( speye(m) + 1/rho*(outerprod(A,'AAt',blkSz)), 'lower' );
    end
    
    % force matlab to recognize the upper / lower triangular structure
    L = sparse(L);
    U = sparse(L');
end

function C = outerprod(A,mode,blkSz)
    if ~isempty(blkSz)
        % A consists of identical blocks on main diag
        % use a more efficient method for product
        blkrows = blkSz(1);
        blkcols = blkSz(2);
        numblks = size(A,1)/blkrows;
        
        % extract block
        Q = full(A(1:blkrows,1:blkcols));
        
        % compute matrix product for block
        switch mode
            case 'AtA'
                C = Q'*Q;
            case 'AAt'
                C = Q*Q';
        end
            
        % put blocks on main diagonal of sparse zero matrix
        C = blkdiageye(sparse(C),numblks);
    else       
        % A is some other structure, multiply the usual way
        switch mode
            case 'AtA'
                C = A'*A;
            case 'AAt'
                C = A*A';
        end
    end
end

%% pos helper
function q = pos(x)
    q = max(x,0);
end
