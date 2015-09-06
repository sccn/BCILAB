% dal - dual augmented Lagrangian method for sparse learaning/reconstruction
%
% Overview:
%  Solves the following optimization problem
%   xx = argmin f(x) + lambda*c(x)
%  where f is a user specified (convex, smooth) loss function and c
%  is a measure of sparsity (currently L1 or grouped L1)
%
% Syntax:
%  [ww, uu, status] = dal(prob, ww0, uu0, A, B, lambda, <opt>)
%
% Inputs:
%  prob   : structure that contains the following fields:
%   .obj      : DAL objective function
%   .floss    : structure with three fields (p: primal loss, d: dual loss, args: arguments to the loss functions)
%   .fspec    : function handle to the regularizer spectrum function
%               (absolute values for L1, vector of norms for grouped L1, etc.)
%   .dnorm    : function handle to the conjugate of the regularizer function
%               (max(abs(x)) for L1, max(norms) for grouped L1, etc.)
%   .softth   : soft threshold function
%   .mm       : number of samples (scalar)
%   .nn       : number of unknown variables (scalar)
%   .ll       : lower constraint for the Lagrangian multipliers ([mm,1])
%   .uu       : upper constraint for the Lagrangian multipliers ([mm,1])
%   .Ac       : inequality constraint Ac*aa<=bc for the LMs ([pp,mm])
%   .bc       :                                             ([pp,1])
%   .info     : auxiliary variables for the objective function
%   .stopcond : function handle for the stopping condition
%   .hessMult : function handle to the Hessian product function (H*x)
%   .softth : function handle to the "soft threshold" function
%  ww0    : initial solution ([nn,1)
%  uu0    : initial unregularized component ([nu,1])
%  A          : struct with fields times, Ttimes, & slice.
%   .times    : function handle to the function A*x.
%   .Ttimes   : function handle to the function A'*y.
%   .slice    : function handle to the function A(:,I).
%  B      : design matrix for the unregularized component ([mm,nu])
%  lambda : regularization constant (scalar)
%  <opt>  : list of 'fieldname1', value1, 'filedname2', value2, ...
%   aa        : initial Lagrangian multiplier [mm,1] (default zero(mm,1))
%   tol       : tolerance (default 1e-3)
%   maxiter   : maximum number of outer iterations (default 100)
%   eta       : initial barrier parameter (default 1)
%   eps       : initial internal tolerance parameter (default 1e-4)
%   eta_multp : multiplying factor for eta (default 2)
%   eps_multp : multiplying factor for eps (default 0.5)
%   solver    : internal solver. Can be either:
%               'nt'   : Newton method with cholesky factorization
%               'ntsv' : Newton method memory saving (slightly slower)
%               'cg'   : Newton method with PCG (default)
%               'qn'   : Quasi-Newton method
%   display   : display level (0: none, 1: only the last, 2: every
%               outer iteration, (default) 3: every inner iteration)
%   iter      : output the value of ww at each iteration
%               (boolean, default 0)
% Outputs:
%  ww     : the final solution
%  uu     : the final unregularized component
%  status : various status values
%
% Reference:
% "Super-Linear Convergence of Dual Augmented Lagrangian Algorithm
% for Sparse Learning."
% Ryota Tomioka, Taiji Suzuki, and Masashi Sugiyama. JMLR, 2011.
% "Dual Augmented Lagrangian Method for Efficient Sparse Reconstruction"
% Ryota Tomioka and Masashi Sugiyama
% http://arxiv.org/abs/0904.0584
%
% Copyright(c) 2009-2011 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt



function [xx, uu, status]=dal(prob, ww0, uu0, A, B, lambda, varargin)

opt=propertylist2struct(varargin{:});
opt=set_defaults(opt, 'aa', [],...
    'tol', 1e-3, ...    % outer minimization tolerance
    'iter', 0, ...
    'miniter', 5, ...   % minimum number of iterations (before checking stopping conditions...)
    'maxiter', 100,...
    'eta', [],...
    'eps', 1, ...       % inner minimization tolerance
    'eps_multp', 0.99,...
    'eta_multp', 2, ...
    'solver', 'cg', ...
    'boostb', 1,...
    'display',2);


prob=set_defaults(prob, 'll', -inf*ones(prob.mm,1), ...
    'uu', inf*ones(prob.mm,1), ...
    'Ac', [], ...
    'bc', [], ...
    'info', [], ...
    'finddir', []);

% generate L-BFGS options struct...
optlbfgs = lbfgs_defaults(struct('display',opt.display-1));

% sanity check
if opt.miniter > opt.maxiter
    opt.miniter = opt.maxiter; end

% set eta
if isempty(opt.eta)
    opt.eta = 0.01/lambda; end
if ~isempty(uu0) && length(opt.eta)==1
    opt.eta = opt.eta*[1 1]; end
if ~isempty(uu0) && length(opt.eta_multp)<2
    opt.eta_multp = opt.eta_multp*[1 1]; end

% display diagnostics
if opt.display>0
    if ~isempty(uu0)
        nuu = length(uu0);
        vstr=sprintf('%d+%d',prob.nn,nuu);
    else
        vstr=sprintf('%d',prob.nn);
    end
    
    lstr=func2str(prob.floss.p); lstr=lstr(6:end-1);
    fprintf(['DAL ver1.05\n#samples=%d #variables=%s lambda=%g ' ...
        'loss=%s solver=%s\n'],prob.mm, vstr, lambda, lstr, ...
        opt.solver);
end

% init statistics
res    = nan*ones(1,opt.maxiter);
fval   = nan*ones(1,opt.maxiter);
etaout = nan*ones(length(opt.eta),opt.maxiter);
time   = nan*ones(1,opt.maxiter);
xi     = nan*ones(1,opt.maxiter);
num_pcg= nan*ones(1,opt.maxiter);
if opt.iter
    nwu = length(ww0(:))+length(uu0(:));
    xx  = [[ww0(:); uu0(:)], ones(nwu,opt.maxiter-1)*nan];
end

% init main parameters
time0=cputime;
ww   = ww0;
uu   = uu0;
gtmp = zeros(size(ww));

% check feasibility, and correct if necessary
if isempty(opt.aa)
    [ff,gg] = evalloss(prob, ww, uu, A, B);
    aa = -gg;
    if any(aa==prob.ll) || any(aa==prob.uu)
        fprintf('invalid initial solution; using ww=zeros(n,1).\n');
        w0 = zeros(prob.nn,1);
        % Set the initial Lagrangian multiplier as the gradient
        [ff,gg]=evalloss(prob, w0, uu, A, B);
        aa = -gg;
    end
else
    aa = opt.aa;
end

% init a few more parameters
dval = inf;
eta  = opt.eta;
epsl = opt.eps;
info = prob.info;
info.solver=opt.solver;
info.ATaa=[];

% pre-compute a few properties for blks, but only if we don't have it already
persistent known_blks;
for k=1:length(known_blks)
    if isequalwithequalnans(info.blks,known_blks(k).key)
        [info.blkival,info.blkvec,info.blkgrp,info.blksizes] = known_blks(k).val{:};
        break;
    end
end
if ~isfield(info,'blkival')
    % group blocks according to their size
    [blocksizes,dummy,blockindices] = unique_bc(info.blks); %#ok<ASGLU>
    % compute index vectors for each block
    I = cell(length(info.blks),1);
    ix0 = 0;
    for kk=1:length(info.blks)
        I{kk} = ix0+(1:info.blks(kk));
        ix0 = I{kk}(end);
    end
    info.blkival = I;
    % for each group of blocks
    for k=1:length(blocksizes)
        % calc index vector into ss
        info.blkvec{k} = find(blockindices==k);
        % calc index matrix into vv
        info.blkgrp{k} = cell2mat(I(info.blkvec{k}));
    end
    info.blksizes = blocksizes;
    cacheitem = struct('key',{info.blks},'val',{{info.blkival,info.blkvec,info.blkgrp,info.blksizes}});;
    if isempty(known_blks)
        known_blks = cacheitem;
    else
        known_blks(end+1) = cacheitem;
    end
end

% include info as parameter to some functions
if nargin(prob.fspec)==2
    prob.fspec = @(ww)prob.fspec(ww,info); end
if nargin(prob.dnorm)==2
    prob.dnorm = @(ww)prob.dnorm(ww,info); end

spec = prob.fspec(ww);

% for each iteration..
for ii=1:opt.maxiter-1
    ww_old = ww;
    uu_old = uu;
    
    % keep stats
    etaout(:,ii)=eta';
    time(ii)=cputime-time0;
    
    % Evaluate objective and Check stopping condition
    if ii > opt.miniter
        fval(ii) = evalprim(prob, ww, uu, A, B, lambda);
        switch(prob.stopcond)
            case 'pdg'
                % duality gap
                dval    = min(dval,evaldual(prob, aa, A, B, lambda));
                res(ii) = (fval(ii)-(-dval))/fval(ii);
                ret     = (res(ii)<opt.tol);
            case 'fval'
                % relative functional value
                res(ii) = fval(ii)-opt.tol;
                ret     = res(ii)<=0;
        end
    else
        ret = 0;
    end
    if ret~=0
        break;
    end
    
    % Display
    if opt.display>1 || opt.display>0 && ret~=0
        nnz = full(sum(spec>0));
        if length(eta)==1
            fprintf('[[%d]] fval=%g #(xx~=0)=%d res=%g eta=%g \n', ii, ...
                fval(ii), nnz, res(ii), eta);
        else
            fprintf('[[%d]] fval=%g #(xx~=0)=%d res=%g eta=[%g %g] \n', ii, ...
                fval(ii), nnz, res(ii), eta(1), eta(2));
        end
    end    
    
    % Save the original dual variable for daltv2d
    info.aa0 = aa;
    
    % Solve minimization with respect to aa
    if length(opt.eps)>1
        epsl=opt.eps(ii); end
    switch(opt.solver)
        case 'qn'
            optlbfgs.epsginfo = epsl;
            [aa,stat] = lbfgs(@objfunc3,aa,prob.ll,prob.uu,prob.Ac,prob.bc,info,optlbfgs);
            stat.num_pcg = stat.kk;
        case 'cg'
            fun = @(aa,info) prob.obj(aa, info, prob,ww,uu,A,B,lambda,eta);
            funh = @(xx,Hinfo)prob.hessMult(xx,A,eta,Hinfo);
            fh = {fun, funh};
            [aa,dfval,dgg,stat] = newton(fh, aa, prob.ll, prob.uu, prob.Ac, ...
                prob.bc, epsl, prob.finddir, info, opt.display>2); %#ok<ASGLU>
        case {'nt','ntsv'}
            [aa,dfval,dgg,stat] = newton(@objfunc4, aa, prob.ll, prob.uu, prob.Ac, ...
                prob.bc, epsl, prob.finddir, info, opt.display>2); %#ok<ASGLU>
        otherwise
            error('Unknown method [%s]',opt.solver);
    end
    info=stat.info;
    xi(ii)=info.ginfo;
    num_pcg(ii)=stat.num_pcg;
    
    % Update primal variable
    if isfield(prob,'Aeq')
        I1 = 1:mm-prob.meq;
        I2 = mm-prob.meq+1:mm;
        gtmp(:) = A.Ttimes(aa(I1))+prob.Aeq'*aa(I2);
        [ww,spec] = prob.softth(ww+eta(1)*gtmp,eta(1)*lambda,info);
    else
        % this is the default code path
        ww  =info.wnew;
        spec=info.spec;
    end
    
    if ~isempty(uu)
        if isfield(prob,'Aeq')
            uu  = uu+eta(2)*(B'*aa(1:end-prob.meq));
        else
            uu  = uu+eta(2)*(B'*aa);
        end
    end
    
    % Boosting the bias term
    if length(eta)>1
        viol = [norm(ww-ww_old)/eta(1), norm(uu-uu_old)/eta(2)];
        if opt.boostb && ii>1 && viol(2)>viol_old*0.5 && viol(2)>min(0.001,opt.tol)
            eta(2)=eta(2)*20.^(stat.ret==0);
        end
        %if (opt.display>1 || opt.display>0 && ret~=0)
        % fprintf('violation = [%g %g]\n', viol(1), viol(2));
        %end
        viol_old = viol(2);
    end
    
    % Update barrier parameter eta and tolerance parameter epsl
    eta     = eta.*opt.eta_multp.^(stat.ret==0);
    epsl    = epsl*opt.eps_multp^(stat.ret==0);
    if opt.iter
        xx(:,ii+1)=[ww(:);uu(:)];
    end
end

% clean up and return
res(ii+1:end)=[];
fval(ii+1:end)=[];
time(ii+1:end)=[];
etaout(:,ii+1:end)=[];
xi(ii+1:end)=[];
num_pcg(ii+1:end)=[];

if opt.iter
    xx(:,ii+1:end)=[];
else
    xx = ww;
end

status=struct('aa', aa,...
    'niter',length(res),...
    'eta', etaout,...
    'xi', xi,...
    'time', time,...
    'res', res,...
    'opt', opt, ...
    'info', info,...
    'fval', fval,...
    'num_pcg',num_pcg);


    % objective functions, faster than using a lambda func
    function [a,b,c] = objfunc3(aa,info)
        [a,b,c] = prob.obj(aa, info, prob,ww,uu,A,B,lambda,eta);
    end

    function [a,b,c,d] = objfunc4(aa,info)
        [a,b,c,d] = prob.obj(aa, info, prob,ww,uu,A,B,lambda,eta);
    end

end


% evaluate the loss function
function [fval,gg]=evalloss(prob, ww, uu, A, B)
fnc=prob.floss;
if ~isempty(uu)
    zz=A.times(ww)+B*uu;
else
    zz=A.times(ww);
end
[fval, gg] = fnc.p(zz, fnc.args{:});
end

% Evaluate primal objective
function fval = evalprim(prob, ww, uu, A, B, lambda)
spec=prob.fspec(ww);
fval = evalloss(prob,ww,uu,A,B)+lambda*sum(spec);
if isfield(prob,'Aeq')
    fval = fval+norm(prob.Aeq*ww-prob.ceq)^2/tol;
end
end

% Evaluate dual objective
function dval = evaldual(prob, aa, A, B, lambda)
mm=length(aa);
fnc=prob.floss;
if ~isempty(B)
    if isfield(prob,'Aeq')
        aa1=aa(I1);
        aa(I1)=aa1-B*((B'*B)\(B'*aa1));
    else
        aa=aa-B*((B'*B)\(B'*aa));
    end
end
if isfield(prob,'Aeq')
    I1=1:mm-prob.meq;
    I2=mm-prob.meq+1:mm;
    vv = A.Ttimes(aa(I1))+prob.Aeq'*aa(I2);
else
    vv = A.Ttimes(aa);
end
[dnm,ishard] = prob.dnorm(vv);
if ishard && dnm>0
    aa  = min(1, lambda/dnm)*aa;
    dnm = 0;
end
dval = fnc.d(aa, fnc.args{:})+dnm;
end

