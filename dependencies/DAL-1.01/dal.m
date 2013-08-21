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
%   .Aeq    : part of some apparent undocumented equality constraints???
%   .ceq    : part of some apparent undocumented equality constraints??
%   
%  ww0    : initial solution ([nn,1)
%  uu0    : initial unregularized component ([nu,1])
%  A      : function handle to the function A*x.
%  AT     : function handle to the function A'*y
%  B      : design matrix for the unregularized component ([mm,nu])
%  lambda : regularization constant (scalar)
%  <opt>  : list of 'fieldname1', value1, 'filedname2', value2, ...
%   aa        : initial Lagrangian multiplier [mm,1] (default zero(mm,1))
%   tol       : tolerance (default 1e-3)
%   maxiter   : maximum number of outer iterations (default 100)
%   eta       : initial barrier parameter (default 1)
%   eps       : initial internal tolerance parameter (default 1e-4)
%               (can also be a vector with 1 value per iteration)
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
% "Dual Augmented Lagrangian Method for Efficient Sparse Reconstruction"
% Ryota Tomioka and Masashi Sugiyama
% http://arxiv.org/abs/0904.0584
%
% Copyright(c) 2009 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt



function [xx, uu, status] = dal(prob, ww0, uu0, A, AT, B, lambda, varargin)

% get solver options
opt=propertylist2struct(varargin{:});
opt=set_defaults(opt, 'aa', [],...
    'tol', 1e-3, ...
    'iter', 0, ...
    'maxiter', 100,...
    'eta', [],...
    'eps', 1, ...
    'eps_multp', 0.99,...
    'eta_multp', 2, ...
    'solver', 'cg', ...
    'display',2);

% get problem description
prob=set_defaults(prob, ...
    'll', -inf*ones(prob.mm,1), ...
    'uu', inf*ones(prob.mm,1), ...
    'Ac', [], ...
    'bc', [], ...
    'info', [], ...
    'finddir', []);

if isempty(opt.eta)
    opt.eta = 0.01/lambda;
end


if opt.display>0
    if ~isempty(uu0)
        nuu = length(uu0);
        vstr=sprintf('%d+%d',prob.nn,nuu);
    else
        vstr=sprintf('%d',prob.nn);
    end
    
    lstr=func2str(prob.floss.p); lstr=lstr(6:end-1);
    fprintf(['DAL ver1.01\n#samples=%d #variables=%s lambda=%g ' ...
        'loss=%s solver=%s\n'],prob.mm, vstr, lambda, lstr, ...
        opt.solver);
end


if opt.iter
    nwu = length(ww0(:))+length(uu0(:));
    xx  = [[ww0(:); uu0(:)], ones(nwu,opt.maxiter-1)*nan];
end

res    = nan*ones(1,opt.maxiter);
fval   = nan*ones(1,opt.maxiter);
etaout = nan*ones(1,opt.maxiter);
time   = nan*ones(1,opt.maxiter);
xi     = nan*ones(1,opt.maxiter);


time0=cputime;
ww   = ww0;
uu   = uu0;
gtmp = zeros(size(ww));
if isempty(opt.aa)
    %% Set the initial Lagrangian multiplier as the gradient at the initial solution, 
    %% w.r.t. design matrix
    [ff,gg]=evalloss(prob, ww, uu, A, B);
    aa = -gg;
else
    % otherwise take the given initial gradient...
    aa = opt.aa;
end

dval = inf;
eta  = opt.eta;
epsl = opt.eps;
info = prob.info;
info.solver=opt.solver;
% for each outer iteration...
for ii=1:opt.maxiter-1
    etaout(ii)=eta;
    time(ii)=cputime-time0;
    
    %% Evaluate primal objective for current solution
    [fval(ii), spec] = evalprim(prob, ww, uu, A, B, lambda);
    
    % and check stopping condition...
    switch(prob.stopcond)
        case 'pdg'
            % relative duality gap 
            dval    = min(dval,evaldual(prob, aa, AT, B, lambda));
            res(ii) = (fval(ii)-(-dval))/fval(ii);
            ret     = (res(ii)<opt.tol);
        case 'fval'
            % function value vs. tolerance...
            res(ii) = fval(ii)-opt.tol;
            ret     = res(ii)<=0;
    end
    
    %% Display
    if opt.display>1 || opt.display>0 && ret~=0
        nnz = full(sum(spec>0));
        fprintf('[[%d]] fval=%g #(xx~=0)=%d res=%g eta=%g \n', ii, ...
            fval(ii), nnz, res(ii), eta);
    end
    
    if ret~=0
        break;
    end
    
    %% Save the original dual variable for daltv2d
    info.aa0 = aa;
    
    %% Solve minimization with respect to aa; ww & uu are fixed, as is the lambda, eta, and the
    %% auxiliary variables in info & opt...
    %% apparently, the result is directly used as the gradient
    fun  = @(aa,info)prob.obj(aa, info, prob,ww,uu,A,AT,B,lambda,eta);
    if length(opt.eps)>1
        epsl=opt.eps(ii);
    end
    switch(opt.solver)
        case 'cg'
            % optimize the DAL objective function (including its Hessian) over aa
            funh = @(xx,Hinfo)prob.hessMult(xx,A,AT,eta,Hinfo);
            fh = {fun, funh};
            [aa,dfval,dgg,stat] = newton(fh, aa, prob.ll, prob.uu, prob.Ac, ...
                prob.bc, epsl, prob.finddir, info, opt.display>2);
        case 'qn'
            % as above, but using quasi-Newton (which doesn't require an explicit Hessian)
            optlbfgs = struct('epsginfo',epsl,'display',opt.display-1);
            [aa,stat] = lbfgs(fun,aa,prob.ll,prob.uu,prob.Ac,prob.bc,info,optlbfgs);
        case {'nt','ntsv'}
            % apparently obsolete, doesn't use conjugate gradient, apparently
            [aa,dfval,dgg,stat] = newton(fun, aa, prob.ll, prob.uu, prob.Ac, ...
                prob.bc, epsl, prob.finddir, info, opt.display>2);
        case 'fminunc'
            % practically obsolete
            optfm=optimset('LargeScale','on','GradObj','on','Hessian', ...
                'on','TolFun',1e-16,'TolX',0,'MaxIter',1000,'display','iter');
            [aa,fvalin,exitflag]=fminunc(@(xx)objdall1fminunc(xx,prob,ww, ...
                uu,A,B,lambda,eta,epsl), aa, optfm);
            stat.info=info;
            stat.ret=exitflag~=1;
        otherwise
            error('Unknown method [%s]',opt.solver);
    end
    info=stat.info;
    xi(ii)=info.ginfo;
    
    
    %% Update primal variable
    if isfield(prob,'Aeq')
        I1=1:mm-prob.meq;
        I2=mm-prob.meq+1:mm;
        gtmp(:) = AT(aa(I1))+prob.Aeq'*aa(I2);
    else
        % get a temporary gradient which points in the right direction... (as A' * aa)
        gtmp(:) = AT(aa);
    end
    
    % shrink the gradient using the soft-thresholding operation
    ww = prob.softth(ww + eta*gtmp, eta*lambda, info);
    
    if ~isempty(uu)
        % also update the unregularized component (as B' * aa)
        if isfield(prob,'Aeq')
            uu  = uu+eta*(B'*aa(1:end-prob.meq));
        else            
            uu  = uu+eta*(B'*aa);
        end
    end
    
    %% Update barrier parameter eta and tolerance parameter epsl
    eta     = eta*opt.eta_multp^(stat.ret==0);
    epsl    = epsl*opt.eps_multp^(stat.ret==0);
    if opt.iter
        xx(:,ii+1)=[ww(:);uu(:)];
    end
end

res(ii+1:end)=[];
fval(ii+1:end)=[];
time(ii+1:end)=[];
etaout(ii+1:end)=[];
xi(ii+1:end)=[];

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
    'fval', fval);


function [fval,gg]=evalloss(prob, ww, uu, A, B)
%% returns the function value and gradient of the primal loss, at the solution ww/uu, given design
% matrices A and B
fnc=prob.floss;
if ~isempty(uu)
    zz=A(ww)+B*uu;
else
    zz=A(ww);
end

[fval, gg] =fnc.p(zz, fnc.args{:});


function [fval,spec] = evalprim(prob, ww, uu, A, B, lambda)
%% Evaluate primal objective (loss + regularizer)
% obtain the regularizer spectrum function
spec = prob.fspec(ww);
% and calculate the primal loss, and add the regularization term
fval = evalloss(prob,ww,uu,A,B) + lambda*sum(spec);

if isfield(prob,'Aeq')
    fval = fval + norm(prob.Aeq*ww-prob.ceq)^2/tol;
end


function dval = evaldual(prob, aa, AT, B, lambda)
%% evaluate dual objective function (loss + regularizer)
mm=length(aa);

fnc=prob.floss;
if ~isempty(B)
    if isfield(prob,'Aeq')
        aa1=aa(I1)
        aa(I1)=aa1-B*((B'*B)\(B'*aa1));
    else
        aa=aa-B*((B'*B)\(B'*aa));
    end
end

if isfield(prob,'Aeq')
    I1=1:mm-prob.meq;
    I2=mm-prob.meq+1:mm;
    vv = AT(aa(I1))+prob.Aeq'*aa(I2);
else
    vv = AT(aa);
end
% evaluate the dual norm
[dnm,ishard] = prob.dnorm(vv);

% Note: whatever is being done here in the "hard" case...??
if ishard && dnm>0
    aa  = min(1, lambda/dnm)*aa;
    dnm = 0;
end

% also evaluate the dual loss, and add the norm
dval = fnc.d(aa, fnc.args{:}) + dnm;


