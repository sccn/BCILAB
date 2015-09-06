% lbfgs - L-BFGS algorithm
%
% Syntax:
%  [xx, status] = lbfgs(fun, xx, ll, uu, <opt>)
%
% Input:
%  fun     - objective function
%  xx      - Initial point for optimization
%  ll      - lower bound on xx
%  uu      - upper bound on xx
%  Ac      - inequality constraint:
%  bc      -     Ac*xx<=bc
%  opt     - Struct or property/value list of optional properties:
%   .m          - size of limited memory
%   .epsg       - gradient tolerance
%   .maxiter    - maximum number of iterations
%   .display    - display level
%
% Output:
%  xx      - Final point of optimization
%  status  - Various numbers
%
% Copyright(c) 2009 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt

function [xx, status] = lbfgs(fun, xx, ll, uu, Ac, bc, info, opt, varargin)
% init memory (from cache)
persistent buffers;
try
    lm = buffers{opt.m};
catch
    nn = size(xx,1);
    lm = repmat(struct('s',zeros(nn,1),'y',zeros(nn,1),'ys',0,'alpha',0),[1, opt.m]);
    buffers{opt.m} = lm;
end

% perform initial step (gradient)
[fval,gg,info] = fun(xx, info);
dd = -gg;
kk = 1;
stp = 1/norm(dd);

ixend = 1;
bound = 0;

% for each L-BFGS iteration...
while 1
    xxp = xx;
    ggp = gg;
    
    % Perform line search
    [ret, xx,fval,gg,info,stp]=linesearch_backtracking(fun, xx, ll, uu, Ac, bc, fval, gg, dd, stp, info, opt, varargin{:});
    
    % check for stopping conditions, and display output
    if ret < 0
        fprintf('ginfo=%g\n',info.ginfo);
        break;
    end
    if info.ginfo < opt.epsginfo
        if opt.display>1
            fprintf('Optimization success! ginfo=%g\n',info.ginfo); end
        ret = 0;
        break;
    end
    if kk == opt.maxiter
        if opt.display>0
            fprintf('Maximum #iterations=%d reached.\n', kk); end
        ret = -3;
        break;
    end
    
    % L-BFGS update
    if opt.m>0
        lm(ixend).s = xx-xxp;
        lm(ixend).y = gg-ggp;
        ys = lm(ixend).y'*lm(ixend).s; yy = sum(lm(ixend).y.^2);
        lm(ixend).ys  = ys;
    else
        ys = 1; yy = 1;
    end
    
    bound = min(bound+1, opt.m);
    ixend = (opt.m>0)*(mod(ixend, opt.m)+1);
    
    % Initially set the negative gradient as descent direction
    dd = -gg;
    
    jj = ixend;
    for ii=1:bound
        jj = mod(jj + opt.m -2, opt.m)+1;
        lm(jj).alpha = lm(jj).s'*dd/lm(jj).ys;
        dd = dd -lm(jj).alpha*lm(jj).y;
    end
    
    dd = dd *(ys/yy);
    
    for ii=1:bound
        beta = lm(jj).y'*dd/lm(jj).ys;
        dd = dd + (lm(jj).alpha-beta)*lm(jj).s;
        jj = mod(jj,opt.m)+1;
    end
    
    stp = 1.0;
    kk = kk + 1;
end

status=struct('ret', ret,...
    'kk', kk,...
    'fval', fval,...
    'gg', gg,...
    'info', info,...
    'opt', opt);




function [ret, xx, fval, gg, info, step] = linesearch_backtracking(fun, xx, ll, uu, Ac, bc, fval, gg, dd, step, info, opt, varargin)
% initial gradient; sanity check
dginit=gg'*dd;
if dginit>=0
    if opt.display>0
        fprintf('dg=%g is not a descending direction!\n', dginit);
    end
    step = 0;
    ret = -1;
    return;
end

% init params
xx0 = xx;
f0  = fval;
ATaa0 = info.ATaa; % The value of AT(xx0)
info.ATaa  = [];
cc = 0;

% determine initial step size
xx = xx0 + step*dd;
if ~(all(xx>ll) && all(xx<uu))
    % need to account for box constraints
    Ip=find(dd>0);
    In=find(dd<0);
    step = min([step, 0.999*min((xx0(In)-ll(In))./(-dd(In))), 0.999*min((uu(Ip)-xx0(Ip))./dd(Ip))]);
    xx = xx0 + step*dd;
end

% accomodate planar constraints, if any
if ~isempty(Ac) 
    while ~all(Ac*xx <= bc)
        fval = Inf;
        step = step/2;
        xx   = xx0 + step*dd;
        cc = cc+1;
        if cc >= opt.max_linesearch
            if opt.display>0
                fprintf('Maximum linesearch=%d reached\n', cc); end
            xx   = xx0;
            [fval, gg, info]=fun(xx, info);
            step = 0;
            ret = -2;
            return;
        end
    end
end

% now backtrack until we have descent
while 1
    % evaluate objective function
    if all(xx>ll) && all(xx<uu)
        [fval,gg,info] = fun(xx, info);
        if fval <= (f0 + opt.ftol*step*dginit)
            % success
            ret = 0;
            return;
        end
    else
        % can happen...
        fval = Inf;
    end

    % max line search reached?
    cc = cc+1;
    if cc >= opt.max_linesearch
        if opt.display>0
            fprintf('Maximum linesearch=%d reached\n', cc); end
        xx   = xx0;
        [fval, gg, info]=fun(xx, info);
        step = 0;
        ret = -2;
        return;
    end
    
    % revise the step
    step = step/2;
    xx    = xx0 + step*dd;
    try
        info.ATaa = (ATaa0 + info.ATaa)/2; 
    catch
        % this is extremely rare
    end
end
