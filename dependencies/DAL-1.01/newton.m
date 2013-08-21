% newton - a simple implementation of the Newton method
%
% Syntax:
%  [xx,fval,gg,status]=newton(fun, xx, ll, uu, Ac, bc, tol, finddir, info, verbose, varargin);
%
% Copyright(c) 2009 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt
%
function [xx,fval,gg,status]=newton(fun, xx, ll, uu, Ac, bc, tol, ...
    finddir, info, verbose, varargin)
% In:
%   fun : objective function (may also be a cell array, if a Hessian is available)
%   xx  : initial solution (the initial gradient)
%   ll  : lower bound for the Lagrangian multipliers
%   uu  : upper bound for the Lagrangian multipliers
%   Ac  : linear constraints for the Lagrangian multipliers
%   bc  : part of the linear constraints
%   tol : solution error tolerance
%   finddir : this is the direction proposal function as part of the newton method (which accounts
%             for the gradient and the Hessian, as far as available...)
%             * can be a cholesky method (apparently this preconditions the matrix in some way, e.g. PSD?)
%             * or an be preconditioned conjugate gradient
%             after a direction has been chosen, the line search is used to determine a reasonable
%             step size
%   info : auxilliary variables for the objective function
%   verbose : verbosity
%   ... misc stuff...
%
% Out:
%   xx     : optimal value
%   fval   : function value at optimum
%   gg     : gradient at optimum
%   status : some status info...

if isempty(verbose)
    verbose=0;
end

if iscell(fun)
    fcnHessMult = fun{2};
    fun = fun{1};
else
    fcnHessMult = [];
end

if isempty(finddir)
    if isempty(fcnHessMult)
        % cholesky method doesn't require multiplying by the Hessian... (but still requires a point value)
        finddir = @dir_chol;
    else
        % but pcg is a really complex beast in itself!
        finddir = @dir_pcg;
    end
end


n = size(xx,1);
if verbose
    fprintf('n=%d tol=%g\n', n, tol);
end

% for each iteration...
cc = 0;
step = nan;
while 1
    % evaluate objective function (get function value, gradient, hessian, ...)
    [fval,gg,H,info]=fun(xx, info);
    
    % display stuff...
    if verbose
        fprintf('[%d] fval=%g norm(gg)=%g step=%g\n',cc, fval, norm(gg),step);
    end
    if info.ginfo<tol % || norm(gg)<1e-3
        if verbose
            fprintf('Optimization success! ginfo=%g\n',info.ginfo);
        end
        ret = 0;
        status=archive('ret','cc','info');
        break;
    end
    
    if isstruct(H)
        H.fcnHessMult=fcnHessMult;
    end
    
    % call the appropriate finddir function; this requires quite a complex structure in the Cg case...
    % (looks somewhat hacky...)
    dd = finddir(gg, H);
    
    % determine an appropriate step size
    [step,fval,info] = linesearch(fun, fval, xx, ll, uu, Ac, bc, dd, info, varargin{:});
    
    % this apparently means that we cannot move any further in this direction...
    if step==0
        ret = -2;
        fprintf('[newton] max linesearch reached. ginfo=%g\n',info.ginfo);
        status=archive('ret','cc','info');
        break;
    end

    % update the current solution
    xx = xx + step*dd;
    
    % update number of iterations...
    cc = cc+1;
    if cc>1000
        ret = -1;
        fprintf('[newton] #iterations>1000.\n');
        status=archive('ret','cc','info');
        break;
    end
end


function dd = dir_chol(gg, H)
% this direction proposal function requires an explicit hessian, and one can in principle also
% use conjugate gradient with that! (though it's not clear how fast that would be)
% but that is an EXCELLENT baseline to compare how & why this damn Hessian multiply FAILS!!!!
R = chol(H);
dd = -R\(R'\gg);
%dd = pcg(H, -gg, max(1e-6,tol*0.01));


function dd = dir_pcg(gg, Hinfo)
% this function uses the Hessian only implicitly, in the form of that multiply function...
% so the result should be an implicit product between a vector and a matrix.
% thus, if the vector has the right dimension, then any error in this function must be a well-contained
% bug
S=warning('off','all');
[dd,dum1,dum2] = pcg(Hinfo.fcnHessMult, -gg, 1e-2, length(gg),Hinfo.prec,[],[],Hinfo);
warning(S);

% [dd,dum1,dum2] = pcg(@(xx)fcnHessMult(xx,Hinfo), -gg, max(1e-6,tol*0.01), length(xx),Hinfo.prec);




function [step,fval, info] = linesearch(fun, fval0, xx, ll, uu, Ac, bc, dd, info, varargin)
% the line search looks for the appropriate step size (in a bisection manner, it appears)
Ip=find(dd>0);
In=find(dd<0);
% and INTERESTINGLY, this is where the lower & upper bounds enter the game...
step=min([1.0, 0.999*min((xx(In)-ll(In))./(-dd(In))), 0.999*min((uu(Ip)-xx(Ip))./dd(Ip))]);

xx0    = xx;

cc = 1;
while 1
    % walk a step
    xx = xx0 + step*dd;
    
    % check for constraints
    if ~isempty(Ac)
        bineq = all(Ac*xx<=bc);
    else
        bineq = true;
    end
    if bineq && all(ll<=xx) && all(xx<=uu)
        % check if we managed to minimize...
        [fval,info] = fun(xx, info);
        if fval<fval0
            break;
        end
    else
        % keyboard;
    end
    
    %fprintf('step=%g fval=%g (fval0=%g)\n',step, fval, fval0);
    
    % if we didn't minimize, half the step size...
    step = step/2;
    
    % after 30 iterations, stop the line search (we're stuck in here...)
    cc = cc+1;
    if cc>50
        fval=fval0;
        step = 0;
        break;
    end
end
