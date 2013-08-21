function [z,y,rho,history] = consensus_admm(z0, F, param)
% [z,y,history] = consensus_admm(x0, F, param)
% Consensus form of ADMM [1] with optional linear operator for each term (called SDMM in [2]).
%
% This algorithm solves the following optimization problem:
%
%       minimize sum  f_i(L_i*x)
%           s.t. L_i*x_i - z = 0
%
%       where f_i are closed, proper, convex functions
%             L_i are linear operators such that sum L_i'*L_i is invertible
%
% In:
%   z0 : initial value for the solution; column vector
%
%   F  : the terms to optimize over; given as array of structs with field names:
%         .eval(x)          : application of the function f to x
%
%         .prox(x,gamma,x0) : prox-operator for f at x; solution to the optimization problem
%                             argmin_{z} 0.5*||x - z||_2^2 + gamma*f(z)
%                             where x0 is an optional initial guess
%
%         .L                : linear operator to apply prior to f() (default: [])
%
%         .y0               : initial value for the unscaled dual variable from previous solution
%                             (note: this is implicitly transformed by -L)
%
%   param : structure of optional parameters with the following fields:
%            .maxit   : maximum # of iterations; default=1000
%
%            .rho     : initial coupling parameter (scales the l2 term in the prox operators); default=1
%
%            .abs_tol : absolute tolerance for termination; default=1e-4
%
%            .rel_tol : relative tolerance for termination; default=1e-3
%
%            .rho_update : update rho parameter; default: true
%
%            .rho_cutoff : update rho when primal/dual residual ratio exceeds this value; default=10
%
%            .rho_incr : increase rho by this factor upon update; default=2
%
%            .rho_decr : decrease rho by this factor upon update; default=2
%
%            .verbose : verbosity level; default=1
%
% Out:
%   z     : final solution to the optimization problem
%
%   y     : cell array of final values for the unscaled dual variables, one per function
%
%   history : history of objective-function values
%
% Notes:
%   The interface is compatible with the prox-operators included with UnLocBox or the SPAMS toolbox.
%   The linear operators L lead to some differences between this formulation and consensus ADMM:
%   * F(i).x is F(i).L*z and F(i).u is -(F(i).L/rho)*yi
%   * the consensus variable is a more sophisticated average (using L' and 1/Q)
%
% References:
%  [1] S. Boyd, N. Parika, E. Chu, B. Peleato, J. Eckstein, "Distributed Optimization and Statistical Learning
%      via the Alternating Direction Method of Multipliers", Foundations and Trends in Machine Learning 3(1), 2011
%
%   [2] P. Combettes and J. Pesquet. "A douglas-rachford splitting approach to
%       nonsmooth convex variational signal recovery." Selected Topics in Signal
%       Processing, IEEE Journal of, 1(4):564-574, 2007.
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2013-02-09

% Copyright (C) Christian Kothe, SCCN, 2013, christian@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License as published by the Free Software Foundation; either version 2 of the
% License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
% even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License along with this program; if not,
% write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
% USA

% consensus_admm_version<1.0> -- for the cache

z = z0(:);
m = length(F);

% handle optional arguments
if nargin<2
    param=struct(); end
if ~isfield(param,'maxit')
    param.maxit = 1000; end
if ~isfield(param,'rho')
    param.rho = 1; end
if ~isfield(param,'abs_tol')
    param.abs_tol = 10e-4; end
if ~isfield(param,'rel_tol')
    param.rel_tol = 10e-3; end
if ~isfield(param,'verbose')
    param.verbose = 1; end
if ~isfield(param,'rho_update')
    param.rho_update = true; end
if ~isfield(param,'rho_cutoff')
    param.rho_cutoff = 10; end
if ~isfield(param,'rho_incr')
    param.rho_incr = 1.5; end
if ~isfield(param,'rho_decr')
    param.rho_decr = 1.5; end

% handle function parameters
for i=1:m
    if ~isfield(F(i),'L') || isempty(F(i).L)
        F(i).L = speye(numel(z)); end
    F(i).x = F(i).L*z;
    if ~isfield(F(i),'y0') || isempty(F(i).y0)
        F(i).y0 = zeros(size(F(i).x)); end
    F(i).u = F(i).y0/param.rho;
end

% max dimensionality of involved primal and dual variables
np = max([arrayfun(@(f)length(f.x),F),length(z)]);
nd = max(arrayfun(@(f)length(f.u),F));

% initialize and cache the inverse scaling matrix Q (typically diagonal)
persistent Q_cache;
qid = ['x' hlp_cryptohash({F.L})]; % note: use qid = matrix_hash(vertcat(F.L)); if you're using this function standalone
try
    Q_inv = Q_cache.(qid);
catch %#ok<CTCH>
    fprintf('\nInverting and caching operator matrix... ');
    Q = sparse(0);
    for i=1:m
        Q = Q + F(i).L' * F(i).L; end
    Q_inv = inv(Q);
    if nnz(Q_inv) > numel(Q_inv)*0.25
        Q_inv = full(Q_inv); end
    Q_cache.(qid) = Q_inv;
    fprintf('done.\n');
end


if param.verbose
    fprintf('\n%3s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', 'r norm', 'eps pri', 's norm', 'eps dual', 'objective', 'rho'); end
for k=1:param.maxit
    for i=1:m
        % primal variable update
        Lz = F(i).L*z;
        F(i).x = F(i).prox(Lz + F(i).u,1/param.rho,F(i).x);
        % dual variable update
        F(i).u = F(i).u + Lz - F(i).x;        
    end
    
    % consensus variable update
    zold = z;
    z = 0;
    for i=1:m
        z = z + F(i).L'*(F(i).x-F(i).u); end
    z = Q_inv*z;

    % diagnostics, reporting, termination checks
    history.rho(k)     = param.rho;
    history.objval(k)  = sum(arrayfun(@(f)f.eval(z),F));
    history.r_norm(k)  = sqrt(sum(arrayfun(@(f)sum((f.x-f.L*z).^2),F)));
    history.s_norm(k)  = norm(-param.rho*(z - zold),'fro');        
    history.eps_pri(k) = sqrt(np)*param.abs_tol + param.rel_tol*max([arrayfun(@(f)norm(f.x),F),norm(z)]);
    history.eps_dual(k)= sqrt(nd)*param.abs_tol + param.rel_tol*param.rho*sqrt(sum(arrayfun(@(f)sum(f.u(:).^2),F)));
    if param.verbose
        fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\t%10.2f\n', k, ...
            history.r_norm(k), history.eps_pri(k), history.s_norm(k), history.eps_dual(k), history.objval(k), history.rho(k));
    end
    if history.r_norm(k) < history.eps_pri(k) && history.s_norm(k) < history.eps_dual(k)
        break; end
    
    % optional rho update
    if param.rho_update
        if history.r_norm(k) > param.rho_cutoff * history.s_norm(k)
            param.rho = param.rho * param.rho_incr;
            for i=1:length(F)
                F(i).u = F(i).u / param.rho_incr; end
        elseif history.s_norm(k) > param.rho_cutoff * history.r_norm(k)
            param.rho = param.rho / param.rho_incr;
            for i=1:length(F)
                F(i).u = F(i).u * param.rho_incr; end
        end
    end
end

% calculate final unscaled dual variables
y = arrayfun(@(f)f.u*param.rho,F,'UniformOutput',false);
rho = param.rho;


function hash = matrix_hash(M)
% calculate a hash of a given matrix
hasher = java.security.MessageDigest.getInstance('MD5');
hasher.update(uint8(full(M(:))));
hash = dec2hex(typecast(hasher.digest,'uint8'),2);
hash = ['X' hash(:)'];
