function [AR, lambdamax]=sim_genRndVARcoeffs(varargin)
% Generate random coefficients for a VAR model of given order and with a
% specified degree of interaction sparsity.
% Coefficients are drawn from a normal distribution with specified standard
% deviation (sigma).
% 
% The function will repeatedly generate random models until it obtains a
% stable VAR process (all eigenvalues of the system matrix are less than or
% equal to 1 in magnitude).
%
% 
% See Also: sim_genTVARcoeffs()
%
% References:
% [1] This function is adapted from Stefan Haufe's gen_ar_sech.m
%
% Author: Tim Mullen 2011-2013, SCCN/INC, UCSD.
% Email:  tim@sccn.ucsd.edu

% This function is part of the Source Information Flow Toolbox (SIFT)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

g = arg_define(varargin,...
        arg({'nchs','NumChans'},10,[],'Number of channels (variables)'), ...
        arg({'morder','ModelOrder'},10,[],'Model order'), ...
        arg({'sigma','Sigma'},0.2,[0 Inf],'Scale (std) of random VAR parameters. This determines the variance of the determinstic component of the generated process.'), ...
        arg({'maxinter','MaxNumInter'},0.5,[],'Maximum number of interactions. If >= 1, this is the expected number of interacting variables (chosen randomly). If < 1, this represents the proportion of interacting variable (out of nchs^2-nchs max possible interactions). This can also be a vector specifying the linear (column-major) index of each nonzero (non-diagonal) element of the interaction (connectivity) matrix'), ...
        arg({'probinter','PerLagInteractProb'},1,[0 1],'Percent of nonzero AR coefficients for a given interaction. This applies to auto-interactions as well as cross-interactions. For instance, a value of 1 will generate nonzero gaussian iid AR filter coefficients. A value of 0.5 will generate a mixed sparse AR filter with 50% zero coefficients and remaining nonzero iid gaussian coefficients'), ...
        arg({'maxtrys','MaxAttempts','maxattempts'},Inf,[],'Max number of retries for stable model gen.'), ...
        arg({'maxlambda','MaxLambdaMag'},1,[0 1],'Largest allowable eigenvalue (magnitude) of system matrix. This determines the stability of the process. The smaller this value the more stable the process. All eigenvalues must be <= 1 in magnitude for a stable process.') ... 
        );
        
% RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));

if g.maxinter < 1
    g.maxinter = ceil(g.maxinter*(g.nchs^2-g.nchs));
end
    
% generate canonical interaction indices
if length(g.maxinter) == 1
    inddiag     = linspace(1, g.nchs^2, g.nchs);
    indndiag    = setdiff(1:g.nchs^2, inddiag);
    per         = randperm(length(indndiag));
    indndiag    = indndiag(per);
    indndiag    = indndiag(1:g.maxinter);
    ind         = [inddiag indndiag];
else
    inddiag     = linspace(1, g.nchs^2, g.nchs);
    ind         = unique([inddiag g.maxinter]);
end
ind = sort(ind); 
% generate interaction indices for all model orders
ind = tensorsum(0:g.nchs^2:g.nchs^2*(g.morder-1),ind);

% generate random coefficients for interacting variables
lambdamax = Inf; ntrys = 0;
while lambdamax > g.maxlambda && ntrys < g.maxtrys
    AR      = zeros(g.nchs,g.nchs*g.morder);
    AR(ind) = double(rand(length(ind), 1) <= g.probinter).*randn(length(ind), 1)*g.sigma;
    
%     for k=1:g.morder
%       aloc = zeros(g.nchs);
%       aloc(ind) = double(rand(length(ind), 1) < g.probinter).*randn(length(ind), 1)*g.sigma;
%       AR(:,:,(k-1)*g.nchs+1:k*g.nchs) = aloc;
%     end

    % check model stability
    E         = eye(g.nchs*g.morder);
    AA        = [AR;E(1:end-g.nchs,:)];
    lambda    = eig(AA);
    lambdamax = max(abs(lambda));
    
    ntrys = ntrys + 1;
end
if ntrys>=g.maxtrys
    fprintf('Unable to find a stable model. Maximum number of attempts exceeded\n');
    AR = [];
    return;
end




