function [sol,iter,objectiv] = ppxa(x_0,F,param)
%PPXA Parallel Proximal algorithm
%   Usage: sol = ppxa(x_0,F,param);
%          sol = ppxa(x_0,F);
%          [sol,iter,objectiv] = ppxa(...);
%
%   Input parameters:
%         x_0   : Starting point of the algorithm
%         F     : Array of function to minimize
%         param : Optional parameter
%   Output parameters:
%         sol   : Solution
%         iter  : Current iteration
%         objectiv: vector (evaluation of the objectiv function each iteration)
%
%   PPXA, derived from the Douglas-Rachford algorithm, solves
% 
%      sol = argmin sum(W_i*f_i(x))
%
%
%   for x in R^N, where x is the variable and x_0 is the starting point.
%
%    F is an array of structures representing functions containing 
%     operators inside and eventually the norm. If no norm is defined, the 
%     l^1 norm is used. The prox: F(i).prox and norm: F(i).eval are defined 
%     in the same way as in the Forward-backward and Douglas-Rachford algorithms
%
%    param a Matlab structure containing the following fields:
%
%     General parameters:
%
%      param.W : the weight (all equal by default)
%
%      param.tol : is stop criterion for the loop. The algorithm stops if
%
%           (  n(t) - n(t-1) )  / n(t) < tol,
%      
%
%       where  n(t) = sum W_i*f_i(x) is the objective function at iteration t*
%       by default, tol=10e-4.
%
%      param.maxit : is the maximum number of iteration. By default, it is 200.
% 
%      param.verbose : 0 no log, 1 print main steps, 2 print all steps.
%
%      param.lambda : is the weight of the update term.
%
%      param.gamma : convergence parameter (default 1)
%
%      param.abs_tol : If activated, this stopping critterion is the
%       objectiv function smaller than param.tol. By default 0.
%
%   See also:  sdmm, admm, generalized_forward_backward
%
%   Demos:  demo_ppxa
%
%   References:
%     P. Combettes and J. Pesquet. Proximal splitting methods in signal
%     processing. Fixed-Point Algorithms for Inverse Problems in Science and
%     Engineering, pages 185-212, 2011.
%
%   Url: http://unlocbox.sourceforge.net/doc/solver/ppxa.php

% Copyright (C) 2012 LTS2-EPFL, by Nathanael Perraudin, Gilles Puy,
% David Shuman, Pierre Vandergheynst.
% This file is part of UnLocBoX version 1.1.70
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% Author:  Nathanael Perraudin
% Date: Oct 19 2012

% test the evaluate function
[F] = test_eval(F);

% number of function
m = size(F,2);

    % Optional input arguments
    if nargin<3, param=struct; end

    if ~isfield(param, 'gamma'), param.gamma = 1; end
    if ~isfield(param, 'tol'), param.tol=10e-4 ; end
    if ~isfield(param, 'maxit'), param.maxit=200; end
    if ~isfield(param, 'verbose'), param.verbose=1 ; end
    if ~isfield(param, 'lambda'), param.lambda=0.99 ; end
    if ~isfield(param, 'W'), param.W=ones(m,1)/m ; end
    
    W=param.W;
    
    
   % Reshape x if vector line
    if (size(W,2)>size(W,1))
        W=W';
    end
    
    test_gamma(param.gamma);
    test_sum(W);
    
    % Create a tabke of scructure containing data
    
    
    data=[];
    for ii = 1:m
        data(ii).y=x_0;     %#ok<AGROW>
        data(ii).p=zeros(size(x_0));     %#ok<AGROW>
    end
    

    x = w_sum(W,data,'y');

    
    % outerloop
    curr_norm = 0;
    for i = 1:m
       curr_norm = F(i).eval(x)+curr_norm;
    end
    [~,~,prev_norm,iter,objectiv,~] = convergence_test(curr_norm);
    
    while 1
        
        if param.verbose >= 1
            fprintf('Iteration %i:\n', iter);
        end
        
        % parallel proximal
        % compute updated prox
        for ii = 1:m
            data(ii).p = F(ii).prox(data(ii).y, param.gamma);%#ok<AGROW>
           
        end
        pn = w_sum(W,data,'p');



        % update y
        
        for ii = 1:m
            data(ii).y  = data(ii).y + param.lambda*(2*pn-x-data(ii).p); %#ok<AGROW>
        end

        % update x
        x = x + param.lambda * (pn - x);
        
        % update solution & relative norm
        sol = x;
        curr_norm = 0;
        for ii = 1:m
            curr_norm = F(ii).eval(sol)+curr_norm;
        end
                
        [stop,rel_norm,prev_norm,iter,objectiv,crit] = convergence_test(curr_norm,prev_norm,iter,objectiv,param);
        if stop
            break;
        end
        if param.verbose >= 1
            fprintf('Current objectiv function : %e \t relative norm : %e \n', curr_norm, rel_norm);
        end
        
    end
    
    if param.verbose>=1
        % Stopping criterion
        fprintf(' %i iterations\n', iter);
        fprintf(' Stopping criterion: %s \n\n', crit);
    end
    
end



function test_gamma(gamma)
    if gamma <= 0 
        fprintf('Warning : gamma is not > 0\n');
    end
end


function test_sum(W)
    if (sum(W) > 1+eps) || (sum(W) < 1-eps)
        fprintf('Warning : sum W is not equal to 1');
    end
end

function s=w_sum(w,data,l)
    s=zeros(size(data(1)));
    
    if l=='p'
        for ii=1:length(data)
            s=s+w(ii)*data(ii).p;
        end
    elseif l=='y'
        for ii=1:length(data)
            s=s+w(ii)*data(ii).y;
        end
    else
        error('Fatal error! Unknown field in data!')
    end
end