function sol = prox_l1inf(x, gamma , param)
%PROX_L1inf Proximal operator with L1inf norm
%   Usage:  sol=prox_linf(x, gamma, param)
%
%   Input parameters:
%         x     : Input signal.
%         gamma : Regularization parameter.
%         param : Structure of parameters.
%   Output parameters:
%         sol   : Solution.
%
%   prox_L1inf(x, gamma, param) solves:
%
%      sol = argmin_{z} 0.5*||x - z||_2^2 + gamma  |A X|1inf
%
%
%   param is a Matlab structure containing the following fields:
%
%    param.weights1 : weights for a weighted L1inf-norm works on the norm L1 (default = 1)
%
%    param.weights2 : weights for a weighted L1inf-norm works on the sup nom (default = 1)
%
%    param.g_d, param.g_t are the group vectors. 
%     (default param.g_d=1:length(x);, param.g_t=ones(size(x));, param.g_b=back_perm(g_d); )
%
%     param.g_d contains the indices of the elements to be grouped and param.g_t the size of the different groups.
%
%     Warning: param.g_d and param.g_t have to be row vector!     
%     
%     Example: suppose x=[x1 x2 x3 x4 x5 x6] and Group 1: [x1 x2 x4 x5] group 2: [x3 x6]
%              
%     In matlab: 
%
%           param.g_d=[1 2 4 5 3 6]; param.g_t=[4 2];
%
%     Or this is also possible:
%
%           param.g_d=[4 5 3 6 1 2]; param.g_t=[2 4]; 
%
%    param.multi_group*: in order to group component in a not disjoint manner, it
%     is possible to use the multi_group option. param.multi_group 
%     is now set automatically by the function. 
%
%     Overlaping group:
%     In order to make overlapping group just give a vector of g_d, g_b
%     and g_t. Example:
%       
%           param.g_d=[g_d1; g_d2; ...; g_dn];
%           param.g_t=[g_t1; g_t2; ...; g_tn];
%
%     Warning! There must be no overlap in g_d1, g_d2,... g_dn
%
%    param.verbose : 0 no log, 1 a summary at convergence, 2 print main
%     steps (default: 1)
%
%   See also:  prox_l1 prox_l12 proj_b1 prox_sumg
%
%   Demos: demo_compress_sensing4
%
%   References:
%     F. Bach, R. Jenatton, J. Mairal, and G. Obozinski. Optimization with
%     sparsity-inducing penalties. arXiv preprint arXiv:1108.0775, 2011.
%     
%
%   Url: http://unlocbox.sourceforge.net/doc/prox/prox_l1inf.php

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


% Author: Nathanael Perraudin
% Date: November 2012
%

% Reshape x if not a row vector
t=size(x);
x=x(:);

% Optional input arguments
if nargin<3, param=struct; end

if ~isfield(param, 'g_d'),    param.g_d = 1:length(x); end
if ~isfield(param, 'g_t'),    param.g_t = ones(length(x),size(param.g_d,2)); end
if ~isfield(param, 'weights1'), param.weights1=ones(length(param.g_t),1) ; end
if ~isfield(param, 'weights2'), param.weights2=ones(length(x),1) ; end
if ~isfield(param, 'verbose'), param.verbose=1 ; end
if ~isfield(param, 'multi_group'), param.multi_group=0 ; end

% Test of gamma
gamma=test_gamma(gamma);
% Test of the weights
param.weights1=test_weights(param.weights1);
param.weights2=test_weights(param.weights2);

test_multigroup(x,param.g_d,param.g_t);

% test if there are more than one group
if size(param.g_d,1)>1
    param.multi_group=1;
end



if param.multi_group==0;
    % Number of group
    l=length(param.g_t);
   
    % Test if all the group have the same size
    if max(param.g_t)==min(param.g_t),
        
        % Useful functions
        softB = @(z, T) sign(z).*max(abs(z)-T.*abs(z), 0);
        
        % reshape x in a useful manner
        X=x;
        X=X(param.g_d);
        X=reshape(X,length(x)/l,l)';
        W2=reshape(param.weights2(param.g_d),length(x)/l,l)';
        
        % soft thresholding
        temp=W2.*X;
        S=gamma./max(abs(temp),[],2).*param.weights1;
        
        sol=softB( temp,repmat(S,1,max(param.g_t)));
        
        %reconstruct the solution
        sol=reshape(sol',length(x),1);
        sol(param.g_d)=sol;
        
        
        

        
    else % group of different size
        sol=zeros(size(x));
        indice=0;
        % Useful functions
        softB = @(z, T) sign(z).*max(abs(z)*(1-T), 0);
        
        xp=x;
        xp=xp(param.g_d);
        W2=param.weights2;
        W2=W2(param.g_d);
        for i=1:l
           temp=xp(indice+1:indice+param.g_t(i));
           w2=W2(indice+1:indice+param.g_t(i));
           s=softB( temp, ...
               gamma/norm(w2.*temp,Inf)*param.weights1(i));
           sol(indice+1:indice+param.g_t(i))=s;
           indice=indice+param.g_t(i);
           
        end
        %reconstruct the solution
        sol(param.g_d)=sol;
    end
    
    
    
    norm_L1inf=norm_l1inf(x,param.g_d,param.g_t,param.weights2,param.weights1);
    % Log after the calculous of the prox
    if param.verbose >= 1
        fprintf('  prox_L12: ||x||_12 = %e\n', norm_L1inf);
    end

else
    r=size(param.g_t,1);
    
    % Parameter for the prox

    G=[];
    
    for k=1:r
        param3.g_t=param.g_t(k,:);
        param3.g_d=param.g_d(k,:);
        param3.multi_group=0;
        param3.verbose=0;
        g.prox=@(x,T) prox_l1inf(x,T,param3);
        g.eval=@(x) norm_l1inf(x,param.g_d(k,:),param.g_t(k,:));
        G=[G,g];
    end
    
    param4.G=G;
    param4.verbose=param.verbose-1;
    sol=prox_sumg(x,gamma,param4);
    
end

%reshape the solution
sol=reshape(sol,t);

end