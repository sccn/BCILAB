function n1inf = norm_l1inf(x, g_d,g_t, weights2,weights1)
%NORM_L1INF Norm L1inf
%
%   norm_l1inf(x, g_d,g_t, weights2,weights1) gives back the norm L1inf of x
%
%
%   The input arguments are:
%
%   - x : the vector which we want the norm
%
%   - g_d, g_t: the group vectors. g_d contain the indices of the
%       element to be group and g_t the size of different group g_b is
%       obtained by calling the function back_perm(g_d).
%       
%       Example: x=[x1 x2 x3 x4 x5 x6] 
%                and Group 1: [x1 x2 x4 x5] group 2: [x3 x6]
%           
%               => g_d=[1 2 4 5 3 6] and g_t=[4 2]
%               Or this is also possible
%               => g_d=[4 5 3 6 1 2] and g_t=[2 4]          
%
%   - weights2: weights for a weighted L12-norm works on the norm L2 
%       (default = 1)
%
%   - weights1: weights for a weighted L12-norm works on the norm L1
%       (default = 1)
%
%   Url: http://unlocbox.sourceforge.net/doc/misc/norm_l1inf.php

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
% Date: October 2011
%

% Reshape x in a row
x=reshape(x,size(x,1)*size(x,2),1);


% Optional input arguments
if nargin<2, g_d=1:length(x); end
if nargin<3, g_t=ones(length(x),1); end
if nargin<5, weights1=ones(length(g_t),1); end
if nargin<4, weights2=ones(length(x),1); end

l=length(g_t);

% Compute the norm
if max(g_t)==min(g_t),
    X=x;
    X=X(g_d);
    X=reshape(X,length(x)/l,l)';
    W2=reshape(weights2(g_d),length(x)/l,l)';
    normX2=max(abs(W2.*X),[],2);
    n1inf=sum(weights1.*normX2);
else % group of different size
    
    n1inf=0;
    indice=0;
    X=x;
    X=X(g_d);
    W2=weights2;
    W2=W2(g_d);
    for i=1:l
    
        n1inf=n1inf+weights1(i)*norm(W2(indice+1:indice+g_t(i)).*X(indice+1:indice+g_t(i)),Inf);
        
        indice=indice+g_t(i);
    end

end


end
