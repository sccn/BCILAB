function [yinds,dists] = mindist(X,Y)
%MINDIST           Finds indices of closest vectors.
%   YINDS = MINDIST(XLIST,YLIST), where XLIST is an M x D matrix and Y is
%   N x D, returns an M x 1 vector YINDS such that YINDS(i) is the row
%   number in YLIST that has the minimum (Euclidean) distance to the i^th
%   row in XLIST.
%
%   [YINDS,DISTS] = MINDIST(XLIST,YLIST) also returns the M x 1 vector
%   DISTS such that DISTS(i) is the (Euclidean) distance between
%   YLIST(YINDS(i),:) and XLIST(i,:).
%
%   [...] = MINDIST(XLIST) is equivalent to MINDIST(XLIST,XLIST).

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parse Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nargin < 2),  Y = X;  end;
[M,D1] = size(X);
[N,D2] = size(Y);
if (D1~=D2),  error('X and Y must have the same number of columns.');  end;
if (~isreal(X) || ~isreal(Y) || ~isa(X,'double') || ~isa(Y,'double'))
    error('Input matrices must be real-valued matrices of type double.');
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%% Do the Computation %%%%%%%%%%%%%%%%%%%%%%%%%%
[yinds,dists] = CORE_mindist(X,Y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Clean Up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nargin < 2),    clear dists;  
else                dists = sqrt(dists);
end;