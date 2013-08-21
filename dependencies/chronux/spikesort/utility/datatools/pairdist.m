function dists = pairdist(X,varargin)
%PAIRDIST          Computes a pairwise Euclidean distance matrix.
%   DISTS = PAIRDIST(X) returns a symmetric positive definite matrix DISTS
%   such that the DISTS(i,j) is the Euclidean distance between the vectors
%   X(i,:) and X(j,:).  X must be a real-valued matrix.
%
%   DISTS = PAIRDIST(X,Y), where X is an M x D matrix and Y is N x D,
%   returns an M x N matrix DISTS such that DISTS(i,j) is the Euclidean
%   distance between X(i,:) and Y(j,:).  X and Y must be real-valued.
%
%   The algorithm used by default utilizes the expansion
%                    (x-y)'(x-y) == x'x + y'y - 2x'y,
%   because the right hand side uses fewer floating point operations.
%   However, round-off error with this technique may produce small but
%   non-zero (even negative) values when x == y.
%
%   DISTS = PAIRDIST(X,Y,'safe') indicates that it is possible for
%   the ith row of X to be identical to the jth row of Y.  This uses a
%   slower algorithm, but guarantees that DISTS(i,j) == 0 in these cases.
%   For the special case of PAIRDIST(X), it is not necessary to specify
%   the 'safe' flag if the only potential zero distances are those along
%   the diagonal of DISTS.
%
%   DISTS = PAIRDIST(..., 'nosqrt') returns DISTS as above, except that
%   the squared distances are returned.  The algorithm to compute the
%   distances finds the squared distances in an intermediate step, so this
%   calculation is faster than returning the Euclidean distance proper.
%
%   DISTS = PAIRDIST(..., 'reuse') attempts to reuse memory from an
%   earlier call to PAIRDIST.  When DISTS is a large matrix, the resulting
%   time savings can be significant.  However, this option can result in
%   unexpected Matlab behavior for the returned variable DISTS.  In
%   particular, clearing DISTS in the calling workspace will not release
%   the memory associated with that variable; to release the memory, you
%   must run the command 'clear functions'.  Further, DISTS may be
%   unexpectedly modified by later calls to PAIRDIST.  For example,
%
%       A = PAIRDIST(X1,Y1);         % Look at the value of A(1,2)
%       B = PAIRDIST(X2,Y2,'reuse'); % The value of A(1,2) will change
%                                    % here even though A was not assigned
%
%   If this behavior is undesirable, make a copy of the returned distance
%   matrix before calling PAIRDIST again.

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parse Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = [];   takeSqrt = 1;   reuseMem = 0;  useSafe = 0;
while (length(varargin) > 0)
    tail = varargin{end};  varargin = varargin(1:end-1);  % chomp last arg
    if (ischar(tail) && strcmpi(tail, 'nosqrt'))
        takeSqrt = 0;
    elseif (ischar(tail) && strcmpi(tail, 'reuse')),
        reuseMem = 1;
	elseif (ischar(tail) && strcmpi(tail, 'safe')),
		useSafe = 1;
    elseif (isa(tail, 'double'))
        Y = tail;   break;    % this should be the leftmost argument
    else
        error('Unknown argument.');
    end
end
if (length(varargin) > 0),  error('Too many arguments.');  end;
if (isempty(Y)),  Y = X;  end;   % self-distance computation

[M,D1] = size(X);
[N,D2] = size(Y);
if (D1~=D2),  error('X and Y must have the same number of columns.');  end;
if (~isreal(X) || ~isreal(Y) || ~isa(X,'double') || ~isa(Y,'double'))
    error('Input matrices must be real-valued matrices of type double.');
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%% Do the Computation %%%%%%%%%%%%%%%%%%%%%%%%%%
dists = CORE_pairdist(X,Y,takeSqrt,reuseMem,useSafe);