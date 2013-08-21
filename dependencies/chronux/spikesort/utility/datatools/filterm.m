function [Y, carry] = filterm(B, A, X, carry)
%FILTERM           One-dimensional digital filter (memory efficient).
%   Y = FILTERM(B,A,X) filters the data in vector X with the filter
%   described by vectors A and B to create the filtered data Y.  The
%   filter is a "Direct Form II Transposed" implementation of the standard
%   difference equation: 
%
%   a(1)*y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
%                           - a(2)*y(n-1) - ... - a(na+1)*y(n-na)
%  
%   If a(1) is not equal to 1, FILTER normalizes the filter
%   coefficients by a(1).  
%
%   If X is a matrix, FILTERM operates along the columns.  X is not
%   allowed to be of dimension greater than 2.
%
%   [Y,Zf] = FILTERM(B,A,X,Zi) gives access to initial and final
%   conditions, Zi and Zf, of the delays.  Zi is a vector of length Q if
%   X is a vector or of length Q x N if X is an M x N matrix, where
%   Q = MAX(LENGTH(A),LENGTH(B))-1.  The returned Zf is always of type
%   double, although Zi can be of any numeric data type.
%
%   FILTERM is similar to Matlab's built-in FILTER except that it operates
%   on non-double data, avoiding the memory penalty of converting an
%   entire data array to double precision.  (The only difference in
%   calling the two functions is that FILTERM will not operate on arrays
%   of dimension greater than 2).
%
%   If X is of type DOUBLE, INT8, INT16 or INT32, Y will be of the same
%   type, with values clipped to the range of the original data type.  If
%   X is of an unsigned integer type (UINT8, UINT16, UINT32), Y will be of
%   the corresponding signed type.  For example, if X is of type UINT8, Y
%   will be of type INT8 and values larger than [-128,127] will be
%   clipped.  If the data X uses the full range of an unsigned data type,
%   it should be first converted to a data type with larger dynamic range.
%
%   Note that if A = 1, it is more efficient to use FILTERZM.
%
%   See also FILTER, FILTERZM.

% Based on from TMW's FILTER.
%  The following commented code tests FILTERM against FILTER:
%   (also, monitor memory usage to see memory advantage)
% X = int8(randn(2^22,2)*32);   [B,A] = butter(6,0.1);
% disp('Running FILTER ...');   drawnow;  tic;
% Y  = int8(filter(B, A, double(X)));  t(1) = toc;
% disp('Running FILTERM ...');  drawnow;  tic;
% YM = filterm(B, A, X);               t(2) = toc;
% check = ceil(rand(1e6,1)*length(X));
% mse = sqrt(sum(double(Y(check))-double(YM(check)).^2));
% printf('FILTER: %4.3f sec   FILTERM: %4.3f sec     (APPROX) MSE: %5.3f', t(1), t(2), mse);


% # elements to process at once: for some reason ~2^11-2^13 seems to work
% best (probably cpu cache, although the range is similar on computers
% with different cache sizes ...) & exact powers of 2 sometimes cause
% hiccups ... so heuristically, 2000 seems reasonable
chunksize = 2000; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Check Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filter coefficients check
if (A(1) ~= 1),  A = A ./ A(1);  end;   % normalize b/c FILTER does

% input data check
if (ndims(X) > 2), error('FILTERM is undefined for arrays of dimension greater than 2.');  end;
if (isvectord(X) == 2), T = 1;  X = X';  else  T = 0;  end;  % col vect for now
[M,N] = size(X);

% initial conditions check
Q = max(length(A), length(B)) - 1;
if (nargin > 3)
    if (isvectord(carry)), carry = carry(:);  end;  % force column vector
    [Qc,Nc] = size(carry);
    if (Q~=Qc || N~=Nc)
        error('Initial conditions array is not of the correct size.');
    end
    carry = double(carry);
else  carry = zeros(Q,N);
end

% data type check
switch (class(X)),
	case 'double',	          convert = @double;
	case {'uint8','int8'},    convert = @int8;
	case {'uint16','int16'},  convert = @int16;
	case {'uint32','int32'},  convert = @int32;
	otherwise, error('X must be a numeric data type.');
end


%%%%%%%%%%%%%%%%%%%%%%%%%% Do the Filtering %%%%%%%%%%%%%%%%%%%%%%%%%%
if (isa(X, 'double') || M < chunksize),  % no point chunking if double
    [Y,carry] = filter(B, A, double(X), carry);
    Y = feval(convert, Y);
else
    % Perform the filtering in chunks so that only a portion of the data
    % needs to be expanded to double at any given time.
    Y = repmat(feval(convert,0), size(X));
	for start = 1:chunksize:M
		finish = min(M, start+chunksize-1);
		[chunk,carry] = filter(B,A,double(X(start:finish,:)),carry);
		Y(start:finish,:) = feval(convert, chunk);
	end
end


%%%%%%%%%%%%%%%%%%%%%%%%%% Clean Up Outputs %%%%%%%%%%%%%%%%%%%%%%%%%%
if (T),  Y = Y';  end;   % if input was row vect, keep same orientation
if (nargout < 2), clear carry; end;