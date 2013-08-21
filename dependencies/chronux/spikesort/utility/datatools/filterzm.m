function Y = filterzm(B, X)
%FILTERZM          Zero-phase 1-D digital filter (memory efficient).
%   Y = FILTERZ(B, X) filters the data in vector X with the FIR filter
%   given by the vector B to create the filtered data Y, shifting the
%   filter to give zero-phase delay.  The process is equivalent to
%   convolution with the time-reversed and shifted filter B:
%   For P = ceil((length(B)+1)/2),
%     y(n) = ... + x(n-1)*b(P+1) + x(n)*b(P) + x(n+1)*b(P-1) + ...
%
%   If X is a matrix, FILTERZM operates along the columns.  X is not
%   allowed to be of dimension greater than 2.
%
%   FILTERZM(B, X) is similar to FILTER2(B,X,'same'), except that FILTERZM
%   operates on non-double data, avoiding the memory penalty of converting
%   an entire data array to double precision.  It is analogous in this
%   respect to FILTERM and follows the same datatype conventions as that
%   function.  However, by restricting usage only to FIR filters, FILTERZM
%   provides zero-phase delay filtering and is faster than FILTERM(B,1,X).
%
%   See also FILTERM, FILTER2.

% Based on from TMW's FILTER.
%  The following commented code tests FILTERMZ against FILTER:
%   (also, monitor memory usage to see memory advantage)
% X = int8(rand(2^22,2)*64);   B = fir1(6, 0.1);
% disp('Running FILTER ...');   drawnow;  tic;
% Y  = int8(filter(B, 1, double(X)));  t(1) = toc;
% disp('Running FILTERZM ...');  drawnow;  tic;
% YM = filterzm(B, X);                 t(2) = toc;
% check = ceil(rand(1e6,1)*length(X));
% mse = sqrt(sum(double(Y(check))-double(YM(check)).^2));
% printf('FILTER: %4.3f sec   FILTERZM: %4.3f sec    (APPROX) MSE: %5.3f', t(1), t(2), mse);


% # elements to process at once: for some reason ~2^11-2^13 seems to work
% best (probably cpu cache, although the range is similar on computers
% with different cache sizes ...) & exact powers of 2 sometimes cause
% hiccups ... so heuristically, 2000 seems reasonable
chunksize = 2000; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Check Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input data check
if (ndims(X) > 2), error('FILTERZM is undefined for arrays of dimension greater than 2.');  end;
if (isvectord(X) == 2), T = 1;  X = X';  else  T = 0;  end;  % col vect for now
[M,N] = size(X);

% Massage filter into proper shape
if (~isvectord(B)),  error('The FIR filter B must be 1-dimensional.');  end;
if (isvectord(B) == 2), B = B';  end;  % force column vector
B = flipud(B); 
P = (length(B)-1)/2;  Pt = floor(P);  Pb = ceil(P);  % top/bot filter overhang

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
    Y = conv2(double(X), B, 'same');
    Y = feval(convert, Y);
else
    % Perform the filtering in chunks so that only a portion of the data
    % needs to be expanded to double at any given time.
	Y = repmat(feval(convert,0), size(X));
	X = [repmat(feval(convert,0), [Pt+1, N]);  X;  repmat(feval(convert,0), [Pb+1, N])];
	for start = 1:chunksize:M
		finish = min(M, start+chunksize-1);
		chunk = conv2(double((X(start:(finish+Pt+Pb+2),:))),B,'same');
		Y(start:finish,:) = feval(convert, chunk((Pt+2):(end-Pb-1),:));
	end
end


%%%%%%%%%%%%%%%%%%%%%%%%%% Clean Up Outputs %%%%%%%%%%%%%%%%%%%%%%%%%%
if (T),  Y = Y';  end;   % if input was row vect, keep same orientation
