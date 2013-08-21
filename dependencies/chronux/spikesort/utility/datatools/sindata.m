function [Y,t] = sindata(dur, Fs, W, A, P, N)
%SINDATA           Generate noisy sinuisoidal data.
%   [Y,t] = SINDATA(DUR,Fs) returns vectors t and Y such that plot(t,Y)
%   draws the sinusoid sin(2*pi*t) on the interval 0..DUR with sampling
%   rate Fs.
%
%   [Y,t] = SINDATA(DUR,Fs,W), for W scalar, instead returns Y as the
%   sinusoid sin(2*pi*W*t).  If W is a vector, t is unchanged, but Y will
%   be a matrix with j'th row given by the sinusoid sin(2*pi*W(j)*t).
%
%   [Y,t] = SINDATA(DUR,Fs,W,A,P,N) further specifies amplitude, phase,
%   and noise parameters.  Here, when arguments W, A, P or N all have the
%   the same number of elements M, Y will then have M rows with j'th row
%   of form 
%                 A(j) * sin(2*pi*W(j)*t + P(j)) + N(j)*WGN
%   where WGN is white gaussian noise.  When any of W, A, P or N is of
%   length M > 1, the other arguments can be the empty matrix, scalar or a
%   vectors of length M.  If the empty matrix is given for an argument,
%   the default value is assumed (W = 1, A = 1, P = 0, N = 0).  If any
%   argument is scalar when another argument is of length M, the scalar
%   argument is treated as a repeated value.  E.g., the following are
%   equivalent:
%                  SINDATA(10,100,[1 2], 3, [pi/2 0], [])
%                  SINDATA(10,100,[1 2], [3 3], [pi/2 0], [0 0])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Defaults %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
freq = 1;      amp = 1;      phase = 0;      noiseamp = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%% Argument Checking %%%%%%%%%%%%%%%%%%%%%%%%%
% first, fill in missing parameters
if ((nargin < 6) || isempty(N)),  N = noiseamp; end;
if ((nargin < 5) || isempty(P)),  P = phase;    end;
if ((nargin < 4) || isempty(A)),  A = amp;      end;
if ((nargin < 3) || isempty(W)),  W = freq;     end;

% next, confirm sizes are consistent
W = col(W);  A = col(A);  P = col(P);  N = col(N);
sz = [numel(W)  numel(A)  numel(P)  numel(N)];
M = unique(sz);  if (any(M==1)), M(M==1) = [];  end;  % (ignore scalars for now)
if ((length(M) > 1) || (any(sz==0)))  % any empties at this pt were caused by col
	error('Sinusoidal parameters must all be vectors with the same number of elements.'); 
end
if (isempty(M)),  M = 1;  end;

% finally, expand scalars
expand = ones(M,1);
if (M > 1),
	if (numel(W)==1),  W = W(expand);  end;
	if (numel(A)==1),  A = A(expand);  end;
	if (numel(P)==1),  P = P(expand);  end;
	if (numel(N)==1),  N = N(expand);  end;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Make the data %%%%%%%%%%%%%%%%%%%%%%%%%%
t = (0:1/Fs:dur);    Nt = length(t);
Y = diag(A) * sin(2*pi*diag(W)*repmat(t,M,1) + repmat(P,1,Nt)) + diag(N)*randn(M,Nt);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Helper Functions %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = col(v)  % force col vector; or empty if not a scalar or 1-Dvector.
s = isvectord(v);
if (s == 2), v = v';
elseif ((s == 0) && (numel(v) > 1)), v = [];
end
