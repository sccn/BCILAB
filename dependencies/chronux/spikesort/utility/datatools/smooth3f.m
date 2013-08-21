function smoothed = smooth3f(data, filt, sz, arg)
%SMOOTH3F          Smooth 3D data (fast version).
%   W = SMOOTH3F(V) smoothes input data V with a Gaussian kernel.  The
%   smoothed data is returned in W.
%   
%   W = SMOOTH3F(V, 'filter') Filter can be 'gaussian' or 'box' (default)
%   and determines the convolution kernel.
%   
%   W = SMOOTH3F(V, 'filter', SIZE) sets the size of the convolution
%   kernel (default is [3 3 3]). If SIZE is a scalar, the size is 
%   interpreted as [SIZE SIZE SIZE].  Each element of SIZE is required
%   to be an odd integer.
%   
%   W = SMOOTH3F(V, 'filter', SIZE, ARG) sets an attribute of the
%   convolution kernel. When filter is 'gaussian', ARG is the standard
%   deviation (default is .65).  If filter is 'box', ARG has no effect.
%
%   SMOOTH3F is similar to Matlab's built-in SMOOTH3 but uses a more
%   efficient algorithm.  (The only difference in calling the two
%   functions is that SMOOTH3F requires an odd SIZE argument).
%
%   See also SMOOTH3.

% Modified from TMW's SMOOTH3.
%  The following commented code tests SMOOTH3f against SMOOTH3:
% data = randn([50,50,50]);
% tic; orig = smooth3(data);  t(1) = toc;
% tic; modf = smooth3f(data); t(2) = toc;
% mse = sqrt(mean(abs(orig(:)-modf(:)).^2));   % mean squared error
% printf('SMOOTH3: %4.3f sec   SMOOTH3F: %4.3f sec    MSE: %5.3f', t(1), t(2), mse);

%%%%%%%%%%%%%%%%%%%%%%%%%%% Parse Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%
if (nargin==1)      %smooth3(data)
  arg = .65;    sz = 3;   filt = 'b';
elseif (nargin==2)  %smooth3(data, filter)
  arg = .65;    sz = 3;
elseif (nargin==3)  %smooth3(data, filter, sz)
  arg = .65;
elseif ((nargin>4) || (nargin==0))
  error('Wrong number of input arguments.'); 
end

if (ndims(data)~=3), error('V must be a 3D array.');  end;

if (numel(sz)==1)
	sz = [sz sz sz];
elseif (numel(sz)~=3)
	error('SIZE must be a scalar or a 3 element vector.')
end
sz = sz(:)';  % force column vector

szHalf = (sz-1)/2;
if (~isequal(szHalf,round(szHalf))), error('SIZE must contain only odd integers.');  end;

%%%%%%%%%%%%%%%%%%%%%%%%%% Make the kernel %%%%%%%%%%%%%%%%%%%%%%%%
% Make three kernels so that the full convolution kernel is the
% outer product of the three ...
switch(lower(filt(1))),
	case 'g',   % 3gaussian
		kernel{1} = gausskernel(szHalf(1),arg);
		kernel{2} = gausskernel(szHalf(2),arg);
		kernel{3} = gausskernel(szHalf(3),arg);
	case 'b',   % box
		kernel{1} = ones(sz(1),1)./sz(1);
		kernel{2} = ones(sz(2),1)./sz(2);
		kernel{3} = ones(sz(3),1)./sz(3);
	otherwise,
		error('Unknown filter.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% Do the Smooth %%%%%%%%%%%%%%%%%%%%%%%%%
% Its grossly inefficient to do the full convolution since the kernel
% is separable.  We use CONVNSEP to do three 1-D convolutions.
smoothed = convnsep(kernel{:}, padreplicate(data,(sz-1)/2), 'valid');


%%%%% TAKEN FROM TMW's SMOOTH3 rev 1.7 -- pads an array by replicating values.
function b=padreplicate(a, padSize)
numDims = length(padSize);
idx = cell(numDims,1);
for k = 1:numDims
  M = size(a,k);
  onesVector = ones(1,padSize(k));
  idx{k} = [onesVector 1:M M*onesVector];
end
b = a(idx{:});
