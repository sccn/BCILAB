function C = convnsep(varargin)
%CONVNSEP          N-dimensional convolution with separable kernels.
%   C = CONVNSEP(H1,H2,...,A) performs an N-dimensional convolution of A
%   with the separable kernel given by the outer product of the H1,H2,...
%   such that the input vector Hk is convolved along the kth dimension of
%   the real N-D array A.  If the number of kernels Hk is equal to M-1,
%   then the (HM..HN) kernels are taken to be 1.  No convolution occurs
%   for dimensions k with corresponding kernel Hk = 1.
%
%   C = CONVN(H1,H2,...,A,'shape') controls the size of the answer C:
%     'full'   - (default) returns the full N-D convolution
%     'same'   - returns the central part of the convolution that
%                is the same size as A.
%     'valid'  - returns only the part of the result that can be
%                computed without assuming zero-padded arrays.  The
%                size of the result is max(size(A,k)-size(Hk,k)+1,0)
%                in the kth dimension.
%
%   See also CONVN, CONV2.

% Modified from TMW's CONVN for efficiency.
%  The following commented code tests CONVNSEP against CONVN:
% D = 50;  R = 3;  sd = 1;
% data = randn([D,D,D]) + i*randn([D,D,D]);
% gauss1 = gausskernel(R,sd);    gauss3 = gausskernel([R,R,R],sd);
% tic; orig = convn(data,gauss3,'same');  toc
% tic; modf = convnsep(gauss1,gauss1,gauss1,data,'same'); toc
% mse = sqrt(mean(abs(orig(:)-modf(:)).^2))   % mean squared error


%%%%%%%%%%%%%%%%%%%%%%%%%%% Parse Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%
% determine output shape 
if (nargin < 2), error('At least two arguments are required.'); end;
if (ischar(varargin{end})),  
	shape = varargin{end};   varargin = varargin(1:end-1);
	if (length(varargin) == 1), 
		C = varargin{1};   return;     % If no kernels specified, no work to do
	end;
else
	shape = 'full';
end

% get target matrix
A = varargin{end};   varargin = varargin(1:end-1);
if (~isa(A,'double')),   A = double(A);  end;
D = ndims(A);

% get kernels
H = varargin;   N = length(H);
if (N > D),  error('Can not have more kernels than the number of dimensions in A.'); end;
Hisreal = ones(D,1);
for k = 1:D,
	if (k <= N)
		if ((numel(H{k})>1) && ~isvectord(H{k})), error('All kernels Hk must be vectors.');  end;
		if (~isa(H{k},'double')),  H{k} = double(H{k}(:));  end;  % force col/double
		Hisreal(k) = isreal(H{k});
	else
		H{k} = 1;
		Hisreal(k) = 1;
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Do the Conv %%%%%%%%%%%%%%%%%%%%%%%%%%
C = A;
for k = 1:N,
	orient = ones(1,ndims(A));   orient(k) = numel(H{k});
	kernel = reshape(H{k}, orient);

	if (Hisreal(k) && isreal(C))
		C = convnc(kernel,C);  
	elseif (Hisreal(k) && ~isreal(C))
		C = convnc(kernel,real(C)) + j*convnc(kernel,imag(C));
	elseif (~Hisreal(k) && isreal(C))
		C = convnc(real(kernel),C) + j*convnc(imag(kernel),C);  
	else
		Hr = real(kernel);    Hi = imag(kernel);
		Cr = real(C);         Ci = imag(C);
		C = convnc(Hr,Cr) - convnc(Hi,Ci) + j*(convnc(Hi,Cr) + convnc(Hr,Ci)); 
	end
end

%%%%%%%%%%%%%%%%%%%%%%%% Get the right shape %%%%%%%%%%%%%%%%%%%%%%
% nothing more to do for 'full' shape
if (strcmp(shape,'full')),  return;  end;
 
% but for 'same' or 'valid' we need to crop the conv result
subs = cell(1,ndims(C));
if (strcmp(shape,'same'))  
  for k = 1:D
	  subs{k}  = floor(length(H{k})/2) + (1:size(A,k));  % central region
  end
elseif (strcmp(shape,'valid'))
  for k = 1:D
	  validLen = max(size(A,k)-length(H{k})+1,0);
	  subs{k}  = length(H{k})-1 + (1:validLen);
  end
end
C = C(subs{:});
