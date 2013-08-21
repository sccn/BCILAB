function kernel = gausskernel(R,S)
%GAUSSKERNEL       Creates a discretized N-dimensional Gaussian kernel.
%   KERNEL = GAUSSKERNEL(R,S), for scalar R and S, returns a 1-D Gaussian
%   array KERNEL with standard deviation S, discretized on [-R:R].
%   
%   If R is a D-dimensional vector and S is scalar, KERNEL will be a
%   D-dimensional isotropic Gaussian kernel with covariance matrix
%   (S^2)*eye(D), discretized on a lattice with points [-R(k):R(k)] in the
%   kth dimension.
%
%   If R and S are both D-dimensional vectors, KERNEL will be a
%   D-dimensional anisotropic Gaussian kernel on the lattice described
%   above, but with standard deviation S(k) in the kth dimension.
%
%   If R is scalar and S is a D-dimensional vector, R is treated as
%   as R*ones(D,1).
%   
%   KERNEL is always normalized to sum to 1.

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Check Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = numel(R);
D2 = numel(S);
if (((D > 1) && (~isvectord(R))) || ((D2> 1) && (~isvectord(S)))),
	error('Matrix arguments are not supported.');
end
if ((D>1)&&(D2>1)&&(D~=D2)), 
	error('R and S must have same number of elements (unless one is scalar).');
end;

if (D2>D),  D = D2;  R = R * ones(1,D); end;   % force bins/sigmas 
if (D>D2),  S = S * ones(1,D);  end;           %   to be same length,
R = R(:)';   S = S(:)';                        % and force row vectors

S(S==0) = 1e-5;  % std==0 causes problems below, 1e-5 has same effect

%%%%%%%%%%%%%%%%%%%%%%%%%%% Make the Kernel %%%%%%%%%%%%%%%%%%%%%%%%%%
RR = 2*R + 1;
for k = 1:D
	% Make the appropriate 1-D Gaussian
	grid = (-R(k):R(k))';
	gauss = exp(-grid.^2./(2*S(k).^2));  
	gauss = gauss ./ sum(gauss);

	% Then expand it against kernel-so-far ...
	if (k == 1),
		kernel = gauss;
	else	
		Dpast = ones(1,(k-1));
		expand = repmat(reshape(gauss, [Dpast RR(k)]), [RR(1:k-1) 1]);
		kernel = repmat(kernel, [Dpast RR(k)]) .* expand;
	end
end

