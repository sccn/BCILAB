function data = get_kernel_data(X,kernel, kernel_params, Xtrain, datatrain);
% nargin >= 4, when computations of the testing data depends on the training data (i.e., with incomplete Cholesky)

switch kernel
	case { 'polynomial','polynomial-mkl','polynomial-bimkl'}
		[n , p ] = size(X);
		q = kernel_params(1);
		data.n = n ;
		data.p = p;
		data.q = q;
		data.qs = (q+1) * ones(1,p);		  % number of kernels per directions
		data.Xs = zeros(n,sum(data.qs));      % storing all dimensions in the same vector
		data.ind_Xs = cell(1,p);
		data.d_Xs = cell(1,p);
		iK = 1;
		for i=1:p
			data.ind_Xs{i} = cell(1,data.qs(i));
			data.d_Xs{i} = zeros(1,data.qs(i));
			jK = iK;
			for j=1:data.qs(i)
				data.ind_Xs{i}{j} = iK;
				data.d_Xs{i}(j) = 1;
				if j==1
					data.Xs(:,iK) = ones(n,1);		% constant term
				else
					data.Xs(:,iK) = X(:,i).^(j-1) * sqrt( 1 / factorial(j-1) * factorial(q) / factorial(q-j+1) );
				end
				iK = iK + 1 ;
			end
			jK = jK + data.qs(i);
		end


	case {'hermite', 'hermite-mkl', 'hermite-bimkl'}
		[n , p ] = size(X);
		q = kernel_params(2);
		alpha = kernel_params(1);
		data.n = n ;
		data.p = p;
		data.q = q;
		data.qs = (q+1) * ones(1,p);		  % number of kernels per directions
		data.Xs = zeros(n,sum(data.qs));      % storing all dimensions in the same vector
		data.ind_Xs = cell(1,p);
		data.d_Xs = cell(1,p);
		iK = 1;
		for i=1:p
			data.ind_Xs{i} = cell(1,data.qs(i));
			data.d_Xs{i} = zeros(1,data.qs(i));
			jK = iK;
			for j=1:data.qs(i)
				data.ind_Xs{i}{j} = iK;
				data.d_Xs{i}(j) = 1;

				if j==1
					data.Xs(:,iK) = ones(n,1);
				else
					data.Xs(:,iK) = hermite_polynomials(j,X(:,i)) * sqrt( (alpha/2)^(j-1) / factorial(j-1) );
				end
				iK = iK + 1 ;
			end
			jK = jK + data.qs(i);
		end




	case {'gauss-hermite', 'gauss-hermite-mkl', 'gauss-hermite-bimkl'}
		[n , p ] = size(X);
		q = kernel_params(3);
		sigma = kernel_params(1);
		b = kernel_params(2);
		a = 1/4/sigma;
		c = sqrt(a*a+2*a*b);
		A = a + b + c;
		data.n = n ;
		data.p = p;
		data.q = q;
		data.qs = (q+1) * ones(1,p);		  % number of kernels per directions
		data.Xs = zeros(n,sum(data.qs));      % storing all dimensions in the same vector
		data.ind_Xs = cell(1,p);
		data.d_Xs = cell(1,p);
		iK = 1;
		for i=1:p
			data.ind_Xs{i} = cell(1,data.qs(i));
			data.d_Xs{i} = zeros(1,data.qs(i));
			jK = iK;
			for j=1:data.qs(i)
				data.ind_Xs{i}{j} = iK;
				data.d_Xs{i}(j) = 1;

				data.Xs(:,iK) =  ( 1 - (b/A)^2).^.25 * hermite_polynomials(j,X(:,i) * sqrt( 2 * c ) ) .* exp( - ( b/A*(a+c) * X(:,i).^ 2 ) ) ...
					* sqrt( (b/A)^(j-1) / ( 2^(j-1) * factorial(j-1) )) ;
				iK = iK + 1 ;
			end
			jK = jK + data.qs(i);
		end

	case {'gauss-hermite-full', 'gauss-hermite-full-mkl', 'gauss-hermite-full-bimkl'}
		% simply the same as 'gauss-hermite', but with the full kernels at
		% the end
		[n , p ] = size(X);
		q = kernel_params(3);
		sigma = kernel_params(1);
		b = kernel_params(2);
		a = 1/4/sigma;
		c = sqrt(a*a+2*a*b);
		A = a + b + c;
		data.n = n ;
		data.p = p;
		data.q = q+1;
		data.qs = (q+2) * ones(1,p);     % number of kernels per directions
		data.Xs = zeros(n,sum(data.qs));      % storing all dimensions in the same vector
		data.ind_Xs = cell(1,p);
		data.d_Xs = cell(1,p);
		iK = 1;
		for i=1:p
			data.ind_Xs{i} = cell(1,data.qs(i));
			data.d_Xs{i} = zeros(1,data.qs(i));
			jK = iK;
			fulldataloc = [];
			for j=1:data.qs(i)-1
				% first kernels as before
				data.ind_Xs{i}{j} = iK;
				data.d_Xs{i}(j) = 1;

				data.Xs(:,iK) =  ( 1 - (b/A)^2).^.25 * hermite_polynomials(j,X(:,i) * sqrt( 2 * c ) ) .* exp( - ( b/A*(a+c) * X(:,i).^ 2 ) ) ...
					* sqrt( (b/A)^(j-1) / ( 2^(j-1) * factorial(j-1) )) ;
				fulldataloc = [ fulldataloc; iK ];

				iK = iK + 1 ;
			end


			if nargin < 4
				% last kernel: first build the kernel matrix, then takes its
				% incomplete Cholesky decomposition

				K = exp( - b * sqdist(X(:,i)',X(:,i)' ) ) - data.Xs(:,fulldataloc)  * data.Xs(:,fulldataloc)';
				[Gf,Pf,mf,residualf] = icd_general(K,1e-3,kernel_params(6));
				If = Pf(1:mf);
				data.Is{i} = If;					% these two for test data

				[temp,Pif] = sort(Pf);
				Gf = Gf(Pif,1:mf);
				data.d_Xs{i}(end) = mf;
				data.ind_Xs{i}{end} = iK:iK+mf-1;
				data.Xs(:,data.ind_Xs{i}{end}) = Gf;
				iK = iK+mf;
			else
				% last kernel for test data! This compute testing data
				K = exp( - b * sqdist(X(:,i)',Xtrain(datatrain.Is{i},i)' ) ) - data.Xs(:,fulldataloc)  * datatrain.Xs(datatrain.Is{i},fulldataloc)';
				data.d_Xs{i}(end) = datatrain.d_Xs{i}(end);
				data.ind_Xs{i}{end} = datatrain.ind_Xs{i}{end};
				data.Xs(:,data.ind_Xs{i}{end}) = K * inv( datatrain.Xs(datatrain.Is{i},data.ind_Xs{i}{end})' ) ;
				iK = iK+data.d_Xs{i}(end);

			end
		end






	case {'anova', 'anova-mkl', 'anova-bimkl'}
		% simply the same as 'gauss-hermite', but with the full kernels at
		% the end
		[n , p ] = size(X);
		b = kernel_params(1);

		data.n = n ;
		data.p = p;
		data.q = 1;
		data.qs = 2 * ones(1,p);     % number of kernels per directions
		data.Xs = zeros(n,sum(data.qs));      % storing all dimensions in the same vector
		data.ind_Xs = cell(1,p);
		data.d_Xs = cell(1,p);
		iK = 1;
		for i=1:p
			% constant term
			data.ind_Xs{i} = cell(1,2);
			data.d_Xs{i} = zeros(1,2);
			data.ind_Xs{i}{1} = iK;
			data.d_Xs{i}(1) = 1;
			data.Xs(:,iK) = ones(size(X,1),1);
			iK = iK + 1 ;

			if nargin < 4
				% first build the kernel matrix, then takes it
				% incomplete Cholesky decomposition

				K = exp( - b * sqdist(X(:,i)',X(:,i)' ) );
				[Gf,Pf,mf,residualf] = icd_general(K,1e-3,kernel_params(4));
				If = Pf(1:mf);

				data.Is{i} = If;					% these two for test data
				[temp,Pif] = sort(Pf);
				Gf = Gf(Pif,1:mf);
				data.d_Xs{i}(2) = mf;
				data.ind_Xs{i}{2} = iK:iK+mf-1;
				data.Xs(:,data.ind_Xs{i}{2}) = Gf;
				iK = iK+mf;
			else
				% last kernel for test!
				K = exp( - b * sqdist(X(:,i)',Xtrain(datatrain.Is{i},i)' ) );
				data.d_Xs{i}(2) = datatrain.d_Xs{i}(2);
				data.ind_Xs{i}{2} = datatrain.ind_Xs{i}{2};
				data.Xs(:,data.ind_Xs{i}{2}) = K * inv( datatrain.Xs(datatrain.Is{i},data.ind_Xs{i}{2})' ) ;
				iK = iK+data.d_Xs{i}(2);

			end
		end





	case {'base kernels', 'base kernels-mkl', 'base kernels-bimkl'}
		p = kernel_params(1);
		if nargin<4
			% train
			n = size(X,1);
			n = floor( sqrt(2 * n) );
		else
			% test
			ntrain = size(Xtrain,1);
			ntrain = floor( sqrt(2 * ntrain) );
			n = size(X,1) / ntrain;
		end
		q = kernel_params(2);
		data.n = n ;
		data.p = p;
		data.q = kernel_params(2);
		data.qs = (q+1) * ones(1,p);			% number of kernels per directions
		data.Xs = zeros(n,sum(data.qs));		% storing all dimensions in the same vector
		data.ind_Xs = cell(1,p);
		data.d_Xs = cell(1,p);
		iK = 1;
		for i=1:p
			% constant term
			data.ind_Xs{i} = cell(1,2);
			data.d_Xs{i} = zeros(1,2);
			data.ind_Xs{i}{1} = iK;
			data.d_Xs{i}(1) = 1;
			data.Xs(:,iK) = ones(n,1);
			iK = iK + 1 ;

			if nargin < 4
				% first build the kernel matrix, then takes it
				% incomplete Cholesky decomposition
				for j=1:q
					K = double( devectorize_single(X(:,i,j)) );

					[Gf,Pf,mf,residualf] = icd_general(K,1e-3,kernel_params(5));
					If = Pf(1:mf);

					data.Is{i}{j+1} = If;					% these two for test data
					[temp,Pif] = sort(Pf);
					Gf = Gf(Pif,1:mf);
					data.d_Xs{i}(j+1) = mf;
					data.ind_Xs{i}{j+1} = iK:iK+mf-1;

					data.Xs(:,data.ind_Xs{i}{j+1}) = Gf;
					iK = iK+mf;
				end
			else
				for j=1:q

					K = double( reshape(X(:,i,j) ,n,ntrain ) );
					K = K(:,datatrain.Is{i}{j+1});
					data.d_Xs{i}(j+1) = datatrain.d_Xs{i}(j+1);
					data.ind_Xs{i}{j+1} = datatrain.ind_Xs{i}{j+1};
					data.Xs(:,data.ind_Xs{i}{j+1}) = K * inv( datatrain.Xs(datatrain.Is{i}{j+1},data.ind_Xs{i}{j+1})' ) ;
					iK = iK+data.d_Xs{i}(j+1);

				end
			end
		end

	case {'spline', 'spline-mkl', 'spline-bimkl'}
		% simply the same as 'gauss-hermite', but with the full kernels at
		% the end
		[n , p ] = size(X);
		data.n = n ;
		data.p = p;
		data.q = 2;
		data.qs = 3 * ones(1,p);     % number of kernels per directions
		data.Xs = zeros(n,sum(data.qs)*30);      % storing all dimensions in the same vector
		data.ind_Xs = cell(1,p);
		data.d_Xs = cell(1,p);
		iK = 1;
		for i=1:p
			% constant term
			data.ind_Xs{i} = cell(1,2);
			data.d_Xs{i} = zeros(1,2);

			data.ind_Xs{i}{1} = iK;
			data.d_Xs{i}(1) = 1;
			data.Xs(:,iK) = ones(size(X,1),1);
			iK = iK + 1 ;

			data.ind_Xs{i}{2} = iK;
			data.d_Xs{i}(2) = 1;
			data.Xs(:,iK) = X(:,i);
			iK = iK + 1 ;


			if nargin < 4
				% first build the kernel matrix, then takes it
				% incomplete Cholesky decomposition
				K = compute_cubic_spline_kernel(X(:,i),X(:,i));
				[Gf,Pf,mf,residualf] = icd_general(K,1e-3,kernel_params(3));
				If = Pf(1:mf);

				data.Is{i} = If;					% these two for test data
				[temp,Pif] = sort(Pf);
				Gf = Gf(Pif,1:mf);
				data.d_Xs{i}(3) = mf;
				data.ind_Xs{i}{3} = iK:iK+mf-1;
				data.Xs(:,data.ind_Xs{i}{3}) = Gf;
				iK = iK+mf;
			else
				% last kernel for test!
				K = compute_cubic_spline_kernel(X(:,i),Xtrain(datatrain.Is{i},i));
				data.d_Xs{i}(3) = datatrain.d_Xs{i}(3);
				data.ind_Xs{i}{3} = datatrain.ind_Xs{i}{3};
				data.Xs(:,data.ind_Xs{i}{3}) = K * inv( datatrain.Xs(datatrain.Is{i},data.ind_Xs{i}{3})' ) ;
				iK = iK+data.d_Xs{i}(3);

			end
		end

		data.Xs(:,iK:end) = [];

end



