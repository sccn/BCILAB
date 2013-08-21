function [outputs,model,accuracies] = hkl(X,Y,lambdas,loss,kernel,kernel_params,varargin);
%%%%%%%%%%%%%%%%%%%%%%%
% HKL
%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% required parameter
% X                     input data: n x p matrix (n=number of observations)
%						or kernel matrices ( n(n+1)/2 x p x q single)
% Y                     responses ( n x 1 matrix )
%                       NB: for classification in {0,1} (and not in {-1,1})
% lambdas               regularization parameters (might be vector or single number)
%						better to start first with large values of lambdas
% loss                  loss function ('square','logistic')
% kernel                kernel to be decomposed (see list below)
% kernel_params         kernel parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KERNELS
%
% directed grids:  note that adding '-mkl' or '-bimkl' reduces the directed
% grids to depth 1 (regular MKL) or depth 2 (almost MKL on all pairs)
%
%
%						'base kernels': given basis kernels
%						1: number of dimensions
%						2: maximal depth
%						3: weight of the root kernel (constant)
%						4: weight w of the kernels of depth 1 (all other kernels have weight w^d where d is depth)
%						5: maximal order of the incomplete Cholesky decomposition
%
%  						
%  						'polynomial': polynomial kernel. Good values: [  3 .1 4 ]
%						1: order of the polynomial
%						2: weight of the root kernel (constant)
%						3: weight w of the kernels of depth 1 (all other kernels have weight w^d where d is depth)
%
%  						'hermite': hermite polynomial kernel. Good values: [ .5 3 .1 4 ]
%						1: scale of hermite polynomial
%						2: order of the polynomial
%						3: weight of the root kernel (constant)
%						4: weight w of the kernels of depth 1 (all other kernels have weight w^d where d is depth)
%
%  						'gauss-hermite': hermite expansion of Gaussian kernel (finite order). Good values: [ 1 .05 3 .1 .5 ]
%						1: scale of Gaussian kernel expansion (corresponds to P(x)~N(0,sigma^2)
%						2: weight of the Gaussian kernel
%						3: order of the expansion
%						4: weight of the root kernel (constant)
%						5: weight w of the kernels of depth 1 (all other kernels have weight w^d where d is depth)
%
%  						'gauss-hermite-full': hermite expansion of Gaussian kernel (last order includes the entire remainder). Good values: [ 1 .05 3 .1 .5 30 ]
%						1: scale of Gaussian kernel expansion (corresponds to P(x)~N(0,sigma^2)
%						2: weight of the Gaussian kernel
%						3: order of the finite expansion
%						4: weight of the root kernel (constant)
%						5: weight w of the kernels of depth 1 (all other kernels have weight w^d where d is depth)
%						6: maximal order of the incomplete Cholesky decomposition
%
%  						'anova': anova kernel (all subsets). Good values: [ .05 .1 4 30]
%						1: weight of the Gaussian kernel
%                       2: weight of the root kernel (constant)
%						3: weight w of the kernels of depth 1 (all other kernels have weight w^d where d is depth)
%						4: maximal order of the incomplete Cholesky decomposition
%
%  						'spline': spline kernel. Good values: [ .1 4 40]
%						1: weight of the root kernel (constant)
%						2: weight w of the kernels of depth 1 (all other kernels have weight w^d where d is depth)
%						3: maximal order of the incomplete Cholesky decomposition
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% optional parameters
%
% mingap				duality gap for checking optimality (default=1e-3)
% gapeta				smoothing of the reduced problem (default=1e-3)
% maxactive				maximum number of selected kernel (degault = 100)
% display				1 if displaying progress of algorithm, 0 otherwise
% Xtest					test data (input)
% Ytest					test data (output), to compute predictive performance
% data_normalization    'center': original vectorial data is centered
%						'scale': data centered and scaled to unit variance (default)
% kernel_normalization  'center': feature space data of based kernel is centered (default)
%						'scale': feature space data of based kernel is scaled
% alpha					dual parameters
% eta					kernel weights
% b						constant term
% memory_cache			available memory in bytes (default=1e9(
%						this has a significant impact on performance, i.e., do use your memory!
% conjgrad              1: uses conjugate gradient
%                       0: does not (default)
%
%
%
%
%

mingap = 1e-3;                      % required duality gap
maxactive = 100;                    % maximum number of selected kernels
display = 0;                        % 1: if display, 0 otherwise
Xtest = [];                         % test data (input)
Ytest = [];                         % test data (response)
data_normalization = 'scale';        % type of normalization for the input data
kernel_normalization = 'scale-root';      % type of normalization for the feature space
hull = [];                          % initialization of hull
alpha = [];                         % initialization of alpha
b = [];                             % initialization of b
eta = [];                             % initialization of eta
gapeta = 1e-3;
memory_cache = 1e9;
conjgrad = 0;

% READ OPTIONAL PARAMETERS
args = varargin;
argstopass = []; % arguments to pass to later functions (solve_hkl)
nargs = length(args);
for i=1:2:nargs
	switch args{i},
		case 'display',         display = args{i+1}; argstopass = [ argstopass, i, i+1];
		case 'Xtest',           Xtest = args{i+1};
		case 'Ytest',           Ytest = args{i+1};
		case 'mingap',          mingap = args{i+1};  argstopass = [ argstopass, i, i+1];
		case 'maxactive',       maxactive = args{i+1};  argstopass = [ argstopass, i, i+1];
		case 'data_normalization',            data_normalization = args{i+1};
		case 'kernel_normalization',            kernel_normalization = args{i+1};
		case 'hull',            hull = args{i+1};
		case 'alpha',        	alpha = args{i+1};
		case 'b',				b = args{i+1};
		case 'eta',				eta = args{i+1};
		case 'memory_cache',   	memory_cache = args{i+1};
		case 'gapeta',        	gapeta = args{i+1};  argstopass = [ argstopass, i, i+1];
		case 'conjgrad',        	conjgrad = args{i+1};  argstopass = [ argstopass, i, i+1];
	end
end
argstopass = args(argstopass);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NORMALIZE DATA (in input space)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~strcmp(kernel,'base kernels') || ~strcmp(kernel,'base kernels-mkl') || ~strcmp(kernel,'base kernels-bimkl')
	switch data_normalization
		case  'scale'
			[n , p ] = size(X);
			model.mean = mean(X,1);	% stored to allow testing on new data points
			X = X - repmat(model.mean,n,1);
			model.std = sqrt( mean(X.^2,1) );
			X = X ./ repmat(model.std,n,1);
			if ~isempty(Xtest)
				[ntest , p ] = size(Xtest);
				Xtest = Xtest - repmat(model.mean,ntest,1);
				Xtest = Xtest ./ repmat(model.std,ntest,1);
			end
		case  'center'
			[n , p ] = size(X);
			model.mean = mean(X,1);
			X = X - repmat(model.mean.X,n,1);
			if ~isempty(Xtest)
				[ntest , p ] = size(Xtest);
				Xtest = Xtest - repmat(model.mean,ntest,1);
			end
	end
end

model.X = X;		% this is the centered data (or the kernels
model.Y = Y;
model.loss = loss;
model.kernel = kernel;
model.kernel_params = kernel_params;
model.data_normalization = data_normalization;
model.kernel_normalization = kernel_normalization;







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREPROCESS DATA
% the structure data store all data necessary to compute the kernel for a given node
%
% for directed grids where each kernel in the graph is the product of base kernels
% each base kernels is represented through a Cholesky decomposition
%
% data.n				number of observations
% data.p				dimension of the grid
% data.q				maximal depth of the grid
% data.qs				number of kernels in each direction of the grid
% data.Xs				storing all dimensions in the same vector
% data.ind_Xs			indices corresponding to each base kernel
% data.d_Xs				rank of each base kernel
%
% for directed grids, the (unique root) is always the constant kernel (and thus zeroed out by centering)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch kernel

case {'base kernels','base kernels-mkl','base kernels-bimkl'}
	
		p = kernel_params(1);
		n = length(Y);
		dag_type = 'grid + input_space';
		data = get_kernel_data(X,kernel,kernel_params);
		data.weight0 = kernel_params(3);
		data.weightA = kernel_params(4);
		for i=1:p
			data.weights{i} =  data.weightA .^ ( 0:kernel_params(2) );
		end
		data.dag_type = dag_type;
		model.dag_type = dag_type;
		if ~isempty(Xtest)
			datatest = get_kernel_data(Xtest,kernel,kernel_params,X,data);
			datatest.dag_type = dag_type;
		end

	case {'polynomial','polynomial-mkl','polynomial-bimkl'}
		[n , p ] = size(X);
		dag_type = 'grid + input_space';
		data = get_kernel_data(X,kernel,kernel_params);
		data.weight0 = kernel_params(2);
		data.weightA = kernel_params(3);
		for i=1:p
			data.weights{i} =  data.weightA .^ ( 0:kernel_params(1) );
		end
		data.dag_type = dag_type;
		model.dag_type = dag_type;
		if ~isempty(Xtest)
			datatest = get_kernel_data(Xtest,kernel,kernel_params);
			datatest.dag_type = dag_type;
		end

	case {'hermite', 'hermite-mkl', 'hermite-bimkl'}
		[n , p ] = size(X);
		dag_type = 'grid + input_space';
		data = get_kernel_data(X,kernel,kernel_params);
		data.weight0 = kernel_params(3);
		data.weightA = kernel_params(4);
		for i=1:p
			data.weights{i} =  data.weightA .^ ( 0:kernel_params(2) );
		end
		data.dag_type = dag_type;
		model.dag_type = dag_type;
		if ~isempty(Xtest)
			datatest = get_kernel_data(Xtest,kernel,kernel_params);
			datatest.dag_type = dag_type;
		end


	case {'anova','anova-mkl','anova-bimkl'}
		[n , p ] = size(X);
		dag_type = 'grid + input_space';
		data = get_kernel_data(X,kernel,kernel_params);
		data.weight0 = kernel_params(2);
		data.weightA = kernel_params(3);
		for i=1:p
			data.weights{i} =  data.weightA .^ ( 0:1 );
		end
		data.dag_type = dag_type;
		model.dag_type = dag_type;
		if ~isempty(Xtest)
			datatest = get_kernel_data(Xtest,kernel,kernel_params,X,data);
			datatest.dag_type = dag_type;
		end



	case {'spline','spline-mkl','spline-bimkl'}
		[n , p ] = size(X);
		dag_type = 'grid + input_space';
		data = get_kernel_data(X,kernel,kernel_params);
		data.weight0 = kernel_params(1);
		data.weightA = kernel_params(2);
		for i=1:p
			data.weights{i} =  data.weightA .^ ( 0:2 );
		end
		data.dag_type = dag_type;
		model.dag_type = dag_type;
		if ~isempty(Xtest)
			datatest = get_kernel_data(Xtest,kernel,kernel_params,X,data);
			datatest.dag_type = dag_type;
		end




	case {'gauss-hermite','gauss-hermite-mkl','gauss-hermite-bimkl'}
		[n , p ] = size(X);
		dag_type = 'grid + input_space';
		data = get_kernel_data(X,kernel,kernel_params);
		data.weight0 = kernel_params(4);
		data.weightA = kernel_params(5);
		for i=1:p
			data.weights{i} =  data.weightA .^ ( 0:kernel_params(3) );
		end
		data.dag_type = dag_type;
		model.dag_type = dag_type;
		if ~isempty(Xtest)
			datatest = get_kernel_data(Xtest,kernel,kernel_params);
			datatest.dag_type = dag_type;
		end



	case {'gauss-hermite-full','gauss-hermite-full-mkl','gauss-hermite-full-bimkl'}
		[n , p ] = size(X);
		dag_type = 'grid + input_space';
		data = get_kernel_data(X,kernel,kernel_params);
		data.weight0 = kernel_params(4);
		data.weightA = kernel_params(5);
		for i=1:p
			data.weights{i} =  data.weightA .^ ( 0:(1+kernel_params(3)) );
		end
		data.dag_type = dag_type;
		model.dag_type = dag_type;
		if ~isempty(Xtest)
			datatest = get_kernel_data(Xtest,kernel,kernel_params,X,data);
			datatest.dag_type = dag_type;
		end


end
switch kernel

	case {'base-kernels-mkl','hermite-mkl', 'gauss-hermite-mkl', 'gauss-hermite-full-mkl', 'polynomial-mkl', 'anova-mkl',  'spline-mkl' }
		dag_type = 'mkl + input_space';
		data.dag_type = dag_type;
		datatest.dag_type = dag_type;
	case {'base-kernels-bimkl','hermite-bimkl', 'gauss-hermite-bimkl', 'gauss-hermite-full-bimkl', 'polynomial-bimkl', 'anova-bimkl',  'spline-bimkl' }
		dag_type = 'bimkl + input_space';
		data.dag_type = dag_type;
		datatest.dag_type = dag_type;
end


 
 

% NORMALIZATION OF BASE KERNELS
switch kernel_normalization
	case 'center'  % center before taking product of base kernels
		switch dag_type
			case { 'grid + input_space', 'mkl + input_space', 'bimkl + input_space'}
				for i=1:p
					for j=2:data.q+1
						model.meanXs{i}{j} =  mean(data.Xs(:,data.ind_Xs{i}{j}),1);
						data.Xs(:,data.ind_Xs{i}{j}) = data.Xs(:,data.ind_Xs{i}{j}) - repmat( model.meanXs{i}{j},n,1);
						if ~isempty(Xtest)
							datatest.Xs(:,data.ind_Xs{i}{j}) = datatest.Xs(:,data.ind_Xs{i}{j}) - repmat( model.meanXs{i}{j},n,1);
						end
					end
				end

		end

	case 'scale'  % scale (and center) before taking product of base kernels
		switch dag_type
			case { 'grid + input_space', 'mkl + input_space', 'bimkl + input_space'}
				for i=1:p
					for j=2:data.q+1
						model.meanXs{i}{j} =  mean(data.Xs(:,data.ind_Xs{i}{j}),1);
						data.Xs(:,data.ind_Xs{i}{j}) = data.Xs(:,data.ind_Xs{i}{j}) - repmat( model.meanXs{i}{j},n,1);
						model.stdXs{i}{j} =  sqrt(sum(mean((data.Xs(:,data.ind_Xs{i}{j})).^2,1)));
						data.Xs(:,data.ind_Xs{i}{j}) = data.Xs(:,data.ind_Xs{i}{j}) / model.stdXs{i}{j};

						if ~isempty(Xtest)
							datatest.Xs(:,data.ind_Xs{i}{j}) = datatest.Xs(:,data.ind_Xs{i}{j}) - repmat( model.meanXs{i}{j},ntest,1);
							datatest.Xs(:,data.ind_Xs{i}{j}) = datatest.Xs(:,data.ind_Xs{i}{j}) / model.stdXs{i}{j};
						end
					end
				end

		end


end
model.data = data;		% after normalization of the kernel




% INITIALIZE CACHE OF KERNEL MATRICES
switch dag_type
	case 'grid + input_space'
		global K_cache;			% cached values (single)
		global max_cache;		% maximum number of cached matrices
		global n_cache;			% current number of cached matrices
		global Kdir_caches;		% cache of kernels reused many times in checking sufficiet conditions
		global source_cache;	% caching which sources are used for each cached kernel
		n_cache = 0;

		source_cache = zeros(0,0,'int32');
		max_cache = round(memory_cache / data.n / (data.n+1)*2 / 4);
		K_cache = zeros(data.n*(data.n+1)/2,max_cache,'single');
		Kdir_cache = ones(data.n*(data.n+1)/2,1,'single');

		% compute each based kernel
		for i=1:data.p
			for jj=1:data.q+1
				temp = zeros(data.n,data.n,'single');
				cumweight = single( 1./cumsum(data.weights{i}(jj:end)).^2 )';

				for j=jj:data.q+1
					temp = temp   + cumweight(j-jj+1) * single( data.Xs(:,data.ind_Xs{i}{j}) * data.Xs(:,data.ind_Xs{i}{j})' );
				end
				Kdir_caches{i}{jj} = symmetric_vectorize_single(temp);
			end

			% caching the kernel for the root
			Kdir_cache = Kdir_caches{i}{1} .* Kdir_cache ;
		end
		source_cache(:,1) = ones(1,p,'int32');
		K_cache(:,1) = Kdir_cache;
		n_cache = 1;
end


outputs  = cell(1,length(lambdas));

accuracies.predtrain = NaN*ones(length(Y),length(lambdas));
accuracies.predtrain_unreg = NaN*ones(length(Y),length(lambdas));
accuracies.training_error = NaN*ones(1,length(lambdas));
accuracies.training_error_unreg = NaN*ones(1,length(lambdas));

switch loss
	case 'logistic'
		accuracies.training_error_class = NaN*ones(1,length(lambdas));
		accuracies.training_error_class_unreg = NaN*ones(1,length(lambdas));

end


if ~isempty(Xtest)

	accuracies.predtest = NaN*ones(length(Ytest),length(lambdas));
	accuracies.predtest_unreg = NaN*ones(length(Ytest),length(lambdas));

	if ~isempty(Ytest)
		accuracies.testing_error = NaN*ones(1,length(lambdas));
		accuracies.testing_error_unreg = NaN*ones(1,length(lambdas));
		switch loss
			case 'logistic'
				accuracies.testing_error_class = NaN*ones(1,length(lambdas));
				accuracies.testing_error_class_unreg = NaN*ones(1,length(lambdas));
		end
	end
end

for ilambda = 1:length(lambdas)
	% run HKL with warm restarts
	output = solve_hkl(data,Y,lambdas(ilambda),loss,'hull',hull,'alpha',alpha,'eta',eta,'b',b,argstopass{:});
	hull = output.hull;
	alpha = output.alpha;
	eta = output.eta;
	b = output.b;

	outputs{ilambda} = output;


	if ~isempty(Xtest)
		[predtrain,predtest] = predict_train_test(datatest,data,output.hull,output.alpha,output.b,output.zeta);
		[predtrain_unreg,predtest_unreg] = predict_train_test_unreg(datatest,data,Y,Ytest,output.hull);
		accuracies.predtest(:,ilambda) = predtest;
		accuracies.predtrain(:,ilambda) = predtrain;
		accuracies.predtest_unreg(:,ilambda) = predtest_unreg;
		accuracies.predtrain_unreg(:,ilambda) = predtrain_unreg;
		if ~isempty(Ytest)
			switch loss
				case 'square'
					accuracies.training_error(ilambda) =  sum( ( Y - predtrain ).^2 ) /length(Y);
					accuracies.testing_error(ilambda) =  sum( ( Ytest - predtest ).^2 ) /length(Ytest);
					accuracies.training_error_unreg(ilambda) =  sum( ( Y - predtrain_unreg ).^2 ) /length(Y);
					accuracies.testing_error_unreg(ilambda) =  sum( ( Ytest - predtest_unreg ).^2 ) /length(Ytest);

				case 'logistic'
					accuracies.training_error_class(ilambda) =  sum( abs( Y - (sign(predtrain-.5)+1)/2 ) ) /length(Y);
					accuracies.testing_error_class(ilambda) =  sum( abs( Ytest - (sign(predtest-.5)+1)/2 ) ) /length(Ytest);
					accuracies.training_error_class_unreg(ilambda) =  sum( abs( Y - (sign(predtrain_unreg-.5)+1)/2 ) ) /length(Y);
					accuracies.testing_error_class_unreg(ilambda) =  sum( abs( Ytest - (sign(predtest_unreg-.5)+1)/2 ) ) /length(Ytest);
					accuracies.training_error(ilambda) =  sum(  Y .* log( 1 + exp( -predtrain) ) + ( 1 - Y ) .* log( 1 + exp(predtrain) )  ) /length(Y);
					accuracies.testing_error(ilambda) =  sum(  Ytest .* log( 1 + exp( -predtest) ) + ( 1 - Ytest ) .* log( 1 + exp(predtest) )  ) /length(Ytest);
					accuracies.training_error_unreg(ilambda) =  sum(  Y .* log( 1 + exp( -predtrain_unreg) ) + ( 1 - Y ) .* log( 1 + exp(predtrain_unreg) )  ) /length(Y);
					accuracies.testing_error_unreg(ilambda) =  sum(  Ytest .* log( 1 + exp( -predtest_unreg) ) + ( 1 - Ytest ) .* log( 1 + exp(predtest_unreg) )  ) /length(Ytest);
			end
			
		end


	end
	active = find(output.eta .* output.weights.^2 - gapeta/  length(output.eta) >1e-12);


	if display
		fprintf('\nlambda = %f - # of active kernels = %d\n\n',lambdas(ilambda),length(active));
	end
end