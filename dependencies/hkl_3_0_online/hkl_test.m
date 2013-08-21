function accuracies = hkl_test(model,outputs,Xtest,varargin);
%%%%%%%%%%%%%%%%%%%%%%%
% HKL
%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% required parameters
%
% Xtest					test data (input)
%
% optional parameters
%
% Ytest					test data (output), to compute predictive performance

Ytest = [];                         % test data (response)

% READ OPTIONAL PARAMETERS
args = varargin;
nargs = length(args);
for i=1:2:nargs
	switch args{i},
		case 'Ytest',           Ytest = args{i+1};
	end
end

data = model.data;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NORMALIZE DATA (in input space)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~strcmp(model.kernel,'base kernels') || ~strcmp(model.kernel,'base kernels-mkl') || ~strcmp(model.kernel,'base kernels-bimkl')
    switch model.data_normalization
        case  'scale'
            [ntest , p ] = size(Xtest);
            Xtest = Xtest - repmat(model.mean,ntest,1);
            Xtest = Xtest ./ repmat(model.std,ntest,1);
        case  'center'
            [ntest , p ] = size(Xtest);
            Xtest = Xtest - repmat(model.mean,ntest,1);
    end
end




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
switch model.kernel

	case {'polynomial','polynomial-mkl','polynomial-bimkl'}

		datatest = get_kernel_data(Xtest,model.kernel,model.kernel_params);
		datatest.dag_type =  'grid + input_space';

	case {'hermite', 'hermite-mkl', 'hermite-bimkl'}

		datatest = get_kernel_data(Xtest,model.kernel,model.kernel_params);
		datatest.dag_type =  'grid + input_space';

	case {'anova','anova-mkl','anova-bimkl'}

		datatest = get_kernel_data(Xtest,model.kernel,model.kernel_params,model.X,model.data);
		datatest.dag_type =  'grid + input_space';

	case {'spline','spline-mkl','spline-bimkl'}

		datatest = get_kernel_data(Xtest,model.kernel,model.kernel_params,model.X,model.data);
		datatest.dag_type =  'grid + input_space';
	
    case {'base kernels','base kernels-mkl','base kernels-bimkl'}

		datatest = get_kernel_data(Xtest,model.kernel,model.kernel_params,model.X,model.data);
		datatest.dag_type =  'grid + input_space';

	case {'gauss-hermite','gauss-hermite-mkl','gauss-hermite-bimkl'}

		datatest = get_kernel_data(Xtest,model.kernel,model.kernel_params);
		datatest.dag_type =  'grid + input_space';

	case {'gauss-hermite-full','gauss-hermite-full-mkl','gauss-hermite-full-bimkl'}

		datatest = get_kernel_data(Xtest,model.kernel,model.kernel_params,model.X,model.data);
		datatest.dag_type =  'grid + input_space';



end

switch model.kernel

	case {'hermite-mkl', 'gauss-hermite-mkl', 'gauss-hermite-full-mkl', 'polynomial-mkl', 'anova-mkl',  'spline-mkl' }
		dag_type = 'mkl + input_space';
		data.dag_type = dag_type;
		datatest.dag_type = dag_type;
	case {'hermite-bimkl', 'gauss-hermite-bimkl', 'gauss-hermite-full-bimkl', 'polynomial-bimkl', 'anova-bimkl',  'spline-bimkl' }
		dag_type = 'bimkl + input_space';
		data.dag_type = dag_type;
		datatest.dag_type = dag_type;
end

% NORMALIZATION OF BASE KERNELS
switch model.kernel_normalization
	case 'center'  % center before taking product of base kernels
		switch datatest.dag_type
			case { 'grid + input_space', 'mkl + input_space', 'bimkl + input_space'}
				for i=1:p
					for j=2:data.q+1
						datatest.Xs(:,data.ind_Xs{i}{j}) = datatest.Xs(:,data.ind_Xs{i}{j}) - repmat( model.meanXs{i}{j},n,1);
					end
				end

		end

	case 'scale'  % scale (and center) before taking product of base kernels
		switch datatest.dag_type
			case { 'grid + input_space', 'mkl + input_space', 'bimkl + input_space'}
				for i=1:p
					for j=2:data.q+1
						datatest.Xs(:,data.ind_Xs{i}{j}) = datatest.Xs(:,data.ind_Xs{i}{j}) - repmat( model.meanXs{i}{j},ntest,1);
						datatest.Xs(:,data.ind_Xs{i}{j}) = datatest.Xs(:,data.ind_Xs{i}{j}) / model.stdXs{i}{j};
					end
				end

		end


end



accuracies.predtest = NaN*ones(length(Ytest),length(outputs));
accuracies.predtest_unreg = NaN*ones(length(Ytest),length(outputs));

if ~isempty(Ytest)
	accuracies.testing_error = NaN*ones(1,length(outputs));
	accuracies.testing_error_unreg = NaN*ones(1,length(outputs));
	switch model.loss
		case 'logistic'
			accuracies.testing_error_class = NaN*ones(1,length(outputs));
			accuracies.testing_error_class_unreg = NaN*ones(1,length(outputs));
	end
end

for ilambda = 1:length(outputs)

	[predtrain,predtest] = predict_train_test(datatest,model.data,outputs{ilambda}.hull,outputs{ilambda}.alpha,outputs{ilambda}.b,outputs{ilambda}.zeta);
	[predtrain_unreg,predtest_unreg] = predict_train_test_unreg(datatest,model.data,model.Y,Ytest,outputs{ilambda}.hull);
	accuracies.predtest(:,ilambda) = predtest;
	accuracies.predtrain(:,ilambda) = predtrain;
	accuracies.predtest_unreg(:,ilambda) = predtest_unreg;
	accuracies.predtrain_unreg(:,ilambda) = predtrain_unreg;
	if ~isempty(Ytest)
		switch model.loss
			case 'square'
				accuracies.training_error(ilambda) =  sum( ( model.Y - predtrain ).^2 ) /length(model.Y);
				accuracies.testing_error(ilambda) =  sum( ( Ytest - predtest ).^2 ) /length(Ytest);
				accuracies.training_error_unreg(ilambda) =  sum( ( model.Y - predtrain_unreg ).^2 ) /length(model.Y);
				accuracies.testing_error_unreg(ilambda) =  sum( ( Ytest - predtest_unreg ).^2 ) /length(Ytest);

			case 'logistic'
				accuracies.training_error_class(ilambda) =  sum( abs( model.Y - (sign(predtrain-.5)+1)/2 ) ) /length(model.Y);
				accuracies.testing_error_class(ilambda) =  sum( abs( Ytest - (sign(predtest-.5)+1)/2 ) ) /length(Ytest);
				accuracies.training_error_class_unreg(ilambda) =  sum( abs( model.Y - (sign(predtrain_unreg-.5)+1)/2 ) ) /length(model.Y);
				accuracies.testing_error_class_unreg(ilambda) =  sum( abs( Ytest - (sign(predtest_unreg-.5)+1)/2 ) ) /length(Ytest);
				accuracies.training_error(ilambda) =  sum(  model.Y .* log( 1 + exp( -predtrain) ) + ( 1 - model.Y ) .* log( 1 + exp(predtrain) )  ) /length(model.Y);
				accuracies.testing_error(ilambda) =  sum(  Ytest .* log( 1 + exp( -predtest) ) + ( 1 - Ytest ) .* log( 1 + exp(predtest) )  ) /length(Ytest);
				accuracies.training_error_unreg(ilambda) =  sum(  model.Y .* log( 1 + exp( -predtrain_unreg) ) + ( 1 - model.Y ) .* log( 1 + exp(predtrain_unreg) )  ) /length(model.Y);
				accuracies.testing_error_unreg(ilambda) =  sum(  Ytest .* log( 1 + exp( -predtest_unreg) ) + ( 1 - Ytest ) .* log( 1 + exp(predtest_unreg) )  ) /length(Ytest);
		end

	end



end