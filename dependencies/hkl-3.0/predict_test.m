function [predtrain,predtest] = predict_test(datatest,hull,alpha,b,zeta);



switch datatest.dag_type
	case {'grid + input_space' , 'mkl + input_space', 'bimkl + input_space' }

 		% prepare data
		predtest = zeros(datatest.n,1)+b;
		ai1=1;
		for i1=1:size(hull,1)
			Xloctest = get_data_reduced(hull(i1,:),datatest);
			predtest = predtest + zeta(i1) * Xloctest * ( Xloctrain' * alpha );
		end
end




















