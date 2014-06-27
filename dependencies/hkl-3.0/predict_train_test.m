
function          [predtrain,predtest] = predict_train_test(datatest,datatrain,hull,alpha,b,zeta);



switch datatrain.dag_type
    case {'grid + input_space' , 'mkl + input_space', 'bimkl + input_space' }





        % prepare data
        predtrain = zeros(datatrain.n,1)+b;
        predtest = zeros(datatest.n,1)+b;
        ai1=1;
        for i1=1:size(hull,1)
            Xloctrain = get_data_reduced(hull(i1,:),datatrain);
            Xloctest = get_data_reduced(hull(i1,:),datatest);
            predtrain = predtrain + zeta(i1) * Xloctrain * ( Xloctrain' * alpha );
            predtest = predtest + zeta(i1) * Xloctest * ( Xloctrain' * alpha );
        end


end




















