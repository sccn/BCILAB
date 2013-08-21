
function          [predtrain,predtest] = predict_train_test_unreg(datatest,datatrain,Y,Ytest,hull);



switch datatrain.dag_type
    case {'grid + input_space' , 'mkl + input_space' , 'bimkl + input_space'}


        Xloctrain = [];
        Xloctest = [];
        for i1=1:size(hull,1)
            Xloctrain = [ Xloctrain , get_data_reduced(hull(i1,:),datatrain)];
            Xloctest = [ Xloctest,  get_data_reduced(hull(i1,:),datatest) ];
        end

        if size(Xloctrain,2)>size(Xloctrain,1)/2
            % if too many features, do not do unregularized problem
            predtrain = Y*0 + mean(Y);
            predtest =  Ytest*0 + mean(Y);
        else
            meantrain = mean(Xloctrain,1);

            Xloctrain = Xloctrain - repmat(meantrain,size(Xloctrain,1),1);
            Xloctest = Xloctest - repmat(meantrain,size(Xloctest,1),1);


            W = (Xloctrain' * Xloctrain + 1e-10 * size(Xloctrain,2) * eye(size(Xloctrain,2)) ) \ ( Xloctrain' * (Y-mean(Y)) );
            predtrain = mean(Y) + Xloctrain * W;
            predtest = mean(Y) + Xloctest * W;
        end


end




















