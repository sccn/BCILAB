function results = utl_run_batchjob(opts,d,appname,setnames)
% Internal: the actual processing function of bci_batchtrain.
% Results = ult_run_batchjob(Options,DatasetIndex,ApproachName,DatasetNames)

try
    results = [];
    storename = env_translatepath(strrep(strrep(opts.storepatt,'%set',setnames{d}),'%approach',appname));
    if opts.reuse && ~isempty(storename) && exist(storename,'file')
        fprintf('Reusing existing result for approach "%s" on set "%s".\n',appname,setnames{d});
        fprintf('Trying to load file %s.\n',storename);
        io_load(storename);
    elseif ~opts.loadonly
        % train a model on the Dataset
        [res.loss,res.model,res.stats] = bci_train(opts.trainargs{:}, 'data',opts.datasets{d}, 'approach',opts.approaches.(appname), 'markers',opts.markers);
        
        if ~isempty(opts.predictsets) && ~isempty(opts.predictsets{d})
            % optionally run bci_predict on the PredictSets
            for k=1:length(opts.predictsets{d})
                try
                    [res.pred_predictions{k},res.pred_loss(k),res.pred_stats{k},res.pred_targets{k}] = bci_predict(opts.predictargs{:},'data',opts.predictsets{d}{k},'model',res.model, 'markers',opts.markers);
                catch e
                    [res.pred_predictions{k},res.pred_loss(k),res.pred_stats{k},res.pred_targets{k}] = deal([],NaN,struct(),[]);
                    fprintf('Error computing predictions for set "%s", prediction set #%i with approach "%s".\n',setnames{d},k,appname);
                    opts.handler(e);
                end
            end
            res.pred_loss = res.pred_loss';
        end
            
        % save results
        if ~isempty(storename)
            io_save(storename,opts.saveargs{:},'res'); end
    end
    if ~isempty(opts.resultpatt)
        try
            statement = [strrep(strrep(opts.resultpatt,'%num',num2str(d)),'%approach',appname) ' res;'];
            eval(statement); 
        catch e
            fprintf('Failed to evaluate ResultPattern (statement was "%s") with error: %s.\n',statement,e.message);
        end
    end
    try
        if isfield(res,'pred_loss')
            fprintf('%s@%s: train loss: %.4f, prediction loss: %.4f\n',appname,setnames{d},res.loss,res.pred_loss);
        else
            fprintf('%s@%s: train loss: %.4f\n',appname,setnames{d},res.loss);
        end
    catch e
        fprintf('Failed to display loss estimates with error: %s\n',e.message);
    end
catch e
    fprintf('Error processing data set "%s" with approach "%s".\n',setnames{d},appname);
    try
        opts.handler(e);
    catch he
       fprintf('Failed to evaluate result handler with error: %s\n.',he.message);
    end
end
