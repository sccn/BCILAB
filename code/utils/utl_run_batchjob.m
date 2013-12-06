function results = utl_run_batchjob(opts,d,appname,setnames)
% Internal: the actual processing function of bci_batchtrain
% Results = ult_run_batchjob(Options,DatasetIndex,ApproachName,DatasetNames)
try
    results = [];    
    storename = env_translatepath(strrep(strrep(opts.storepatt,'%set',setnames{d}),'%approach',appname));
    if opts.reuse && ~isempty(storename) && exist(storename,'file')
        disp(['Reusing existing result for approach "' appname '" on set "' setnames{d} '".']);
        disp(['Trying to load file ' storename]);
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
                    disp(['Error computing predictions for set "' setnames{d} '", prediction set #' num2str(k) ' with approach "' appname '".']);
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
        eval([strrep(strrep(opts.resultpatt,'%num',num2str(d)),'%approach',appname) ' res;']); end
    try
        if isfield(res,'pred_loss')
            fprintf('%s@%s: train loss: %.4f, prediction loss: %.4f\n',appname,setnames{d},res.loss,res.pred_loss);
        else
            fprintf('%s@%s: train loss: %.4f\n',appname,setnames{d},res.loss);
        end
    catch
        disp('Could not display loss estimates.');
    end
catch e
    disp(['Error processing data set "' setnames{d} '" with approach "' appname '".']);
    opts.handler(e);
end
