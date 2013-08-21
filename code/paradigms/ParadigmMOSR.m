classdef ParadigmMOSR < ParadigmBase
    % Multi-subject Overcomplete Spectral Regression
    %
    % This method is built for oscillations that are stationary in the time window of interest.
    % It relies on an AMICA (or other) decomposition of the data to get spatially filtered source signals,
    % then performs a multi-taper spectral estimation on the resulting components and a per-component PCA
    % of the spectra (to reduce dimensionality). The actual classifier that operates on these features is 
    % by default Least-Angle regression (LARS) due to its speed. Note that this is a sparse feature-selecting
    % classifier which can effectively deal with an arbitrary number of features (as long a subset of these contain
    % the information of interest). The method has not yet been optimized for speed, so is relatively slow to train
    % currently (esp. the ICA part).
    %
    % The multi-subject aspect amounts to using observations from a pool of other subjects and computing a prior distribution
    % over model parameters on this pool. The resulting distribution is then used as a prior when learning from the calibration
    % set for the specific person (or session) of interest. Thus, this is a simple hierarchical Bayesian model.
    %
    % References:
    %  C. A. Kothe and S. Makeig, 
    %  “Estimation of Task Workload from EEG Data: New and Current Tools and Perspectives,” 
    %  IEEE EMBC, vol. 2011, pp. 6547-6551, 2011.
    % 
    % Name:
    %   Multi-subject Overcomplete Spectral Regression, work in progress
    %
    
    methods
        
        function defaults = preprocessing_defaults(self)
            % define the default pre-processing parameters for the M-OSR paradigm
            defaults = { ...
                 'EOGRemoval',{'RemoveEOG',true}, ...
                 'DipoleFitting','on',...
                 'ICA',{{'amica', ...
                         'version','stable11', ...
                         'num_models',3, ...
                         'useqsub','on', ...
                         'max_iter',500, ...
                         'max_init_waiting',2000, ...
                         'max_restarts',20, ...
                         'fallback_reduce',0.3}, 'clean',{'noisy'}}, ...
                 'Projection','on', ...
                 'EpochExtraction',{'TimeWindow',[-15 15]},...
                 'SpectralTransform',{'Representation',{'multitaper','bandwidth',3},'LogTransform',true,'LogSpacing',200}, ...
                 'EpochPCA',{'RetainDimensions',20}};
        end
        
        
        function defaults = machine_learning_defaults(self)
            % set up the default parameters for machine learning
            defaults = {'logreg',[],'variant','lars'};
        end
        
        
        function model = calibrate(self,varargin)
            % calibrate an M-OSR model from a corpus of training sets
            args = arg_define(varargin, ...
                arg_norep({'collection','Collection'}), ...
                arg_norep({'goal_identifier','GoalIdentifier'}), ...
                arg({'variant','Variant'},'pool',{'pool','weighted','jointprior','sparse_jointprob'},'Variant to use.'), ...
                arg({'alpha','ElasticMixing'},1,[0.01 1],'ElasticNet mixing parameter. The default is the lasso penalty.'), ...
                arg({'ref_weight','ReferenceWeight'},1,[],'Reference weight. Relative weight of the reference set.'), ...
                arg({'scaling','Scaling'}, 'std', {'none','center','std','minmax','whiten'}, 'Pre-scaling of the data. For the regulariation to work best, the features should either be naturally scaled well, or be artificially scaled.'),...
                arg({'subsample','Subsampling'}, 3, [], 'Sub-sampling of the data. Larger means more samples/trials skipped.'),...
                arg({'testwindow','TestWindow'}, [-15 15], [], 'Online window length. This is the window length used for test-set prediction.'),...
                arg({'nfolds','NumFolds'},'auto',[],'Cross-validation folds. The cross-validation is used to determine the best regularization parameter.'),...
                arg({'priorsolve','PriorSolver'}, 'lax', {'lax','precise','sparse'}, 'Prior Solver. Determines how the prior is derived; can be lax, precise, or smooth. Precise is > 30x slower than lax.'),...
                arg_sub({'flt','SignalProcessing'}, self.preprocessing_defaults(), @flt_pipeline, 'Signal processing stages. These parameters control filter stages that run on the signal level; they can be enabled, disabled and configured for the given paradigm. The prediction operates on the outputs of this stage.'), ...
                arg_sub({'ml','MachineLearning'},{'Learner',self.machine_learning_defaults()},@ml_train,'Machine learning stage of the paradigm. Operates on the feature vectors that are produced by the feature-extraction stage.'),...
                arg({'arg_dialogsel','ConfigLayout'},self.dialog_layout_defaults(),[],'Parameters displayed in the config dialog. Cell array of parameter names to display (dot-notation allowed); blanks are translated into empty rows in the dialog. Referring to a structure argument lists all parameters of that struture, except if it is a switchable structure - in this case, a pulldown menu with switch options is displayed.','type','cellstr','shape','row'));
            
            % do a full ICA decomposition on the set from the collection that is closest to the target set
            % and keep all others that are in the "closest set"
            [reference,remaining] = utl_collection_closest(args.collection,args.goal_identifier); 
            remaining_ref = reference(1:end-1); reference = reference{end};
            reference.streams{1} = exp_eval_optimized(flt_pipeline('signal',reference.streams{1}, args.flt)); %#ok<*NODEF>
            
            % apply that pipeline to all other reference sets
            for k=1:length(remaining_ref)
                remaining_ref{k} = self.dummy_preprocess(remaining_ref{k},reference.streams{1}.tracking.online_expression); end
            % ... and to all testing sets in the collection
            for k=1:length(remaining)
                remaining{k} = self.dummy_preprocess(remaining{k},reference.streams{1}.tracking.online_expression); end
            
            % rebuild the reference set
            reference_set = [remaining_ref {reference}];
            
            % pool all sets into one uniform corpus
            corpus = [reference_set remaining];
            
            % extract features from all sets
            for k=1:length(corpus)
                signal = corpus{k}.streams{1};
                features{k} = squeeze(reshape(signal.data,[],1,size(signal.data,3)))';
                targets{k} = set_gettarget(signal);
                if k <= length(reference_set)
                    weights{k} = ones(size(targets{k}));  % ref weights pos
                else
                    weights{k} = -ones(size(targets{k})); % otehr weights neg
                end
            end
            
            % concatenate them, and derive what was what
            features = vertcat(features{:});
            targets = vertcat(targets{:});
            weights = vertcat(weights{:});
            
            % sub-sample them optionally
            features = features(1:args.subsample:end,:);
            targets = targets(1:args.subsample:end);
            weights = weights(1:args.subsample:end);
            
            refs = weights == 1;
            
            switch args.variant
                case 'pool'
                    % throw a machine learning method at them...
                    model.predictivemodel = ml_train('data',{features,targets}, args.ml);
                    
                case 'weighted'
                    % compute proper weights
                    weights = refs/nnz(refs) * args.ref_weight + ~refs/nnz(~refs);                    
            
                    % throw a machine learning method at them...
                    model.predictivemodel = ml_train('data',{features,targets,weights}, 'Learner',{'logreg',0.1,'variant','l1'});
                            
                case 'sparse_jointprior'
                    % find out what the reference set is sorted by
                    if strcmp(args.nfolds,'auto')
                        setfields = cellfun(@(x)fieldnames(x),reference_set,'UniformOutput',false);
                        allfields = setfields{1};
                        for k=2:length(setfields)
                            allfields = intersect(allfields,setfields{k}); end
                        allfields = setdiff(allfields,{'streams'});
                        for f=1:length(allfields)
                            fn = allfields{f};
                            vals = cellfun(@(x)x.(fn),reference_set,'UniformOutput',false);
                            if all(cellfun('isreal',vals))
                                vals = cell2mat(vals); end
                            try
                                if issorted(vals) && length(unique(vals)) > 1
                                    args.nfolds = length(unique(vals)); 
                                    break;
                                end
                            catch
                            end
                        end
                        if ischar(args.nfolds)
                            args.nfolds = min(5,length(reference_set)); end % fall back to 5-fold
                    end
                    
                    % first scale the data
                    model.scaling = hlp_findscaling(features,args.scaling);
                    features = hlp_applyscaling(features,model.scaling);
                    
                    % compute a Gaussian prior over the other data sets' models
                    switch args.priorsolve
                        % on the tractable set, the vb mode takes 44s and the vb-iter mode taketh 1388 s = 23 minutes; this is fast enough for that set, but would be prohib slow on the full set (likely >10x as slow)
                        case 'lax'
                            prior = ml_train('data',{features(~refs,:),targets(~refs)},'Learner',{'logreg',[],'variant','vb','scaling','none'}); 
                        case 'precise'
                            prior = ml_train('data',{features(~refs,:),targets(~refs)},'Learner',{'logreg',[],'variant','vb-iter','scaling','none'});
                        case 'sparse'
                            prior = ml_train('data',{features(~refs,:),targets(~refs)},'Learner',{'logreg',[],'variant','vb-ard','scaling','none'});
                    end
                    mu0 = prior.model.w(1:end-1);
                    sig0 = prior.model.V(1:end-1,1:end-1);

                    % get the correct penalty factor (including the prior variance, bias and the prior mean)
                    penalty_factor = [1./sqrt(diag(sig0)); 0; 0];

                    % expand the features by -features*mu0 and by a bias term (will be an unreg. mean term)
                    features = [features -features*mu0 ones(size(features,1),1)];
                    
                    % compute the model given the prior
                    model.predictivemodel = ml_train('data',{features,targets}, 'Learner',{'logreg',1/args.ref_weight,'variant',{'lars','nfolds',args.nfolds,'alpha',args.alpha,'penalty_factor',penalty_factor,'nlambda',100,'maxit',300},'scaling','none'});
                    model.mu0 = mu0;
            end
            
            % set the filter graph
            model.tracking.filter_graph = reference;
            model.variant = args.variant;
        end
        
        
        function outputs = predict(self,bundle,model)
            % predict given the extracted features and the model
            features = squeeze(reshape(bundle.streams{1}.data,[],1,size(bundle.streams{1}.data,3)))';
            
            if isfield(model,'variant') && (strcmp(model.variant,'jointprior') || strcmp(model.variant,'sparse_jointprior'))
                % apply the data standardization
                features = hlp_applyscaling(features,model.scaling);
                % expand feature space
                features = [features -features*model.mu0 ones(size(features,1),1)];
            end
            
            outputs = ml_predict(features, model.predictivemodel);
        end
        
        
        % temporary preprocessing function
        function bundle = dummy_preprocess(self,bundle,onlexp)
            % handle only the first stream: replace the @rawdata node by the actual data
            bundle.streams{1} = utl_replacerepeated(onlexp,{exp_rule(exp_blank(@rawdata),bundle.streams{1})});
            % then evaluate it with some level of caching enabled..
            bundle.streams{1} = exp_eval_optimized(bundle.streams{1});
        end
      
        function layout = dialog_layout_defaults(self)
            % define the default configuration dialog layout 
            layout = {'SignalProcessing.EOGRemoval.RemoveEOG', '', 'SignalProcessing.ICA.CleaningLevel.DataSetting', '',...
                'SignalProcessing.ICA.Variant.AmicaVersion', 'SignalProcessing.ICA.Variant.NumModels', ...
                'SignalProcessing.ICA.Variant.UseGridEngine', 'SignalProcessing.ICA.Variant.NumProcessors', ...
                'SignalProcessing.EpochExtraction','SignalProcessing.SpectralTransform.Representation.TimeBandwidth', ...
                'SignalProcessing.EpochPCA.RetainDimensions','', ...
                'Variant','PriorSolver','','MachineLearning.Learner'};
        end
                
    end
end
            
% (turn off a few editor warnings because some actual implementations are missing in this file)
%#ok<*INUSD,*STOUT,*MANU>
