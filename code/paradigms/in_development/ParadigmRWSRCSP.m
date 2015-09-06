classdef ParadigmRWSRCSP < ParadigmBase
    % Riemannian Weighted Subjects Regularized CSP (RWSRCSP)
    %
    % This paradigm implements RWSRCSP [1], which is a generalization of the Common Spatial Patterns
    % algorithms to calibration data comprising multiple subjects (or recordings). To train a model
    % for a particular "goal" (or target) subject using auxiliary data from other subjects, this
    % algorithm blends the subject-specific covariance matrices used in standard CSP and LDA with
    % weighted averages of covariance matrices from other subjects. The weighting is determined
    % based on the riemannian distance between the target subject and each respective other subject.
    % Also, the method can optionally average predictions obtained with different values of the 
    % regularization parameter instead of doing a (costly) parameter search as suggested in [1].
    % 
    % References:
    % [1] Lotte, F.
    %     "Signal processing approaches to minimize or suppress calibration time in oscillatory activity-based Brain-Computer Interfaces", 
    %     Proceedings of the IEEE, vol. 103, no. 6, pp. 871-890, 2015
    %
    % Name:
    %   Riemannian Weighted Subjects Regularized Common Spatial Patterns
    %
    %                            Christian Kothe, Syntrogi
    %                            2015-07-24
    
    methods
        
        function defaults = preprocessing_defaults(self)
            defaults = {'FIRFilter',{'Frequencies',[6 8 28 32],'Type','minimum-phase'}, 'EpochExtraction',[0.5 3.5], 'Resampling',100};
        end
                
        function defaults = machine_learning_defaults(self)
            % set up the default parameters for machine learning; not necessary
            defaults = {'lda'};
        end
                
        function model = calibrate(self,varargin)
            % calibrate an SSRCSP model from a corpus of training sets
            args = arg_define(varargin, ...
                arg_norep({'collection','Collection'}), ...
                arg_norep({'goal_identifier','GoalIdentifier'}), ...
                arg({'patterns','PatternPairs'},3,uint32([1 1 64 10000]),'Number of CSP patterns (times two).'),...
                arg({'lambdas','Lambdas'},0.1:0.1:0.9,[0 1],'Covariance shrinkage. A range of shrinkage parameters to run over (classifiers will be averaged).'),...                
                arg({'group_tasks_by','GroupTasksBy'},'subject',{'group','subject','day','montage','session','recording','block'},'Group tasks into. This allows to group the training data into tasks solved in a multi-task manner, e.g., such that data of a given subject forms one task. When hyper-parameters need to be optimized, this would usually be done using a basic blockwise cross-validation within each task.'), ...
                arg({'shrinkage_cov','ShrinkageCovariance','ShrinkageCov'},false,[],'Shrinkage covariance estimator. Whether to use shrinkage to estimate the covariance matrices.'), ...
                arg({'weight_bias','WeightedBias'}, false, [], 'Account for class priors in bias. If you do have unequal probabilities for the different classes, this should be enabled.'), ...
                arg({'weight_cov','WeightedCov'}, false, [], 'Account for class priors in covariance. If you do have unequal probabilities for the different classes, it makes sense to enable this.'), ...
                arg({'apply_to','ApplyTo'},'channels',{'channels','sources','components','full CSD'},'Apply classifier to. Allows to select the type of time series to apply this model to.'), ...
                arg_sub({'flt','SignalProcessing'}, self.preprocessing_defaults(), @flt_pipeline, 'Signal processing stages. These parameters control filter stages that run on the signal level; they can be enabled, disabled and configured for the given paradigm. The prediction operates on the outputs of this stage.'), ...
                arg_sub({'ml','MachineLearning'},{'Learner',self.machine_learning_defaults()},@ml_train,'Machine learning stage of the paradigm. Operates on the feature vectors that are produced by the feature-extraction stage.'),...
                arg({'arg_dialogsel','ConfigLayout'},self.dialog_layout_defaults(),[],'Parameters displayed in the config dialog. Cell array of parameter names to display (dot-notation allowed); blanks are translated into empty rows in the dialog. Referring to a structure argument lists all parameters of that struture, except if it is a switchable structure - in this case, a pulldown menu with switch options is displayed.','type','cellstr','shape','row'));
           
            % if this is run on a worker, we'll set the cache capacity to zero since no machine
            % has enough RAM to hold multiple workers' copies of the corpus in memory
            on_worker = hlp_iscaller('par_worker');
            if on_worker
                global tracking;
                tracking.cache.capacity = 0; 
            end
            
            
            % split data into reference data (of goal subject) and remaining data
            [refsets,remaining] = utl_collection_closest(args.collection,args.goal_identifier);
            
            % recombine and move ref data to the beginning (because this paradigm will later 
            % extract the weights learned for the first task) and turn into struct array for
            % convenience
            corpus = [refsets{:} remaining{:}];
            
            % determine group membership
            group_membership = {corpus.(args.group_tasks_by)};
            if iscellstr(group_membership)
                groups = unique(group_membership);
            else
                groups = num2cell(unique([group_membership{:}]));
            end
            
            % initialize per-class covariance matrices
            CovCSP = deal(cell(2,length(groups)));
            CovLDA = deal(cell(2,length(groups)));
            % for each group...
            for s=length(groups):-1:1
                matches = find(cellfun(@(x)isequal(x,groups{s}),group_membership));
                matchdata = cell(1,length(matches));
                
                % collect all matching data sets...
                for p=1:length(matches)
                    matchdata{p} = corpus(matches(p));
                    if length(matchdata{p}.streams) > 1
                    disp_once('Note: ParadigmRWCSP will use only the first data stream of a recording (no support for multi-modal data).'); end
                    matchdata{p} = matchdata{p}.streams{1};
                end
                
                % concatenate them into a single set
				procdata = set_concat(matchdata{:});
				% and preprocess the result (with further settings/overrides according to args.flt)
                procdata = flt_pipeline('Signal',procdata, args.flt); %#ok<*NODEF>
                if on_worker
                    % if we're running on a worker we don't cache the result due to memory
                    % constraints
                    procdata = exp_eval(procdata);
                else
                    procdata = exp_eval_optimized(procdata);
                end                    
                % extract data
                switch args.apply_to
                    case 'channels'
                        X = procdata.data;
                    case 'components'
                        X = reshape((procdata.icaweights*procdata.icasphere)*procdata.data(procdata.icachansind,:),[],procdata.pnts,procdata.trials);
                    case 'sources'
                        X = procdata.srcpot;
                    case 'full CSD'
                        X = procdata.srcpot_all;
                end
                X(~isfinite(X(:))) = 0;
                for k=1:2
                    trials{k} = exp_eval(set_picktrials(procdata,'rank',k));
                    % calculate the CSP class covariance matrix
                    if args.shrinkage_cov
                        CovCSP{k,s} = hlp_diskcache('featuremodels',@cov_shrink,reshape(trials{k}.data,size(trials{k}.data,1),[])');
                    else
                        CovCSP{k,s} = cov(reshape(trials{k}.data,size(trials{k}.data,1),[])');
                    end
                    CovCSP{k,s}(~isfinite(CovCSP{k,s})) = 0;
                end
                % solve CSP
                [V,D] = eig(CovCSP{1,s},CovCSP{1,s}+CovCSP{2,s}); %#ok<NASGU>
                featuremodel.filters = V(:,[1:args.patterns end-args.patterns+1:end]);
                features = self.feature_extract(procdata,featuremodel); % nT x nF
                % get the LDA covariance matrix per class
                %targets = [signal.epoch.target];
                for k=1:2
                    CovLDA{k,s} = cov(features([trials{k}.epoch.targnum],:)); end
            end
            
%             % check for and remove bad data
%             remove = [];
%             for s=1:length(features)
%                 if size(features{s},3) < 2
%                     fprintf('Encountered bad data at subject %s/%i.\n',corpus(s).streams{1}.parts{2}.parts{1:2}); 
%                     remove(end+1) = s;
%                 end
%             end
%             if ~isempty(remove)
%                 fprintf('Removing bad data...\n'); 
%                 scales(remove) = [];
%                 features(remove) = [];
%                 targets(remove) = [];
%                 transforms(remove) = [];
%                 % data_weights(remove) = [];
%             end

        
            % calculate Riemannian weights for each subject
            Covs = vertcat(CovCSP,CovLDA);
            % for each type of covariance matrix...
            for c=size(Covs,1):-1:1
                Ctarg = Covs{c,1};
                
                % calc distance to all other covs
                for s=size(Covs,2):-1:2
                    Cother = Covs{c,s};
                    [V,D] = eig(Ctarg\Cother);
                    dist(s) = sqrt(sum(log(diag(D)).^2));
                end
                
                % calc weighting
                sumdist = sum(dist(2:end));
                for s=2:size(Covs,2)
                    weights(s) = 1./(dist(s)/sumdist); end
                
                % use it to avg the covs
                Avg = zeros(length(Ctarg));
                for s=2:size(Covs,2)
                    Avg = Avg + weights(s)*Covs{c,s}; end
                AvgCovs{c} = Avg;
            end
            
            % for each reg. param, solve a classifier
            signal = procdata;
            targets = [signal.epoch.target];
            classes = unique(targets);
            for li = length(args.lambdas):-1:1
                lam = args.lambdas(li);
                for c=size(Covs,1):-1:1
                    BlendCovs{c} = lam*Covs{c,1} + (1-lam)*AvgCovs{c}; end
                [V,D] = eig(BlendCovs{1},BlendCovs{1}+BlendCovs{2}); P = inv(V); %#ok<NASGU>
                % train CSP part
                model.featuremodel.filters{li} = V(:,[1:args.patterns end-args.patterns+1:end]);
                model.featuremodel.patterns{li} = P([1:args.patterns end-args.patterns+1:end],:);
                % extract features
                trials = self.feature_extract(signal,model.featuremodel,li);
                % train LDA part
                for c = 1:2
                    X = trials(targets==classes(c),:);
                    n{c} = size(X,1);
                    mu{c} = mean(X,1);
                    sig{c} = BlendCovs{c+2};
                end
                ns = quickif(args.weight_cov,n,{1 1});
                nb = quickif(args.weight_bias,n,{1 1});
                % do the math
                mu_both = (mu{1}*nb{2} + mu{2}*nb{1}) / (nb{1}+nb{2});    
                sig_both = (sig{1}*ns{1} + sig{2}*ns{2}) / (ns{1}+ns{2});
                w = (mu{2} - mu{1}) / sig_both;
                w = w / (mu{2}*w' - mu_both*w');
                model.predictivemodel.model{li} = struct('w',{w}, 'b',{mu_both*w'}, 'classes',{classes});
            end
            
            % set the filter graph based on the reference data
            model.tracking.filter_graph = signal;
            % also store channel locations for model visualization
            model.chanlocs = signal.chanlocs;
            model.classes = classes;
        end
        
        function predictions = predict(self,bundle,model)
            % for each lambda
            for m=length(model.featuremodel.filters):-1:1
                % extract features
                features = self.feature_extract(bundle.streams{1},model.featuremodel.filters{m});
                % apply classifier
                raw_preds(:,m) = features*model.predictivemodel.model{m}.w' - model.predictivemodel.model{m}.b;
            end
            % average all predictions
            raw_labels = mean(raw_preds,2);
            raw_labels = min(+1,max(-1,raw_labels));
            predictions = {'disc', [(1-raw_labels)/2 1-(1-raw_labels)/2], model.classes};
        end
        
        function features = feature_extract(self,signal,featuremodel,lam)
            if isstruct(featuremodel)
                featuremodel = featuremodel.filters; end
            if iscell(featuremodel)
                featuremodel = featuremodel{lam}; end
            % extract log-variance features from an epoched and preprocessed recording
            features = zeros(size(signal.data,3),size(featuremodel,2));
            for t=1:size(signal.data,3)
                features(t,:) = sum((signal.data(:,:,t)'*featuremodel).^2,1); end
            features = log(features/size(signal.data,2));
        end
        
        function visualize(self,varargin) %#ok<*INUSD>
            % visualize an mklCSP model
            args = arg_define(varargin, ...
                arg_norep({'model','Model'},[],[],'BCI Model to visualize.'), ...
                arg({'patterns','PlotPatterns'},true,[],'Plot patterns instead of filters. Whether to plot spatial patterns (forward projections) rather than spatial filters.'), ...
                arg({'paper','PaperFigure'},false,[],'Use paper-style font sizes. Whether to generate a plot with font sizes etc. adjusted for paper.'));

            f = figure;            
            % get number of pairs, and index of pattern per subplot
            np = size(args.model.featuremodel.patterns,1)/2; 
            idx = [1:np 2*np:-1:np+1];
            % for each CSP pattern...
            for p=1:np*2
                subplot(2,np,p,'Parent',f);
                if args.patterns
                    topoplot(args.model.featuremodel.patterns(idx(p),:),args.model.featuremodel.chanlocs);
                else
                    topoplot(args.model.featuremodel.filters(:,idx(p)),args.model.featuremodel.chanlocs);
                end
                t = title(['CSP Pattern ' num2str(idx(p))]);
                if args.paper
                    set(t,'FontUnits','normalized');
                    set(t,'FontSize',0.1);                    
                end
            end
        end
        
        function layout = dialog_layout_defaults(self)
            % define the default configuration dialog layout 
            layout = {'SignalProcessing.Resampling.SamplingRate', 'SignalProcessing.FIRFilter.Frequencies', ...
                'SignalProcessing.FIRFilter.Type', 'SignalProcessing.EpochExtraction', '', ...
                'PatternPairs', 'CovarianceShrinkage', '', 'MachineLearning.Learner'};
        end
                
    end
end
            
% (turn off a few editor warnings because some actual implementations are missing in this file)
%#ok<*INUSD,*STOUT,*MANU>
