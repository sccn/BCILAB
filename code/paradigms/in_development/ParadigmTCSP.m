classdef ParadigmTCSP < ParadigmBase
    % Transfer Common Spatial Patterns (TCSP)
    % 
    % This class is currently a testbed for experimental transfer-learning enabled methods in the
    % CSP family.
    %
    % Name:
    %   Multiple Kernel Learning Common Spatial Patterns
    %
    %                            Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
    %                            2014-02-05
    
    methods
        
        function defaults = preprocessing_defaults(self)
            % most basic version
            defaults = {'FIRFilter',{'Frequencies',[6 8 28 32],'Type','minimum-phase'}, 'EpochExtraction',[0.5 3.5], 'Resampling',100};
        end
                
        function defaults = machine_learning_defaults(self)
            % set up the default parameters for machine learning
            defaults = {'logreg', 'variant','vb'};
        end
                
        function model = calibrate(self,varargin)
            % calibrate an mklCSP model from a corpus of training sets
            args = arg_define(varargin, ...
                arg_norep({'collection','Collection'}), ...
                arg_norep({'goal_identifier','GoalIdentifier'}), ...
                arg({'patterns','PatternPairs'},3,uint32([1 1 64 10000]),'Number of CSP patterns (times two).','cat','Feature Extraction','type','expression','shape','row'),...
                arg({'shrinkage','ShrinkageLevel'},0,[0 1],'Shrinkage level. The amount of shrinkage (regularization) to apply during covariance estimation.'), ...
                arg({'verbose','Verbose'},true,[],'Verbose output.'), ...
                arg_sub({'flt','SignalProcessing'}, self.preprocessing_defaults(), @flt_pipeline, 'Signal processing stages. These parameters control filter stages that run on the signal level; they can be enabled, disabled and configured for the given paradigm. The prediction operates on the outputs of this stage.'), ...
                arg_sub({'ml','MachineLearning'},{'Learner',self.machine_learning_defaults()},@ml_train,'Machine learning stage of the paradigm. Operates on the feature vectors that are produced by the feature-extraction stage.'),...
                arg({'arg_dialogsel','ConfigLayout'},self.dialog_layout_defaults(),[],'Parameters displayed in the config dialog. Cell array of parameter names to display (dot-notation allowed); blanks are translated into empty rows in the dialog. Referring to a structure argument lists all parameters of that struture, except if it is a switchable structure - in this case, a pulldown menu with switch options is displayed.','type','cellstr','shape','row'));

            if args.verbose
                fprintf('Now training model for: %s...\n',hlp_tostring(args.goal_identifier)); end

            % first solve CSP for each subject in the corpus individually and aggregate CSP filters            
            
            % find the unique subjects in the collection
            try
                corpus = [args.collection{:}];
            catch e
                error('The dataset collection must have the same field names for each recording.');
            end
            if ~isfield(corpus,'subject')
                error('The datasets in the collection must each have a .subject field.'); end
            subjects = {corpus.subject};
            if all(cellfun('isclass',subjects,'char'))
                subjects = unique(subjects);
            elseif all(cellfun('isclass',subjects,'double'))
                subjects = unique([subjects{:}]);
            else
                error('The subject identifiers must either be all strings or all doubles');
            end
            
            if args.verbose
                fprintf('Pre-processing each of %i recordings (%i subjects) in the corpus and solving CSP...\n',length(args.collection),length(subjects)); end

            % for each subject in the collection...
            [preproc,weighted_cov,mean_cov] = deal(cell(1,length(subjects)));
            for s=subjects
                % find all recordings that match that subject
                recordings = args.collection(cellfun(@(x)isequal(x,s),{corpus.subject}));
                recordings = cellfun(@(x)x.streams(1),recordings);
                % concatenate and preprocess the data
                concat = set_concat(recordings{:});
                preproc{s} = exp_eval_optimized(flt_pipeline('signal',concat,args.flt));
                % extract class-weighted and average covariance with shrinkage
                weighted_cov{s} = zeros(preproc{s}.nbchan);
                mean_cov{s} = zeros(preproc{s}.nbchan);
                try
                    for k=1:preproc{s}.trials
                        weighted_cov{s} = weighted_cov{s} + ParadigmTCSP.cov_shrinkage(preproc{s}.data(:,:,k),args.shrinkage) * preproc{s}.epoch(k).target; 
                        mean_cov{s} = mean_cov{s} + ParadigmTCSP.cov_shrinkage(preproc{s}.data(:,:,k),args.shrinkage);
                    end
                catch e
                    fprintf('WARNING: could not compute COVs for subject %s: %s',hlp_tostring(s),e.message);
                end
            end
            
            % aggregate the covariance matrices
            weighted_cov = sum(cat(3,weighted_cov{:}),3)/length(weighted_cov);
            mean_cov = sum(cat(3,mean_cov{:}),3)/length(mean_cov);
            
            % solve CSP/SPoC
            [V,D] = eig(weighted_cov,mean_cov); %#ok<ASGLU,NASGU>
            P = inv(V);
            filters = V(:,[1:args.patterns end-args.patterns+1:end]);
            patterns = P([1:args.patterns end-args.patterns+1:end],:);
            
            model.featuremodel = struct('filters',filters,'patterns',patterns);

            % extract features and get target labels for all subjects
            if args.verbose
                fprintf('Training predictive model (this may take a while)...\n'); end
            preproc_concat = exp_eval(set_joinepos(preproc{:}));
            features = self.feature_extract(preproc_concat,model.featuremodel);
            targets = set_gettarget(preproc_concat);
            
            % train predictive model
            if args.verbose
                fprintf('Training predictive model (this may take a while)...\n'); end
            model.predictivemodel = ml_train('data',{features,targets}, args.ml);
            
            % adapt the filter graph, based on the available reference data for the subject
            [reference,remaining] = utl_collection_closest(args.collection,args.goal_identifier); %#ok<ASGLU,NASGU>
            reference = cellfun(@(r)r.streams(1),reference);
            model.tracking.filter_graph = exp_eval_optimized(flt_pipeline('signal',set_concat(reference{:}),args.flt));           
            %model.tracking.filter_graph = exp_eval_optimized(flt_pipeline('signal',args.collection{1}.streams{1},args.flt));
            
            % also store channel locations for model visualization
            model.chanlocs = model.tracking.filter_graph.chanlocs;
        end
        
        function predictions = predict(self,bundle,model)
            % extract features
            features = self.feature_extract(bundle.streams{1},model.featuremodel);
            % apply classifier
            predictions = ml_predict(features, model.predictivemodel);
        end
        
        function features = feature_extract(self,signal,featuremodel)
            % extract log-variance features from an epoched and preprocessed recording
            features = zeros(size(signal.data,3),size(featuremodel.filters,2));
            for t=1:size(signal.data,3)
                features(t,:) = sum((signal.data(:,:,t)'*featuremodel.filters).^2,1); end
            features = log(features/size(signal.data,2));
        end
        
        function visualize(self,varargin) %#ok<*INUSD>
            % visualize an mklCSP model
            args = arg_define(varargin, ...
                arg_norep({'model','Model'},[],[],'BCI Model to visualize.'), ...
                arg({'patterns','PlotPatterns'},true,[],'Plot patterns instead of filters. Whether to plot spatial patterns (forward projections) rather than spatial filters.'), ...
                arg({'paper','PaperFigure'},false,[],'Use paper-style font sizes. Whether to generate a plot with font sizes etc. adjusted for paper.'));
            
            f = figure;            

            % mask out unused filters
            mask = args.model.predictivemodel.model.w(:)' ~= 0;
            args.model.featuremodel.patterns = args.model.featuremodel.patterns(:,mask);
            args.model.featuremodel.filters = args.model.featuremodel.filters(:,mask);
            
            % number of plots, and index of pattern per subplot            
            np = nnz(mask);
            horz = ceil(sqrt(np));
            vert = ceil(np/horz);
                
            % get number of pairs, and index of pattern per subplot
            % for each CSP pattern...
            for p=1:np
                subplot(horz,vert,p,'Parent',f);
                if args.patterns
                    topoplot(args.model.featuremodel.patterns(:,p),args.model.chanlocs);
                else
                    topoplot(args.model.featuremodel.filters(:,p),args.model.chanlocs);
                end
                t = title(['CSP Pattern ' num2str(p)]);
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
                'PatternPairs', '', 'MachineLearning.Learner'};
        end
                
    end
    methods (Static)
        % estimate covariance with shrinkage for some data X
        function V = cov_shrinkage(X,shrinkage)
            V = reshape(X,size(X,1),[])*reshape(X,size(X,1),[])'/(size(X,2)*size(X,3));
            V(~isfinite(V)) = 0;
            V = (1-shrinkage)*V + shrinkage*eye(size(V))*trace(V)/length(V);
        end        
    end
    
end
            
% (turn off a few editor warnings because some actual implementations are missing in this file)
%#ok<*INUSD,*STOUT,*MANU>
