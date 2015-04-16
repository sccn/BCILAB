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
            defaults = {'FIRFilter',{'Frequencies',[6 8 28 32],'Type','minimum-phase'}, 'EpochExtraction',[0.5 3.5], 'Resampling',100};
        end
                
        function defaults = machine_learning_defaults(self)
            % set up the default parameters for machine learning
            defaults = {'dal', 2.^(4:-0.25:-3), 'Scaling','none', 'Regularizer','grouplasso-columns'};
            % defaults = {'logreg', 'variant','lars'};
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
            filters = [];
            patterns = [];
            
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

            % remove actual data from corpus so we can micro-cache it
            for s=1:length(corpus)
                if isfield(corpus(s).streams{1},'tracking')
                    corpus(s).streams{1} = corpus(s).streams{1}.tracking.expression; end
            end
            
            % for each subject...
            for subj=subjects                
                % find all recordings that match that subject
                recordings = corpus(cellfun(@(s)isequal(s,subj),{corpus.subject}));
                % calculate CSP filters
                [newfilters,newpatterns] = hlp_microcache('filters',@ParadigmMKLCSP.filters_for_subject,recordings, args.flt, args.shrinkage, args.patterns);
                % if you get an error here then your data sets had varying number of channels
                filters = [filters newfilters];
                patterns = [patterns newpatterns];
            end
            model.featuremodel = struct('filters',filters,'patterns',patterns);

            if args.verbose
                fprintf('Preprocessing and extracting features for reference data...\n'); end            
            % get the data of the reference subject
            [reference,remaining] = utl_collection_closest(args.collection,args.goal_identifier); %#ok<ASGLU,NASGU>
            % preprocess each recording in the reference collection and concatenate them across epochs into a single set
            for r=1:length(reference)
                refsets{r} = exp_eval_optimized(flt_pipeline('signal',reference{r}.streams{1}, args.flt)); end
            refdata = exp_eval(set_joinepos(refsets{:}));
            % extract features and get target labels
            features = self.feature_extract(refdata,model.featuremodel);
            targets = set_gettarget(refdata);
            
            if args.verbose
                fprintf('Training predictive model (this may take a while)...\n'); end
            % train classifier, overriding with the correct feature shape (based on the group size)
            if isfield(args.ml.learner,'shape')
                args.ml.learner.shape = [2*args.patterns,length(subjects)]; end
            model.predictivemodel = ml_train('data',{features,targets}, args.ml);
            % set the filter graph based on the reference data
            model.tracking.filter_graph = refsets{end};
            % also store channel locations for model visualization
            model.chanlocs = refdata.chanlocs;
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
        function [filters, patterns] = filters_for_subject(recordings, flt, shrinkage, n_patterns)
            % get the CSP filters for a given subject
            
            % preprocess and concatenate across trials
            preproc = {};
            for r=1:length(recordings)
                preproc{r} = exp_eval_optimized(flt_pipeline('signal',recordings(r).streams{1}, flt)); end
            preproc = exp_eval(set_joinepos(preproc{:}));

            for k=1:2
                % filter trial subrange in time and frequency
                classdata = exp_eval(set_picktrials(preproc,'rank',k));
                covar{k} = reshape(classdata.data,size(classdata.data,1),[])*reshape(classdata.data,size(classdata.data,1),[])'/(size(classdata.data,2)*size(classdata.data,3)); % cov(reshape(classdata.data,size(classdata.data,1),[])');
                covar{k}(~isfinite(covar{k})) = 0; %#ok<*AGROW>
                covar{k} = (1-shrinkage)*covar{k} + shrinkage*eye(size(covar{k}))*trace(covar{k})/length(covar{k});
            end            
            try
                [V,D] = eig(covar{1},covar{1}+covar{2}); %#ok<ASGLU,NASGU>
                P = inv(V);
                % if you get an error here then your data sets had varying number of channels
                filters = V(:,[1:n_patterns end-n_patterns+1:end]);
                patterns = P([1:n_patterns end-n_patterns+1:end],:)';
            catch e
                fprintf('Got a degenerate CSP solution, replacing by identity matrix:%s\n',e.message);
                n_chans = preproc.nbchan;
                if ~n_chans
                    % no epochs, need to determine the number of channels in the filter stage prior to epoching
                    raw = utl_get_argument(utl_find_filter(preproc,'set_makepos'),'signal');
                    n_chans = raw.nbchan;
                end
                filters = eye(n_chans,2*n_patterns);
                patterns = eye(n_chans,2*n_patterns);
            end
        end        
    end
    
end
            
% (turn off a few editor warnings because some actual implementations are missing in this file)
%#ok<*INUSD,*STOUT,*MANU>
