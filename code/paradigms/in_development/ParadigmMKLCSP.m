classdef ParadigmMKLCSP < ParadigmBase
    % Multiple Kernel Learning Common Spatial Patterns (mklCSP)
    %
    % This paradigm implements mklCSP [1], which is a generalization of the Common Spatial Patterns
    % algorithms to calibration data comprising multiple subjects (or recordings). The basic
    % approach is to solve CSP for each subject separately and then to concatenate all spatial
    % filters into a filter bank. Then, to train a model for a particular subject, the calibration
    % data for that subject is filtered through the filter bank, log-variance features are extracted
    % (which yields Nsubj x Npatterns features), and a classifier is trained which has a
    % group-sparsity penalty with groups of size Npatterns (one group per subject), which is set up
    % to jointly penalize all filters computed for a particular subject.
    %
    % Thereby only features are retained that stem from subjects whose spatial filters are
    % sufficiently applicable to that of the target subject.
    %
    % References:
    % [1] Samek, W., Binder, A., & Muller, K. R. 
    %     "Multiple kernel learning for brain-computer interfacing."
    %      In Engineering in Medicine and Biology Society (EMBC) pp. 7048-7051 (2013)
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
            defaults = {'dal', 2.^(10:-0.25:-1), 'Scaling','none', 'Regularizer','glc'};
        end
                
        function model = calibrate(self,varargin)
            % calibrate an mklCSP model from a corpus of training sets
            args = arg_define(varargin, ...
                arg_norep({'collection','Collection'}), ...
                arg_norep({'goal_identifier','GoalIdentifier'}), ...
                arg({'patterns','PatternPairs'},3,[],'Number of CSP patterns (times two).','cat','Feature Extraction','type','expression','shape','row'),...
                arg({'shrinkage','ShrinkageLevel'},0,[0 1],'Shrinkage level. The amount of shrinkage (regularization) to apply during covariance estimation.'), ...
                arg({'verbose','Verbose'},true,[],'Verbose output.'), ...
                arg_sub({'flt','SignalProcessing'}, self.preprocessing_defaults(), @flt_pipeline, 'Signal processing stages. These parameters control filter stages that run on the signal level; they can be enabled, disabled and configured for the given paradigm. The prediction operates on the outputs of this stage.'), ...
                arg_sub({'ml','MachineLearning'},{'Learner',self.machine_learning_defaults()},@ml_train,'Machine learning stage of the paradigm. Operates on the feature vectors that are produced by the feature-extraction stage.'),...
                arg({'arg_dialogsel','ConfigLayout'},self.dialog_layout_defaults(),[],'Parameters displayed in the config dialog. Cell array of parameter names to display (dot-notation allowed); blanks are translated into empty rows in the dialog. Referring to a structure argument lists all parameters of that struture, except if it is a switchable structure - in this case, a pulldown menu with switch options is displayed.','type','cellstr','shape','row'));
           
            % first solve CSP for each subject in the corpus individually
            if args.verbose
                fprintf('Now training model for: %s...\n',hlp_tostring(args.goal_identifier)); end
            filters = [];
            patterns = [];
            if args.verbose
                fprintf('Pre-processing each of %i sets in the corpus and solving CSP...\n',length(args.collection)); end
            for s=1:length(args.collection)
                if length(args.collection{s}.streams) > 1
                    disp_once('Note: ParadigmMKLCSP will use only the first data stream of a recording (no support for multi-modal data).'); end
                % preprocess
                procdata = exp_eval_optimized(flt_pipeline('signal',args.collection{s}.streams{1}, args.flt)); %#ok<*NODEF>
                if procdata.nbchan < patterns
                    error('CSP requires at least as many channels as you request output patterns. Please reduce the number of pattern pairs.'); end
                % solve CSP
                for k=1:2
                    classdata = exp_eval(set_picktrials(procdata,'rank',k));
                    covar{k} = reshape(classdata.data,size(classdata.data,1),[])*reshape(classdata.data,size(classdata.data,1),[])'/(size(classdata.data,2)*size(classdata.data,3)); % cov(reshape(classdata.data,size(classdata.data,1),[])');
                    covar{k}(~isfinite(covar{k})) = 0;
                    covar{k} = (1-args.shrinkage)*covar{k} + args.shrinkage*eye(size(covar{k}))*trace(covar{k})/length(covar{k});
                end
                [V,D] = eig(covar{1},covar{1}+covar{2}); %#ok<NASGU>
                P = inv(V);
                % if you get an error here then your data sets had varying number of channels
                filters = [filters V(:,[1:args.patterns end-args.patterns+1:end])];
                patterns = [patterns P([1:args.patterns end-args.patterns+1:end],:)'];
            end
            model.featuremodel = struct('filters',filters,'patterns',patterns);

            if args.verbose
                fprintf('Preprocessing and extracting features for reference data...\n'); end            
            % get the data of the reference subject
            [reference,remaining] = utl_collection_closest(args.collection,args.goal_identifier); %#ok<NASGU>
            % preprocess each recording in the reference collection and concatenate them across epochs into a single set
            for r=1:length(reference)
                refsets{r} = exp_eval_optimized(flt_pipeline('signal',reference{r}.streams{1}, args.flt)); end
            refdata = exp_eval(set_joinepos(refsets{:}));
            % extract features and get target labels
            features = self.feature_extract(refdata,model.featuremodel);
            targets = set_gettarget(refdata);
            if args.verbose
                fprintf('Training predictive model...\n'); end
            % train classifier, using the correct feature shape (based on the group size)
            args.ml.learner.shape = [2*args.patterns,length(args.collection)];
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
                'Prediction.FeatureExtraction.PatternPairs', '', 'Prediction.MachineLearning.Learner'};
        end
                
    end
end
            
% (turn off a few editor warnings because some actual implementations are missing in this file)
%#ok<*INUSD,*STOUT,*MANU>
