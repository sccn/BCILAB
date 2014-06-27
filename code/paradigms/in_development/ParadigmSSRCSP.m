classdef ParadigmSSRCSP < ParadigmBase
    % Selected-Subjects Regularized CSP (SSRCSP)
    %
    % This paradigm implements SSRCSP [1], which is a generalization of the Common Spatial Patterns
    % algorithms to calibration data comprising multiple subjects (or recordings). To train a model
    % for a particular "goal" (or target) subject using auxiliary data from other subjects, this
    % algorithm attempts to find a subset of other subjects such that, when data from those subjects
    % is combined with data from the goal subject, the performance on the goal subject's data is
    % optimal. The combination of data from multiple subjects and from the goal subject is done by
    % shrinking the goal-subject covariance matrix towards the average covariance matrix of the
    % other subjects (within each class), using a regularization parameter. The subset selection
    % algorithm employed is Sequential Floating Forward Selection (SFFS) [2].
    % 
    % References:
    % [1] Lotte, F., & Guan, C. 
    %     "Regularizing common spatial patterns to improve BCI designs: unified theory and new algorithms."
    %     Biomedical Engineering, IEEE Transactions on, 58(2), 355-362, 2011.
    %
    % [2] Pudil, P., Ferri, F. J., Novovicova, J., & Kittler, J. 
    %     "Floating search methods for feature selection with nonmonotonic criterion functions."
    %     In Pattern Recognition, Proceedings of the 12th IAPR International. Conference on (Vol. 2, pp. 279-283). (1994)
    %
    % Name:
    %   Selected-Subjects Regularized Common Spatial Patterns
    %
    %                            Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
    %                            2014-02-06
    
    methods
        
        function defaults = preprocessing_defaults(self)
            defaults = {'FIRFilter',{'Frequencies',[6 8 28 32],'Type','minimum-phase'}, 'EpochExtraction',[0.5 3.5], 'Resampling',100};
        end
                
        function defaults = machine_learning_defaults(self)
            % set up the default parameters for machine learning
            defaults = {'lda',0.01,'regularization','shrinkage'};
        end
                
        function model = calibrate(self,varargin)
            % calibrate an SSRCSP model from a corpus of training sets
            args = arg_define(varargin, ...
                arg_norep({'collection','Collection'}), ...
                arg_norep({'goal_identifier','GoalIdentifier'}), ...
                arg({'patterns','PatternPairs'},3,uint32([1 1 64 10000]),'Number of CSP patterns (times two).'),...
                arg({'beta','CovarianceShrinkage'},0.5,[0 1],'Covariance shrinkage. This is the degree to which data from the goal subject is shrunken towards that of other subjects.'),...
                arg_sub({'flt','SignalProcessing'}, self.preprocessing_defaults(), @flt_pipeline, 'Signal processing stages. These parameters control filter stages that run on the signal level; they can be enabled, disabled and configured for the given paradigm. The prediction operates on the outputs of this stage.'), ...
                arg_sub({'ml','MachineLearning'},{'Learner',self.machine_learning_defaults()},@ml_train,'Machine learning stage of the paradigm. Operates on the feature vectors that are produced by the feature-extraction stage.'),...
                arg({'arg_dialogsel','ConfigLayout'},self.dialog_layout_defaults(),[],'Parameters displayed in the config dialog. Cell array of parameter names to display (dot-notation allowed); blanks are translated into empty rows in the dialog. Referring to a structure argument lists all parameters of that struture, except if it is a switchable structure - in this case, a pulldown menu with switch options is displayed.','type','cellstr','shape','row'));
           
            % get the data of the reference subject
            [reference,remaining] = utl_collection_closest(args.collection,args.goal_identifier); 
            % preprocess each recording in the reference collection and concatenate them across epochs into a single set
            for r=1:length(reference)
                refsets{r} = exp_eval_optimized(flt_pipeline('signal',reference{r}.streams{1}, args.flt)); end
            refdata = exp_eval(set_joinepos(refsets{:}));

            % pre-process data of all other subjects
            otherdata = {};
            for s=1:length(remaining)
                if length(remaining{s}.streams) > 1
                    disp_once('Note: ParadigmMKLCSP will use only the first data stream of a recording (no support for multi-modal data).'); end
                % preprocess
                otherdata{s} = exp_eval_optimized(flt_pipeline('signal',remaining{s}.streams{1}, args.flt)); %#ok<*NODEF>
                if otherdata{s}.nbchan < args.patterns
                    error('CSP requires at least as many channels as you request output patterns. Please reduce the number of pattern pairs.'); end
            end
            
            % get the best subjects
            model.best_subjects = hlp_diskcache('featuremodels',@self.find_best_subjects,otherdata,refdata,args.patterns,args.ml);
                                    
            % calculate composite CSP
            covar = self.class_covariances(refdata);
            other_covar = self.class_covariances(exp_eval(set_joinepos(otherdata{model.best_subjects})));
            for k=1:2
                covar{k} = (1-args.beta)*covar{k} + args.beta*other_covar{k}; end
            [V,D] = eig(covar{1},covar{1}+covar{2}); P = inv(V); %#ok<NASGU>
            model.featuremodel.filters = V(:,[1:args.patterns end-args.patterns+1:end]);
            model.featuremodel.patterns = P([1:args.patterns end-args.patterns+1:end],:);
            
            % train predictive model
            model.predictivemodel = ml_train('data',{self.feature_extract(refdata,model.featuremodel),set_gettarget(refdata)}, args.ml);
            % set the filter graph based on the reference data
            model.tracking.filter_graph = refsets{end};
            % also store channel locations for model visualization
            model.chanlocs = refdata.chanlocs;
        end
        
        function predictions = predict(self,bundle,model)
            % extract features
            features = self.feature_extract(bundle.streams{1},model.featuremodel);
            % apply classifier
            predictions = ml_predict(features,model.predictivemodel);
        end
        
        function best_subjects = find_best_subjects(self,otherdata,refdata,patterns,ml)
            % find set of best subjects to include
            selected = {[]};
            remaining = {1:length(otherdata)};
            accuracy = {-Inf};
            n = 0;
            while n < length(otherdata)
                % find best subject to add
                best_accuracy = -Inf;
                best_index = NaN;
                for k = remaining{1+n}
                    acc = self.evaluate_subset(otherdata([selected{1+n} k]),refdata,patterns,ml);
                    if acc > best_accuracy
                        best_accuracy = acc;
                        best_index = k;
                    end
                end
                selected{1+n+1} = [selected{1+n} best_index];
                remaining{1+n+1} = setdiff(remaining{1+n},best_index);
                accuracy{1+n+1} = best_accuracy;
                n = n+1;
                % remove subjects
                while n > 2
                    % find best subject to remove
                    best_accuracy = -Inf;
                    best_index = NaN;
                    for k=selected{1+n}
                        acc = self.evaluate_subset(otherdata(setdiff(selected{1+n},k)),refdata,patterns,ml);
                        if acc > best_accuracy
                            best_accuracy = acc;
                            best_index = k;
                        end
                    end
                    if best_accuracy > accuracy{1+n-1}
                        selected{1+n-1} = setdiff(selected{1+n},best_index);
                        remaining{1+n-1} = [remaining{1+n} best_index];
                        accuracy{1+n-1} = best_accuracy;
                        n = n-1;
                    else
                        break;
                    end
                end
            end
            best_n = argmax([accuracy{:}]);
            best_subjects = selected{best_n}; 
        end
        
        function accuracy = evaluate_subset(self,trainset,testset,patterns,ml)
            % note: we are here implicitly weighting by the amount of training data per subject
            if iscell(trainset)
                trainset = exp_eval(set_joinepos(trainset{:})); end
            % train CSP on training set
            covar = self.class_covariances(trainset);
            [V,D] = eig(covar{1},covar{1}+covar{2}); %#ok<NASGU>
            featuremodel.filters = V(:,[1:patterns end-patterns+1:end]);
            % train classifier
            classifier = ml_train('data',{self.feature_extract(trainset,featuremodel),set_gettarget(trainset)},ml);
            % test on test set
            accuracy = -ml_calcloss('auto',set_gettarget(testset),ml_predict(self.feature_extract(testset,featuremodel),classifier));
        end
        
        function covar = class_covariances(self,dataset)
            % calculate per-class covariance matrices
            for k=1:2
                classdata = exp_eval(set_picktrials(dataset,'rank',k));
                covar{k} = (classdata.data(:,:) * classdata.data(:,:)') / (size(classdata.data,2)*size(classdata.data,3));
                covar{k}(~isfinite(covar{k})) = 0;
            end        
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
