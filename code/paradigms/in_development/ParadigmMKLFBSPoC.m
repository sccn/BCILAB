classdef ParadigmMKLFBSPoC < ParadigmBase
    % Multiple Kernel Learning Filter-Bank Source Power Comodulation (mklFBSPoC)
    %
    % This paradigm implements a generalization of mklCSP [1] to multiple frequency bands as in 
    % FBCSP [2] and to multiple time windows, and using the SPoC formulation for regression.
    %
    % References:
    % [1] Samek, W., Binder, A., & Muller, K. R. 
    %     "Multiple kernel learning for brain-computer interfacing."
    %     In Engineering in Medicine and Biology Society (EMBC) pp. 7048-7051 (2013)
    % [2] Kai K. Ang, Zhang Y. Chin, Haihong Zhang, Cuntai Guan, 
    %     "Filter Bank Common Spatial Pattern (FBCSP) in Brain-Computer Interface"
    %     In 2008 IEEE International Joint Conference on Neural Networks (IEEE World Congress on Computational Intelligence) (June 2008), pp. 2390-2397.
    % [3] Daehne, S., Meinecke, F. C., Haufe, S., Hoehne, J., Tangermann, M., Mueller, K. R., & Nikulin, V. V.
    %     "SPoC: A novel framework for relating the amplitude of neuronal oscillations to behaviorally relevant parameters."
    %     NeuroImage 86 (2014), 111-122.
    %
    % Name:
    %   Multiple Kernel Learning Filter-Bank Source Power Comodulation
    %
    %                            Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
    %                            2015-04-15
    
    methods
        
        function defaults = preprocessing_defaults(self)
            defaults = {'EpochExtraction',[0.5 3.5],'Resampling',200};
        end
        
        function defaults = machine_learning_defaults(self)
            % set up the default parameters for machine learning
            defaults = 'ridge';
        end
                
        function model = calibrate(self,varargin)
            % calibrate an mklCSP model from a corpus of training sets
            args = arg_define(varargin, ...
                arg_norep({'collection','Collection'}), ...
                arg_norep({'goal_identifier','GoalIdentifier'}), ...
                arg({'patterns','PatternPairs'},3,uint32([1 1 64 10000]),'Number of SPoC patterns (times two).','cat','Feature Extraction','type','expression','shape','row'),...
                arg({'shrinkage','ShrinkageLevel'},0,[0 1],'Shrinkage level. The amount of shrinkage (regularization) to apply during covariance estimation.'), ...
                arg({'freqwnds','FreqWindows'},[0.5 3; 4 7; 8 12; 13 30; 31 42],[0 0.5 200 1000],'Frequency bands of interest. Matrix containing one row for the start and end of each frequency band from which CSP patterns shall be computed. Values in Hz.','cat','Feature Extraction'), ...
                arg({'timewnds','TimeWindows'},[],[],'Time windows of interest. Matrix containing one row for the start and end of each time window from which CSP patterns shall be computed. Values in seconds. If both this and the freqwnds parameter are non-empty, they should have the same number of rows.','cat','Feature Extraction'), ...
                arg({'winfunc','WindowFunction'},'rect',{'barthann','bartlett','blackman','blackmanharris','bohman','cheb','flattop','gauss','hamming','hann','kaiser','nuttall','parzen','rect','taylor','triang','tukey'},'Type of window function. Typical choices are rect (rectangular), hann, gauss, blackman and kaiser.'),...
                arg({'winparam','WindowParameter','param'},[],[],'Parameter of the window function. This is mandatory for cheb, kaiser and tukey and optional for some others.','shape','scalar'), ...
                arg({'verbose','Verbose'},true,[],'Verbose output.'), ...
                arg({'hotpatching','HotPatching'},false,[],'Hot-patch the data. This can be enabled to ensure that a long-running computation survives bad data.'), ...
                arg_sub({'flt','SignalProcessing'}, self.preprocessing_defaults(), @flt_pipeline, 'Signal processing stages. These parameters control filter stages that run on the signal level; they can be enabled, disabled and configured for the given paradigm. The prediction operates on the outputs of this stage.'), ...
                arg_sub({'ml','MachineLearning'},{'Learner',self.machine_learning_defaults()},@ml_train,'Machine learning stage of the paradigm. Operates on the feature vectors that are produced by the feature-extraction stage.'),...
                arg({'arg_dialogsel','ConfigLayout'},self.dialog_layout_defaults(),[],'Parameters displayed in the config dialog. Cell array of parameter names to display (dot-notation allowed); blanks are translated into empty rows in the dialog. Referring to a structure argument lists all parameters of that struture, except if it is a switchable structure - in this case, a pulldown menu with switch options is displayed.','type','cellstr','shape','row'));

            if ~isempty(args.freqwnds) && ~isempty(args.timewnds) && size(args.freqwnds,1) ~= size(args.timewnds,1)
                error('If both time and frequency windows are specified, both arrays must have the same number of rows (together they define the windows in time and frequency).'); end
            if isempty(args.timewnds)
                args.timewnds = zeros(size(args.freqwnds,1),0); end
            if isempty(args.freqwnds)
                args.freqwnds = zeros(size(args.timewnds,1),0); end
            
            % pre-parse arguments for flt_window and flt_spectrum (for fast subsequent online use)
            for w = 1:max(size(args.freqwnds,1),size(args.timewnds,1))                
                time_args{w} = arg_report('vals',@flt_window,{'time',{args.timewnds(w,:),args.winfunc,args.winparam}});
                freq_args{w} = arg_report('vals',@flt_spectrum,{'freq',args.freqwnds(w,:)});
            end
            
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
                % calculate FBSPoC filters
                [newfilters,newpatterns,chanlocs] = hlp_microcache('filters',@ParadigmMKLFBSPoC.filters_for_subject,recordings, args.flt, time_args, freq_args, args.shrinkage, args.patterns);
                % if you get an error here then your data sets had varying number of channels
                filters = [filters newfilters];
                patterns = [patterns newpatterns];
            end
            model.featuremodel = struct('filters',{filters},'patterns',{patterns}, ...
                'n_subjects',length(subjects),'time_args',{time_args},'freq_args',{freq_args}, ...
                'chanlocs',{chanlocs}, 'hotpatching', {args.hotpatching});
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
     
            if args.hotpatching && size(features,1) < 5
                fprintf('You have too few trials in the data; hot-patching.\n');
                features = features(1+mod(0:4,size(features,1)),:);
                features = features+0.1*randn(size(features));
                targets = 1+mod(0:size(features,1)-1,2);
                targets = targets(:);
            end
                            
            if args.hotpatching && length(targets) ~= size(features,1)
                fprintf('Your # of target markers does not match the # of extracted features; hot-patching.\n');
                if isempty(targets)
                    targets = 1+mod(0:size(features,1)-1,2);
                else
                    targets = targets(1+mod(0:size(features,1)-1,length(targets)));
                end
                targets = targets(:);
            end
            
            if args.hotpatching && length(unique(targets))==1
                fprintf('Your reference data has only one class; hot-patching the data.\n');
                for ii=1:min(length(targets),max(2,round(length(targets)/10)))
                    targets(ii) = 3-targets(ii); end
            end
                                   
            if args.hotpatching && any(~isfinite(features(:)))
                fprintf('Some of your features are non-finite; hot-patching the data.\n');
                tofix = find(~isfinite(features(:)));
                features(tofix) = randn(1,length(tofix));
            end
            
            if args.verbose
                fprintf('Training predictive model (this may take a while)...\n'); end
            % train classifier, overriding with the correct feature shape (based on the group size)
            if isfield(args.ml.learner,'shape')
                args.ml.learner.shape = [2*args.patterns,length(subjects)]; end
            try
                model.predictivemodel = ml_train('data',{features,targets}, args.ml);
            catch e
                if args.hotpatching && ~isempty(strfind(e.message,'Null probability for class'))
                    fprintf('One of the classes has a probability of 0; hot-patching the data.\n');
                    targets = 1+mod(0:length(targets)-1,2); targets = targets(:);
                    model.predictivemodel = ml_train('data',{features,targets}, args.ml);
                else
                    rethrow(e);
                end
            end
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
            try
                N = featuremodel.n_subjects;
                W = length(featuremodel.freq_args);
                F = size(featuremodel.filters,2)/N;
                T = size(signal.data,3);
                features = zeros(T,F*N);
                for w = 1:W
                    % filter data in time & frequency
                    data = exp_eval(flt_spectrum('signal',flt_window('signal',signal,featuremodel.time_args{w}),featuremodel.freq_args{w}));
                    mask = false(1,F); mask((w-1)*(F/W)+(1:(F/W))) = true; mask = repmat(mask,1,N);
                    for t=1:size(signal.data,3)
                        features(t,mask) = sum((data.data(:,:,t)'*featuremodel.filters(:,mask)).^2,1); end
                end
                features = log(features/size(signal.data,2));
            catch e
                if featuremodel.hotpatching
                    fprintf('Trying to prevent error during feature extraction: %s\n',hlp_handleerror(e));
                    fprintf('size(featuremodel.filters): %s\n',hlp_tostring(size(featuremodel.filters)));
                    fprintf('size(signal.data): %s\n',hlp_tostring(size(signal.data)));
                    fprintf('trying to hot-fix the issue...');
                    featuremodel.filters = featuremodel.filters(1+mod(0:size(signal.data,1)-1,size(featuremodel.filters,1)),:);
                    features = self.feature_extract(signal,featuremodel);
                    fprintf('succeeded.\n');
                else
                    rethrow(e);
                end
            end                
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
        function [filters, patterns, chanlocs] = filters_for_subject(recordings, flt, time_args, freq_args, shrinkage, n_patterns)
            % get the CSP filters for a given subject
            
            % preprocess and concatenate across trials
            preproc = {};
            for r=1:length(recordings)
                preproc{r} = exp_eval_optimized(flt_pipeline('signal',recordings(r).streams{1}, flt)); end
            preproc = exp_eval(set_joinepos(preproc{:}));

            % standardize target values
            targets = [preproc.epoch.target];
            targets = (targets - mean(targets)) / std(targets);

            % for each window...
            filters = [];
            patterns = [];            
            for w = 1:length(time_args)
                signal = exp_eval(flt_spectrum('signal',flt_window('signal',preproc,time_args{w}),freq_args{w}));
                weighted_cov = zeros(signal.nbchan);
                mean_cov = zeros(signal.nbchan);
                for k=signal.trials:-1:1
                    trialcov = reshape(signal.data(:,:,k)',size(signal.data,1),[])*reshape(signal.data(:,:,k)',size(signal.data,1),[])'/size(signal.data,2);
                    trialcov(~isfinite(trialcov)) = 0; 
                    trialcov = (1-shrinkage)*trialcov + shrinkage*eye(size(trialcov))*trace(trialcov)/length(trialcov);
                    weighted_cov = weighted_cov + trialcov * targets(k);
                    mean_cov = mean_cov + trialcov;
                end
                try
                    [V,D] = eig(weighted_cov,mean_cov); %#ok<NASGU>
                    P = inv(V);
                    % if you get an error here then your data sets had varying number of channels                
                    filters = [filters V(:,[1:n_patterns end-n_patterns+1:end])];
                    patterns = [patterns P([1:n_patterns end-n_patterns+1:end],:)'];                    
                catch e
                    fprintf('Got a degenerate FBSPoC solution, replacing by identity matrix:%s\n',e.message);
                    n_chans = preproc.nbchan;
                    if ~n_chans
                        % no epochs, need to determine the number of channels in the filter stage prior to epoching
                        raw = exp_eval_optimized(utl_get_argument(utl_find_filter(preproc,'set_makepos'),'signal'));
                        n_chans = raw.nbchan;
                    end
                    filters = [filters eye(n_chans,2*n_patterns)];
                    patterns = [patterns eye(n_chans,2*n_patterns)];
                end                
            end
            chanlocs = preproc.chanlocs;
        end        
    end
end
            
% (turn off a few editor warnings because some actual implementations are missing in this file)
%#ok<*INUSD,*STOUT,*MANU>
