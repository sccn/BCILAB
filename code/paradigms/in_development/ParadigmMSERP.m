classdef ParadigmMSERP < ParadigmBase
    % Multi-Subject Event Related-Potentials (MSERP) Method.
    %
    % This paradigm implements MSERP, which is an experimental method to learn event-related
    % potential based BCIs from a group of subjects.
    %
    % Name:
    %   Multi-Subject Event-Related Potentials
    %
    %                            Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
    %                            2014-02-05
    
    methods
        
        function defaults = preprocessing_defaults(self)
            defaults = {'IIRFilter',{[0.05 0.3],'highpass'}, 'FIRFilter',{[14 16],'lowpass'},'EpochExtraction',[-1.5 1.5],'Resampling',60};
        end
                
        function defaults = machine_learning_defaults(self)
            % set up the default parameters for machine learning
            defaults = {'proximal' ...
                'Regularizers', { ...
                    'Term1', 'trace' ...
                    'Term2', {'trace' ...
                        'LinearOperator', '@(x)reshape(x,a*b,[])'}} ...
                'LambdaSearch', { ...
                    'ReturnRegpath', false}};
        end
                
        function model = calibrate(self,varargin)
            % calibrate an MSERP model from a corpus of training sets
            args = arg_define(varargin, ...
                arg_norep({'collection','Collection'}), ...
                arg_norep({'goal_identifier','GoalIdentifier'}), ...
                arg({'spatial_whitening','SpatialWhitening'},1,[0 1],'Degree of spatial whitening. This is a regularization parameter that governs to what extent the data of each subject shall be whitened spatially.'), ...
                arg({'temporal_whitening','TemporalWhitening'},1,[0 1],'Degree of temporal whitening. This is a regularization parameter that governs to what extent the data of each subject shall be whitened temporally.'), ...
                arg({'reference_weight','ReferenceWeight'},0,[0 Inf],'Weight of the reference set. This is the weight that the reference data (of the goal subject) has, while the data of remaining subjects has 1 minus this weight. If set to 0, the reference set is weighted according to the proportion in the corpus. If this is greater than 1, it is taken as a multiplicative factor on top of its proportion in the corpus (e.g., 3 would weigh as much as 3 other subjects).'), ...
                arg({'cov_type','CovarianceType'},'full',{'diag','full','shrink'},'Covariance estimator. The covariance estimator to use; can be diagonal, full covariance, or shrinkage covariance.'), ...
                arg({'normalizers','NormalizationExponents'},[-0.25,-0.25],[],'Normalization exponents [lhs, rhs]. Two-element array of powers for the left-hand-side and right-hand-side normalization matrices that are applied to the data from the region.','guru',true), ...
                arg({'apply_to','ApplyTo'},'channels',{'channels','sources','components','full CSD'},'Apply classifier to. Allows to select the type of time series to apply this model to.'), ...
                arg_sub({'flt','SignalProcessing'}, self.preprocessing_defaults(), @flt_pipeline, 'Signal processing stages. These parameters control filter stages that run on the signal level; they can be enabled, disabled and configured for the given paradigm. The prediction operates on the outputs of this stage.'), ...
                arg_sub({'ml','MachineLearning'},{'Learner',self.machine_learning_defaults()},@ml_train,'Machine learning stage of the paradigm. Operates on the feature vectors that are produced by the feature-extraction stage.'),...
                arg({'arg_dialogsel','ConfigLayout'},self.dialog_layout_defaults(),[],'Parameters displayed in the config dialog. Cell array of parameter names to display (dot-notation allowed); blanks are translated into empty rows in the dialog. Referring to a structure argument lists all parameters of that struture, except if it is a switchable structure - in this case, a pulldown menu with switch options is displayed.','type','cellstr','shape','row'));

            
            % split data into reference data (of goal subject) and remaining data
            [refsets,remaining] = utl_collection_closest(args.collection,args.goal_identifier);
            
            % determine the reference weight
            num_reference = length(refsets);
            num_remaining = length(remaining);            
            if args.reference_weight == 0
                args.reference_weight = num_reference/(num_reference+num_remaining); end
            if args.reference_weight > 1
                args.reference_weight = args.reference_weight * num_reference / (num_reference+num_remaining); end
            args.reference_weight = min(1,args.reference_weight);
            
            % determine the weight of each data set
            remaining_weight = 1-args.reference_weight;
            data_weights = [args.reference_weight remaining_weight*ones(1,num_remaining)/num_remaining];
            
            % concatenate all the reference sets across time
            refsets = cellfun(@(s)s.streams(1),refsets);
            reference = set_concat(refsets{:});
            % first pre-process each data set in the corpus and extract features
            corpus = [{struct('streams',{{reference}})}, remaining];
            [features,targets,transforms] = deal(cell(1,length(corpus)));
            for s=length(corpus):-1:1
                if length(corpus{s}.streams) > 1
                    disp_once('Note: ParadigmMKLCSP will use only the first data stream of a recording (no support for multi-modal data).'); end
                % preprocess
                procdata = exp_eval_optimized(flt_pipeline('signal',corpus{s}.streams{1}, args.flt)); %#ok<*NODEF>
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
                X = num2cell(X,[1 2]);
                % calc spatial and temporal pre-processing matrices
                switch args.cov_type
                    case 'shrink'
                        P = {hlp_diskcache('featuremodels',@cov_shrink,cat(2,X{:})')^args.normalizers(1),hlp_diskcache('featuremodels',@cov_shrink,cat(1,X{:}))^args.normalizers(2)}; 
                    case 'full'
                        P = {cov(cat(2,X{:})')^args.normalizers(1),cov(cat(1,X{:}))^args.normalizers(2)}; 
                    case 'diag'                        
                        P = {diag(var(cat(2,X{:})'))^args.normalizers(1),diag(var(cat(1,X{:})))^args.normalizers(2)};
                    otherwise
                        error('Unsupported covariance type requested.');
                end                
                % apply regularization
                P = {args.spatial_whitening*P{1} + (1-args.spatial_whitening)*eye(procdata.nbchan)*trace(P{1})/procdata.nbchan, args.temporal_whitening*P{2} + (1-args.temporal_whitening)*eye(procdata.pnts)*trace(P{2})/procdata.pnts};
                % extract pre-processed features
                tmpfeatures = zeros([size(X{1}),length(X)]);
                for t=1:length(X)
                    tmpfeatures(:,:,t) = P{1}*X{t}*P{2}; end
                features{s} = tmpfeatures;
                % extract target values
                targets{s} = set_gettarget(procdata);
                % save the transforms
                transforms{s} = P;
            end
            
            % train classifier using multi-task learning
            args.ml.learner.data_weights = data_weights;
            model.predictivemodel = ml_train('data',{features,targets}, args.ml);
            % store some more model parameters
            model.featuremodel.P = transforms{1};
            model.featuremodel.apply_to = args.apply_to;
            model.times = procdata.xmin + (0:procdata.pnts-1)/procdata.srate;
            model.cov = cov(procdata.data(:,:)');            
            % set the filter graph based on the last reference data set
            model.tracking.filter_graph = exp_eval_optimized(flt_pipeline('signal',refsets{end}, args.flt));
            model.chanlocs = procdata.chanlocs;
        end
        
        function predictions = predict(self,bundle,model)
            % extract features
            features = self.feature_extract(bundle.streams{1},model.featuremodel);
            % apply classifier
            predictions = ml_predict(features, model.predictivemodel);
        end
        
        function features = feature_extract(self,signal,featuremodel)
            switch featuremodel.apply_to
                case 'channels'
                    features = signal.data;
                case 'sources'
                    features = signal.srcpot;
                case 'full CSD'
                    features = signal.srcpot_all;
                case 'components'
                    features = reshape((procdata.icaweights*procdata.icasphere)*procdata.data(procdata.icachansind,:),[],procdata.pnts,procdata.trials);
            end
            for t=1:size(features,3)
                features(:,:,t) = featuremodel.P{1}*features(:,:,t)*featuremodel.P{2}; end
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
                'Prediction.FeatureExtraction.PatternPairs', '', 'Prediction.MachineLearning.Learner'};
        end
                
    end
end
            
% (turn off a few editor warnings because some actual implementations are missing in this file)
%#ok<*INUSD,*STOUT,*MANU>
