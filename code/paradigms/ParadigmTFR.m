classdef ParadigmTFR < ParadigmDataflowSimplified
    % Time/Frequency Regression. This is a new approach that is not yet published.
    %
    % This method learns a linear classifier (or regressor) of second-order dynamics in the EEG over
    % a particular set of time and frequency bins. The key trick lies in a set of extra assumptions
    % that make the (very high-dimensional) solution statistically tractable. The algorithm amounts
    % to a single convex optimization problem dependent on a small number of hyper-parameters that can
    % be optimized using grid search. For a reasonably exhaustive grid the solution is the globally
    % optimal second-order dynamics solution.
    %
    % Notes:
    %  This implementation is not yet optimally tuned -- the spectral estimation should be replaced by a 
    %  similar approach to multi-taper CSP. 
    %
    % Examples:
    %   % learn a spectral classifier for features within -2 to 2 seconds relative to some marker in some standard EEG frequency bands
    %   % here using area under curve (AUC) to optimize the regularization parameter (assuming that the number of exemplars per class is imbalanced)
    %   % and setting the set of relative regularization term weights that should be searched over to {[1 1 1 1]}, i.e., not searching over those, for speed.
    %   myapproach = {'TFR' 'SignalProcessing',{'EpochExtraction',[-2 2]},'Prediction',{'FeatureExtraction',{'Times',[-1.5,-1,-0.5,0,0.5,1,1.5],'Frequencies',[4 7 11 14 20 30]}, ...
    %       'MachineLearning',{'Learner',{'proximal','LambdaSearch',{'ParameterMetric','auc'},'TermWeights',{[1 1 1 1]}}}}};
    %
    % Name:
    %   Time/Frequency Regression
    %
    %                            Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
    %                            2013-03-08    
    
    methods
      
        function defaults = preprocessing_defaults(self)
            % define the default pre-processing parameters of this paradigm
            defaults = {'FIRFilter',{[0.5 1],'highpass'}, 'EpochExtraction',[-3 3], 'Resampling',100};
        end
        
        function defaults = machine_learning_defaults(self)
            % global sharing approach
            defaults = {'proximal', 'lambdaSearch',{'lambdas',2.^(5:-0.1:-4),'foldmargin',0},'regularizers',{ ...
                    'term1',{'trace'}, ...                                                                      % the weights are a small sum of second-order spatial filters per T/F resel
                    'term2',{'trace','LinearOperator','@(x)reshape(x,a*b,c*d)'}, ...                            % across all of time/frequency we learn combinations of few latent spatial filters
                    'term3',{'l2', 'LinearOperator','@(x)vec(diff(x,[],4))','NonorthogonalTransform',true}, ... % temporal smoothness
                    'term4',{'l2', 'LinearOperator','@(x)vec(diff(x,[],3))','NonorthogonalTransform',true}, ... % spectral smoothness
                }, 'regweights',{[1 1 1 1],[1 1 0.5 0.5],[1 1 2 2],[1 1 0.25 0.25],[1 1 4 4],[1 2 1 2],[2 1 1 2],[1 2 2 1],[2 1 2 1]}};
        end
                
        function model = feature_adapt(self,varargin)
            % adapt a feature representation using the CSP algorithm
            args = arg_define(varargin, ...
                arg_norep('signal'), ...
                arg({'freqpoints','Frequencies'},[4 7 11 14 20 30],[0 0.2 200 1000],'Frequency points to consider.'),...
                arg({'timepoints','Times'},-2.5:0.5:2.5,[],'Time points to consider.'), ...
                arg({'from_edges','FromEdges'},true,[],'Run up to epoch edges. This is both in time and frequency.'), ...
                arg({'vectorize_features','VectorizeFeatures'},false,[],'Vectorize the features. For compatibility with basic classifiers.'));
            model.args = rmfield(args,'signal');
            model.chanlocs = args.signal.chanlocs;
        end
        
        function [features,shape] = feature_extract(self,signal,featuremodel)
            % first perform the time/freq decomposition
            args = featuremodel.args;
            if args.from_edges
                args.timepoints = [signal.xmin args.timepoints signal.xmax]; 
                args.freqpoints = [0 args.freqpoints signal.srate/2]; 
            end
            time2idx = @(t) min(signal.pnts,max(1,1+round((t-signal.xmin)*signal.srate)));
            half_hann = @(a,b) 0.5*(1-cos(pi*((((a+1):b) - a) / (b - a))));
            shape = [size(signal.data,1),size(signal.data,1),length(args.freqpoints)-2,length(args.timepoints)-2];
            features = zeros([shape signal.trials]);
            for t=2:length(args.timepoints)-1
                % calculate time window function
                last_tp = time2idx(args.timepoints(t-1));
                cur_tp = time2idx(args.timepoints(t));
                next_tp = time2idx(args.timepoints(t+1));
                wndrange = last_tp+1:next_tp;
                wndfunc = [half_hann(last_tp,cur_tp) 1-half_hann(cur_tp,next_tp)];
                % window and fourier-transform the data
                X = fft(bsxfun(@times,wndfunc,signal.data(:,wndrange,:)),[],2);
                % multiply out the cross-spectral covariance
                tmp = zeros(size(X,1),size(X,1),ceil(size(X,2)/2),size(X,3));
                freq2idx = @(f) min(size(tmp,3),max(1,1+round(f*length(wndrange)/signal.srate)));
                for f=1:size(tmp,3)
                    for n=1:size(tmp,4)
                        tmp(:,:,f,n) = 2*real(X(:,f,n)*X(:,f,n)'); 
                    end
                end
                % average around frequency centers
                for f=2:length(args.freqpoints)-1
                    last_fp = freq2idx(args.freqpoints(f-1));
                    cur_fp = freq2idx(args.freqpoints(f));
                    next_fp = freq2idx(args.freqpoints(f+1));
                    wndrange = last_fp+1:next_fp;
                    wndfunc = [half_hann(last_fp,cur_fp) 1-half_hann(cur_fp,next_fp)];
                    features(:,:,f-1,t-1,:) = mean(bsxfun(@times,tmp(:,:,wndrange,:),reshape(wndfunc,1,1,[])),3);
                end
            end            
            % do final vectorization if desired
            if featuremodel.args.vectorize_features
                features = reshape(features,[],signal.trials)'; end            
        end
        
        function layout = dialog_layout_defaults(self)
            layout = {'SignalProcessing.Resampling.SamplingRate','SignalProcessing.FIRFilter.Frequencies', ...
                '', 'SignalProcessing.EpochExtraction','', ...
                'Prediction.FeatureExtraction.Frequencies', 'Prediction.FeatureExtraction.Times', ...
                'Prediction.FeatureExtraction.FromEdges','', ...
                'Prediction.MachineLearning.Learner.LossType', '', ...
                'Prediction.MachineLearning.Learner.LambdaSearch.Lambdas', ...
                'Prediction.MachineLearning.Learner.TermWeights', ...
                'Prediction.MachineLearning.Learner.LambdaSearch.ParameterMetric'};
        end
        
    end
end

