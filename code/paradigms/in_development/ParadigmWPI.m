classdef ParadigmWPI < ParadigmDataflowSimplified
    % Experimental paradigm for ERPs; Assumes that ERPs are composed of localized propagating wavelets.
    %
    % This is a new approach to ERP analysis which does not require any parameters to be set except for the 
    % overall time window of interest -- all the rest is learned as a solution to a jointly convex optimization
    % problem.
    %
    % The model is laid out as follows. The features of interest are an over-complete continuous wavelet transform
    % of each channel (using complex Gaussian wavelets of a specific order) in the overall time window of interest. 
    % This yields a 3d weight tensor of space (channel), time, and wavelet scale. The model is a logistic regression 
    % from these features onto the output variable (which is a generalized linear model, meaning that it is invariant
    % under volume conduction). The model is group sparse with groups over channels (so it is effectively sparse in the
    % wavelet decomposition space), and smooth in time and scale by means of a Gaussian Markov Random field imposed
    % between any two neighbouring sets of weights. The solution is not a MAP estimate but instead a posterior mean 
    % approximation obtained using a relatively new convex reformulation of variational inference [1].
    %
    % The way in which the sparsity/accuracy tradeoff is implemented is currently not yet ideal, which is why this 
    % method is declared work in progress.
    %
    % Name:
    %   Wave Propagation Imaging, work in progress
    %
    %                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
    %                           2011-10-08
    
    methods
      
        function defaults = preprocessing_defaults(self)
            defaults = {'EpochExtraction',[-1 1],'Resampling',70};
        end
        
        function defaults = machine_learning_defaults(self)
            defaults = {'glm' 'Solver','Conjugate Gradients'};
        end
        
        function [featuremodel,predictivemodel] = calibrate_prediction_function(self,varargin)
            args = arg_define(varargin, ...
                arg_norep('signal'), ...
                arg_sub({'fex','FeatureExtraction'},{},...
                    {arg({'sparsity','ModelSparsity'},1,[0 Inf],'Model sparsity. The higher this value, the fewer wavelet components will be used.') ...
                    arg({'smoothness','ModelSmoothness'},1,[0 Inf],'Model smoothness. The higher this value, the more parameter coupling between time-shifted wavelets is encouraged.') ...
                    arg({'normalizers','NormalizationExponents'},[-0.25,-0.25],[],'Normalization exponents [lhs, rhs]. Two-element array of powers for the left-hand-side and right-hand-side normalization matrices that are applied to the data from the region.','guru',true), ...
                    arg({'subsampling','TimeSubsampling'},1,[1 100],'Temporal subsampling factor. Can speed up the learning.'), ...
                    arg({'finestscale','FinestScale'},2,[0 Inf],'Finest wavelet scale.'), ...
                    arg({'coarsestscale','CoarsestScale'},0.5,[0 Inf],'Coarsest wavelet scale. Relative to the size of the epoch.'), ...
                    arg({'scalesteps','NumScales'},20,uint32([1 100]),'Number of scales. Number of wavelet scales to consider between finest and coarsest; log-spaced.'), ...
                    arg({'waveorder','Wigglyness'},3,[0 Inf],'ERP Wigglyness. Assumed wigglyness of the underlying ERP consituents; integer between 1 to 8.') ...
                    }, 'Parameters for the feature-adaptation function. These parameters control how features are statistically adapted and extracted from the filtered data before they are passed int othe machine learning stage.'), ...
                arg_sub({'ml','MachineLearning'},{'Learner',self.machine_learning_defaults()},@ml_train,'Machine learning stage of the paradigm. Operates on the feature vectors that are produced by the feature-extraction stage.'));
            
            % precompute a few wavelet properties
            featuremodel.scales = logspace(log10(args.fex.finestscale),log10(size(args.signal.data,2)*args.fex.coarsestscale),args.fex.scalesteps);
            featuremodel.subsampling = args.fex.subsampling;
            featuremodel.family = sprintf('cgau%i',args.fex.waveorder);
            [featuremodel.wavelet,featuremodel.waveletX] = intwave(featuremodel.family,10);
            featuremodel.waveletStep = featuremodel.waveletX(2)-featuremodel.waveletX(1);
            featuremodel.waveletMax = featuremodel.waveletX(end)-featuremodel.waveletX(1);
            featuremodel.wavelet = conj(featuremodel.wavelet);
            % do a decomposition
            X = self.wpi_decompose(args.signal,featuremodel);
            X = num2cell(X,[1 2]);
            featuremodel.P = {cov(cat(2,X{:})')^args.fex.normalizers(1),var(cat(1,X{:})).^args.fex.normalizers(2)};
            featuremodel.chanlocs = args.signal.chanlocs;
            featuremodel.times = args.signal.xmin + (0:args.signal.pnts-1)/args.signal.srate;
            featuremodel.sparsity = args.fex.sparsity;
            featuremodel.smoothness = args.fex.smoothness;
            % track some properties for inspection
            global tracking; tracking.inspection.wpi_model = featuremodel;
            
            % update machine learning settings
            args.ml.learner.setupfcn = @(varargin)self.wpi_setup(varargin{:},featuremodel.scales,args.fex);
            args.ml.learner.scaling = 'none';
            % extract features & target labels
            features = self.feature_extract(args.signal, featuremodel);
            targets = set_gettarget(args.signal);
            % run the machine learning stage
            predictivemodel = ml_train('data',{features,targets}, args.ml);
        end        
        
        function features = feature_extract(self,signal,featuremodel)
            features = self.wpi_decompose(signal,featuremodel);
            for t=1:size(features,3)
                features(:,:,t) = featuremodel.P{1}*bsxfun(@times,features(:,:,t),featuremodel.P{2}); end
        end
        
        function X = wpi_decompose(self,data,mdl)
            N = length(mdl.scales);
            [C,S,T] = size(data.data);
            % stack all time courses as an array of columns
            X = zeros(S,C,T);
            for t = 1:T
                X(:,:,t) = data.data(:,:,t)'; end
            X = reshape(X,S,C*T);
            Y = zeros(S-1,C*T,N);
            for s = 1:N
                a = mdl.scales(s);
                j = 1+floor((0:a*mdl.waveletMax)/(a*mdl.waveletStep));
                if isscalar(j) j = [1 1]; end
                f = fliplr(mdl.wavelet(j));
                Y(:,:,s) = -sqrt(a) * diff(conv2(X,f','same'));
            end
            Y = reshape(Y,S-1,C,T,N);
            Y = permute(Y,[2 1 3 4]);
            X = zeros(C,2*(S-1)*N,T);
            for t=1:size(data.data,3)
                tmp = reshape(Y(:,:,t,:),C,[]);
                X(:,:,t) = [real(tmp) imag(tmp)];
            end
        end
        
        function [X,y,B,pot,tau,G] = wpi_setup(self,trials,targets,shape,opts,scales,args)
            [n,f] = size(trials);   % number of trials, number of features
            c = shape(1);           % number of channels
            w = shape(2);           % number of wavelet coeffs
            s = length(scales);     % number of scales
            t = w/s/2;              % number of time points (per scale); note that this is to be taken *2 (as we have real & imag)
            
            if strcmp(opts.ptype,'regression')
                error('Currently, only the classification model is implemented.'); end
            if t-round(t) > eps
                error('Some of the design parameters (trials,shape,scales) do not match up as they should.'); end
            
            % inverse link variances for each scale of wavelets
            vars = 1 ./ ((2./scales).^2);
            % variances expanded for all time points (setting the inv variance for the last time point's edge, going to the next row's time point, to 0)
            timevars = [repmat(vars,t-1,1); zeros(1,length(vars))]; timevars = timevars(:)';
            % variances expanded for all channels
            chantimevars = repmat(timevars,c,1); chantimevars = chantimevars(:)';
            % and replicated for the real/imag blocks
            fulltimevars = [chantimevars chantimevars]';
            
            % X is a finite-difference operator on features for smoothing...
            X = spdiags([fulltimevars [zeros(c,1); -fulltimevars(1:end-c,1)]],[0,c]',f,f)*args.smoothness;
            % y is the mean of the Gaussian prior on finite feature differences
            y = zeros(f,1);
            
            % B is the matrix that maps onto the factorial super-Gaussian potentials
            % first, we carry over the features as they are, followed by the projected trials...
            B = [speye(f); double(trials)];
            % we then impose the grouping structure over channels using matrix G...
            % with an additional identity block of size NxN to carry through the features for the logreg potential
            grp = repmat({sparse(ones(1,c))},1,t*s*2);
            G = blkdiag(grp{:},speye(n));
            % now the corresponding tau's for the group sparsity & logistic regression
            tau  = [ones(t*s*2,1)/args.sparsity; targets];
            % finally the potential function, which is a concatentation of Laplacian and logistic
            pot = @(ss,varargin) potCat(ss,varargin{:},{@potLaplace,@potLogistic},{1:(t*s*2),(t*s*2)+(1:n)});            
        end
        
        function layout = dialog_layout_defaults(self)
            layout = {'SignalProcessing.Resampling.SamplingRate','SignalProcessing.EpochExtraction','', ...
                'Prediction.FeatureExtraction.TimeSubsampling','Prediction.FeatureExtraction.FinestScale', ...
                'Prediction.FeatureExtraction.CoarsestScale','Prediction.FeatureExtraction.NumScales','', ...
                'Prediction.MachineLearning.Learner.Lambdas','Prediction.MachineLearning.Learner.Type', ...
                'Prediction.MachineLearning.Learner.Solver'};
        end
        
    end
end

