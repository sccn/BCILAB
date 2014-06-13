classdef ParadigmSpectralmeans < ParadigmDataflowSimplified
    % Conventional paradigm for stationary oscillatory processes, using per-channel frequency band averages.
    %
    % This method is essentially the Fourier domain equivalent of para_windowmeans. Since spectral
    % power is not a linear measure of the signal, a spatial filter can significantly improve the
    % performance of this method over simple channel-space band power. Some of the applicable spatial
    % filters are the surface Laplacian and the ICA.
    %
    %
    % Name:
    %   Spectral Means
    %
    %                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
    %                               2011-05-26
    
    methods
      
        function defaults = preprocessing_defaults(self)
            defaults = {'FIRFilter',{[0.5 2],'highpass'}, 'EpochExtraction',[0.5 3.5], 'SpectralTransform',{'multitaper',true,false}, 'Resampling',100};
        end
                
        function model = feature_adapt(self,varargin)
            arg_define(varargin, ...
                arg_norep('signal'), ...
                arg({'wnds','FreqWindows'},[0.5 3; 4 7; 8 12; 13 30; 31 42],[0 0.5 200 1000],'Frequency bands to take as features. Matrix containing one row for the start and end of each frequency band over which the signal mean (per every channel) is taken as a feature. Values in Hz.','cat','Feature Extraction'), ...
                arg({'logtransform','LogTransform'},false,[],'Log-transform output. Log-transformed spectra are more likely to be separable by a linear classifier.'));
            model.wnds = wnds;
            model.chanlocs = signal.chanlocs;
            model.logtransform = logtransform;
        end
        
        function features = feature_extract(self,signal,featuremodel)
            [dummy, idx] = utl_nearest(signal.freqs, featuremodel.wnds); %#ok<ASGLU>
            features = reshape(utl_picktimes(signal.data,idx),[],size(signal.data,3))';
            if featuremodel.logtransform
                features = log(features); end
        end
        
        function visualize_model(self,varargin) %#ok<*INUSD>
            args = arg_define([0 3],varargin, ...
                arg_norep({'parent','Parent'},[],[],'Parent figure.'), ...
                arg_norep({'fmodel','FeatureModel'},[],[],'Feature model. This is the part of the model that describes the feature extraction.'), ...
                arg_norep({'pmodel','PredictiveModel'},[],[],'Predictive model. This is the part of the model that describes the predictive mapping.'), ...
                arg({'paper','PaperFigure'},false,[],'Use paper-style font sizes. Whether to generate a plot with font sizes etc. adjusted for paper.'));
            arg_toworkspace(args);
            parent = args.parent;
            
            % no parent: create new figure
            if isempty(parent)
                parent = figure('Name','Per-window weights'); end
            
            % number of pairs, and index of pattern per subplot
            np = size(fmodel.wnds,1);
            horz = ceil(sqrt(np));
            vert = ceil(np/horz);
            
            % get the weights
            if isfield(pmodel.model,'w')
                weights = pmodel.model.w;
            elseif isfield(pmodel.model,'W')
                weights = pmodel.model.W;
            elseif isfield(pmodel.model,'weights')
                weights = pmodel.model.weights;
            else
                error('Cannot find model weights.');
            end
            
            % check if weights contains a bias value
            if numel(weights)==length(fmodel.chanlocs)*np+1
                weights = weights(1:end-1);
            elseif numel(weights)~=length(fmodel.chanlocs)*np
                error('The model is probably not linear');
            end
            
            % turn into matrix, and optionally convert to forward projections
            weights = reshape(weights,length(fmodel.chanlocs),np);            
            
            % for each window...
            for p=1:np
                subplot(horz,vert,p,'Parent',parent);
                topoplot(weights(:,p),fmodel.chanlocs,'maplimits',[-max(abs(weights(:))) max(abs(weights(:)))]);
                t=title(['Window' num2str(p) ' (' num2str(fmodel.wnds(p,1)) 's to ' num2str(fmodel.wnds(p,2)) 's)']);
                if args.paper
                    set(t,'FontUnits','normalized');
                    set(t,'FontSize',0.1);                    
                    set(gca,'FontUnits','normalized');
                    set(gca,'FontSize',0.1);
                end         
            end
        end
        
        function layout = dialog_layout_defaults(self)
            layout = {'SignalProcessing.Resampling.SamplingRate', 'SignalProcessing.FIRFilter.Frequencies', ...
                'SignalProcessing.EpochExtraction', 'SignalProcessing.SpectralTransform.Representation.TimeBandwidth', ...
                'SignalProcessing.SpectralTransform.Normalized', 'SignalProcessing.SpectralTransform.LogTransform', '', ...
                'Prediction.FeatureExtraction.FreqWindows', '', 'Prediction.MachineLearning.Learner'};
        end
        
    end
end
