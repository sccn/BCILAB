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
                arg({'wnds','FreqWindows'},[0.5 3; 4 7; 8 12; 13 30; 31 42],[],'Frequency bands to take as features. Matrix containing one row for the start and end of each frequency band over which the signal mean (per every channel) is taken as a feature. Values in Hz.','cat','Feature Extraction'), ...
                arg({'logtransform','LogTransform'},false,[],'Log-transform output. Log-transformed spectra are more likely to be separable by a linear classifier.'));
            model.wnds = wnds;
            model.chanlocs = signal.chanlocs;
            model.logtransform = logtransform;
        end
        
        function features = feature_extract(self,signal,featuremodel)
            [~, idx] = utl_nearest(signal.freqs, featuremodel.wnds);
            features = reshape(utl_picktimes(signal.data,idx),[],size(signal.data,3))';
            if featuremodel.logtransform
                features = log(features); end
        end
        
        function visualize_model(self,parent,fmodel,pmodel) %#ok<*INUSD>
            % no parent: create new figure
            if isempty(parent)
                parent = figure('Name','Per-window weights'); end
            % number of pairs, and index of pattern per subplot
            np = size(fmodel.wnds,1);
            horz = ceil(sqrt(np));
            vert = ceil(np/horz);
            % for each window...
            for p=1:np
                subplot(horz,vert,p,'Parent',parent);
                % get the weights
                if isfield(pmodel.model,'w')
                    weights = pmodel.model.w;
                elseif isfield(pmodel.model,'W')
                    weights = pmodel.model.W;
                elseif isfield(pmodel.model,'weights')
                    weights = pmodel.model.weights;
                else
                    title('Cannot find model weights.');
                    weights = [];
                end
                % extract appropriate weights portion
                if ~isempty(weights)
                    if length(weights) == np*length(fmodel.chanlocs) || length(weights) == np*length(fmodel.chanlocs)+1
                        subweights = weights(1+(p-1)*length(fmodel.chanlocs) : p*length(fmodel.chanlocs));
                    else
                        title('Model is probably not linear.');
                        subweights = [];
                    end
                end
                % display
                if ~isempty(weights) && ~isempty(subweights)
                    topoplot(subweights,fmodel.chanlocs,'maplimits',[-max(abs(weights)) max(abs(weights))]);
                    title(['Window' num2str(p) ' (' num2str(fmodel.wnds(p,1)) 's to ' num2str(fmodel.wnds(p,2)) 's)']);
                end
            end
        end
        
        function layout = dialog_layout_defaults(self)
            layout = {'SignalProcessing.Resampling.SamplingRate', 'SignalProcessing.IIRFilter.Frequencies', ...
                'SignalProcessing.EpochExtraction', 'SignalProcessing.SpectralTransform.Representation.TimeBandwidth', ...
                'SignalProcessing.SpectralTransform.Normalized', 'SignalProcessing.SpectralTransform.LogTransform', '', ...
                'Prediction.FeatureExtraction.FreqWindows', '', 'Prediction.MachineLearning.Learner'};
        end
        
    end
end
