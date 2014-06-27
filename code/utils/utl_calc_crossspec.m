function [result,result2] = utl_calc_crossspec(varargin)
% Compute (average) multi-taper cross-spectrum (frequencies x channels x channels) across trials
% [CrossSpec,WeightedCrossSpec] = utl_calc_crossspec(Signal,FrequencyRange,TimeBandwidth,Tapers,Padding,SubsampleSpectrum,RobustEstimation,SumWeights,FilteredSubsampling,FeatureFilters)
%
% In:
%   Signal : epoched EEGLAB data set struct to analyze
%
%   FrequencyRange : Frequency range of interest. This is the overall frequency range within which
%                    to compute the cross spectrum. (default: [1 45])
%
%   TimeBandwidth : Spectral smoothing. Controls the bias vs. variance of the spectral estimation.
%                   Reasonable values are 1 to 3 (1 being fairly noisy, and 3 being fairly smooth 
%                   but 5x slower) (default: 5)
%
%   Tapers : Number of tapers. Should be an integer smaller than 2*TimeBandwith; 
%            (default: 2*TimeBandwidth-1)
%
%   Padding : FFT padding factor. Controls the oversampling of the spectrum; 0 is the next largest
%             power of two, 1 is 2x as much, etc. (default: 0)
%
%   SubsampleSpectrum : Sub-sample the spectrum. This is the factor by which to sub-sample.
%                       (default: 1)
%
%   RobustEstimation : Robust cross-spectral estimation. Whether cross-spectral matrices should be
%                      aggregated across trials in a robust manner (default: false)
%
%   SumWeights : Weights for weighted averaging. If passed, the weighted average cross-spectrum will
%                be returned as a second output. (default: [])
%
%   FilteredSubsampling : Use a filter prior to sub-sampling. Slower but yields a better spectral
%                         estimate. (default: false)
%
%   FeatureFilters : Filter tensor for feature extraction. If empty, average cross-spectra are
%                    computed. If false, the raw cross-spectra are returned.
%
% Out:
%   CrossSpec : the resulting cross-spectrum
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2013-04-15

% utl_calc_crossspec_version<1.2> -- for the cache

args = arg_define(varargin, ...
    arg_norep('signal'), ...
    arg({'freqwnd','FrequencyRange'},[1 45],[],'Frequency range of interest. This is the overall frequency range within which to compute the cross spectrum.'), ...
    arg({'bandwidth','TimeBandwidth'},5,[],'Spectral smoothing. Controls the bias vs. variance of the spectral estimation. Reasonable values are 1 to 3 (1 being fairly noisy, and 3 being fairly smooth but 5x slower)'), ...
    arg({'tapers','Tapers'},[],[],'Number of tapers. Should be an integer smaller than 2*TimeBandwith; default 2*TimeBandwidth-1'), ...
    arg({'padding','Padding'},0,[],'FFT padding factor. Controls the oversampling of the spectrum; 0 is the next largest power of two, 1 is 2x as much, etc.'), ...
    arg({'subsample_spectrum','SubsampleSpectrum'},1,[],'Sub-sample the spectrum. This is the sub-sampling factor.'), ...
    arg({'robust_estimation','RobustEstimation'},false,[],'Robust cross-spectral estimation. Whether cross-spectral matrices should be aggregated across trials in a robust manner.'), ...
    arg({'sum_weights','SumWeights'}, [],[],'Weights for weighted averaging. If passed, the weighted average cross-spectrum will be returned as a second output.'), ...
    arg({'filtered_subsampling','FilteredSubsampling'},false,[],'Use a filter prior to sub-sampling. Slower but yields a better spectral estimate.'), ...
    arg_nogui({'feature_filters','FeatureFilters'}, [],[],'Filter tensor for feature extraction. If empty, average cross-spectra are computed. If false, the raw cross-spectra are returned.'));

signal = args.signal;
signal.data = double(signal.data);
if isempty(args.tapers)    
    args.tapers = round(2*args.bandwidth-1);
elseif ~isscalar(args.tapers) || round(args.tapers) ~= args.tapers || args.tapers < 1
    error('The tapers argument, if given, must be a scalar integer.');
end
if args.subsample_spectrum > 1
    N = 2 * args.subsample_spectrum;
    wnd = 0.5 * (1-cos(2*pi*(0:(N-1))/(N-1)));
end

[C,S,T] = size(signal.data); %#ok<ASGLU>
if ~isempty(args.feature_filters)
    if ~isequal(args.feature_filters,false)
        % apply feature-extraction filters
        flts = permute(args.feature_filters,[2 3 1]);
        for t = T:-1:1
            % calculate single-trial cross-spectrum
            cs = abs(CrossSpecMatc(signal.data(:,:,t)',signal.pnts/signal.srate,struct('tapers',[2*args.bandwidth args.tapers],'pad',args.padding,'Fs',signal.srate,'fpass',args.freqwnd)));
            % filter and sub-sample
            if args.subsample_spectrum > 1 && args.filtered_subsampling
                cs = filter(wnd,1,cs,1); end
            cs = permute(cs(1:args.subsample_spectrum:end,:,:),[2 3 1]);
            % extract spectral features
            for f=size(cs,3):-1:1
                result(f,:,t) = log(diag(flts(:,:,f)' * cs(:,:,f) * flts(:,:,f))); end
        end
    else
        % emit raw cross-spectra
        for t = T:-1:1
            cs = abs(CrossSpecMatc(signal.data(:,:,t)',signal.pnts/signal.srate,struct('tapers',[2*args.bandwidth args.tapers],'pad',args.padding,'Fs',signal.srate,'fpass',args.freqwnd)));
            if args.subsample_spectrum > 1 && args.filtered_subsampling
                cs = filter(wnd,1,cs); end
            spec(:,:,:,t) = cs(1:args.subsample_spectrum:end,:,:);
        end
        result = spec;
    end
else
    if ~args.robust_estimation
        % average the resampled cross-spectra on the fly
        if ~isempty(args.sum_weights) && length(args.sum_weights) ~= T
            error('BCILAB:utl_calc_crossspec:incompatible_lengths','The length of the SumWeights parameter needs to be equal to the number of trials in the data'); end            
        meanspec = [];
        weightedspec = [];
        for t = 1:T
            cs = abs(CrossSpecMatc(signal.data(:,:,t)',signal.pnts/signal.srate,struct('tapers',[2*args.bandwidth args.tapers],'pad',args.padding,'Fs',signal.srate,'fpass',args.freqwnd)));
            % filter and sub-sample
            if args.subsample_spectrum > 1 && args.filtered_subsampling
                cs = filter(wnd,1,cs,1); end
            cs = cs(1:args.subsample_spectrum:end,:,:);
            if isempty(meanspec)
                meanspec = cs;
            else
                meanspec = meanspec + cs;
            end
            if ~isempty(args.sum_weights)
                if isempty(weightedspec)
                    weightedspec = cs*args.sum_weights(t);
                else
                    weightedspec = weightedspec + cs*args.sum_weights(t);
                end
            end
        end
        result = meanspec/size(signal.data,3);
        if ~isempty(args.sum_weights)
            result2 = weightedspec; end
    else
        if ~isempty(args.sum_weights)
            error('Robust weighted averaging is not yet implemented.'); end
        % first generate the cross-spectra, then take a geometric median over them
        for t = T:-1:1
            cs = abs(CrossSpecMatc(signal.data(:,:,t)',signal.pnts/signal.srate,struct('tapers',[2*args.bandwidth args.tapers],'pad',args.padding,'Fs',signal.srate,'fpass',args.freqwnd)));
            % filter and sub-sample
            if args.subsample_spectrum > 1 && args. filtered_subsampling
                cs = filter(wnd,1,cs,1); end
            spec(:,:,:,t) = cs(1:args.subsample_spectrum:end,:,:);
        end
        for f=size(spec,1):-1:1
            fspec = reshape(permute(spec(f,:,:,:),[4 2 3 1]),T,C*C);
            medspec(f,:,:) = reshape(geometric_median(fspec,1.e-5,mean(fspec)),C,C);
        end
        result = medspec;
    end
end
