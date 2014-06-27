function signal = flt_ersp(varargin)
% Calculate the event-related spectral perturbation for an epoched signal.
% Signal = flt_ersp(Signal)
%
% This calculates a time/frequency representation for each channel at the given resolution.
%
% In:
%   Signal : epoched data set to be processed
%
% Out: 
%   Signal : processed data set
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2013-12-06

% flt_ersp_version<1.0> -- for the cache

if ~exp_beginfun('filter') return; end

% requires epoched data, works best on spatially filtered data
declare_properties('name','ERSPTransform', 'depends','set_makepos', 'follows',{'flt_fourier'}, 'independent_channels',true, 'independent_trials',true);

arg_define(varargin,...
    arg_norep({'signal','Signal'}), ...
    arg({'windowlength','WindowLength'},0.5,[0 0.1 5 Inf],'Moving window length. Length of the moving window within which the spectrum shall be calculated (in seconds).'), ...
    arg({'windowstep','WindowStep'},0.1,[0 0.05 1 Inf],'Moving window stepping. Step size by which the moving window shall be moved (in seconds).'), ...
    arg({'freqrange','FrequencyRange'},[2 48],[0 Inf],'Frequency range. The frequency range for which the ERSP shall be calculated.'), ...
    arg({'timebandwidth','TimeBandwidth','bandwidth'},3,[],'Spectral smoothing. Controls the bias vs. variance of the spectral estimation. Reasonable values are 1 to 3 (1 being fairly noisy, and 3 being fairly smooth but 5x slower)'), ...
    arg({'numtapers','Tapers','tapers'},[],uint32([1,1000]),'Number of tapers. Should be an integer smaller than 2*TimeBandwith; default 2*TimeBandwidth-1','guru',true), ...
	arg({'padding','Padding'},0,[],'FFT padding factor. Controls the oversampling of the spectrum; 0 is the next largest power of two, 1 is 2x as much, etc.','guru',true), ...
    arg({'normalized','Normalized'}, true, [], 'Normalize the spectrum by 1/f. Doing this has benefits for classifiers that work best with naturally normalized features (e.g. some regularized classifiers).'), ...
    arg({'spectralmap','SpectralMap','SpectralMapping'}, 'linear', {'linear','sqrt','log'}, 'Spectral mapping. The sqrt and log transformations can help make the features more suitable for linear classifiers.'), ...
    arg({'logspacing','LogSpacing'}, 0, [], 'Log-Spacing. Whether to sub-sample the data in the log domain. If this is a number (>1) it determines the number of samples taken. If this is a fractional number < 1, it is a fraction of the number of trials.'));

if isempty(numtapers)
    numtapers = 2*timebandwidth-1; end

% for each non-empty time-series field...
for f = utl_timeseries_fields(signal)
    if ~isempty(signal.(f{1}))
        X = double(signal.(f{1}));
        dims = size(X);
        % fix NaN's and Inf's
        X(~isfinite(X(:))) = 0;
        % flip dimensions so the time dimension comes first
        X = permute(X,[2 1 3:length(dims)]);
        % calculate time/frequency decomposition with all non-time dimensions collapsed into one
        [X,T,F] = mtspecgramc(X(:,:),[windowlength windowstep],struct('tapers',[2*timebandwidth numtapers],'pad',padding,'Fs',signal.srate,'fpass',freqrange));
        % expand the non-time dimensions again (time x frequency x channels x trials)
        X = reshape(X,[size(X,1) size(X,2) dims([1 3:end])]);
        % permute back into (channels x time x frequency x trials)
        X = permute(X,[3 1 2 4:(1+length(dims))]);
        % normalize spectrum
        if normalized
            X = bsxfun(@times,X,permute(([1 1:size(X,3)-1])/size(X,3),[1 3 2])); end
        % perform spectral mapping
        if strcmp(spectralmap,'log')
            X = log(X+0.001); 
        elseif strcmp(spectralmap,'sqrt')
            X = sqrt(X);
        end
        % perform logspacing
        if logspacing
            if logspacing <= 1
                logspacing = round(logspacing*size(X,3)); end
            idx = unique(round(logspace(log10(3),log10(size(X,3)),logspacing)));
            X = X(:,:,idx,:,:,:,:,:);
            F = F(idx);
        end
        % write back, but permute as channels x time x trials x frequency
        signal.(f{1}) = permute(X,[1 2 4 3]);
    end
end

% update signal meta-data
signal.srate = 1/windowstep;
signal.pnts = size(signal.data,2);
signal.times = signal.xmin+T;
signal.freqs = F;

exp_endfun;
