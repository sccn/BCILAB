function signal = flt_fourier(varargin)
% Transform an epoched data set into a fourier representation.
% Signal = flt_fourier(Signal, Representation, Normalized, LogTransform, Decorrelate)
%
% The fourier representation is practical if there is a highly complex (relevant) spectral structure
% in the data. In any representation except for 'complex', it is necessary that an appropriate
% spatial filter has been applied to the data beforehand (e.g., ICA with the 'transform' option set
% to 1, or the surface laplacian). The complex representation can be spatially transformed after
% running the Fourier filter, since it is a linear operator. If classifiers operate directly on the
% fourier representation, it is usually necessary to chose a representation that can be linearly
% mapped to a prediction, such as 'amplitude', 'polar' or 'phase'; 'complex' requires a non-linear
% classification, which is almost certain to overfit any other random non-linear dependencies in the
% data. Some interesting assumptions that can be imposed on Fourier data include sparsity (l1) and
% group sparsity (if independent component activity was transformed).
%
% The raw fourier spectrum is fairly noisy, but fast. The multitaper approaches yield a smoothed
% (and thus less noisy) estimate -- at significant computational cost. The Welch method gives
% very accurate and low-noise estimates, but is enormously slow.
%
% In:
%   Signal         :   Epoched data set to be processed
%
%   Representation : desired representation; can be one of the following:
%                    'complex'   : the epoch contains the complex spectrum (default)
%                    'power'     : the epoch contains the power spectrum
%                    'amplitude' : the epoch contains the amplitude spectrum (sqrt(of power spectrum)
%                    'phase'     : the epoch contains the phase spectrum
%                    'polar'     : the epoch contains the amplitude/phase spectrum
%                                  (first N channels: amplitude, remaining N channels: phase)
%                    'welch'     : the epoch contains the Welch PSD [1] (which is really slow!); can
%                                  also be specified as a cell array 
%                                  {'welch', 'name',value, 'name',value ...} for names:
%                                  * WindowLength: the WINDOW parameter of pwelch()
%                                  * Overlapping: the NOVERLAP parameter of pwelch()
%                                  * FFTLength: the NFFT parameter of pwelch()
%                    'multitaper': use multi-taper PSD estimation [2,3] (via Chronux); like welch, 
%                                  can be specified with additional parameters: 
%                                  * TimeBandwidth: controls the smoothing
%                                  * Tapers: the number of tapers to use (<= 2*TimeBandwidth-1)
%                                  * Padding: FFT padding (oversampling) factor; 0 is next power of 2
%
%   Normalized     : whether the spectrum should be normalized by 1/f (default: true)
%
%   LogTransform   : whether to transform the result into the log domain (default: false)
%
%
% Out: 
%   Signal  :   processed data set
%
% Notes:
%   Multi-taper estimation is relatively costly; for 32 channel data and 2-second windows, a BCI 
%   could only be updated at ~5 Hz on a typical machine.
%
% Examples:
%   % use default settings (produces frequency-normalized power spectral density features)
%   eeg = flt_fourier(eeg)
%
%   % transform into the complex domain (requires a strong classifier)
%   eeg = flt_fourier(eeg,'complex')
%
%   % transform into the log-power domain (log(abs(fft(X))))
%   eeg = flt_fourier(eeg,'power',true,true)
%
%   % as before, but do not normalize the data by 1/f
%   eeg = flt_fourier(eeg,'power',false,true)
%
%   % as before, but pass all arguments by name
%   eeg = flt_fourier('Signal',eeg,'Representation','power','Normalized',false,'LogTransform',true)
%
%   % transform into amplitude/phase (polar) coordinates
%   eeg = flt_fourier(eeg,'polar')
%
%   % use the Welch method to get more robust estimates
%   eeg = flt_fourier(eeg,'welch')
%
%   % as before, but specify a custom window length
%   eeg = flt_fourier(eeg,{'welch','WindowLength',64})
%
%   % as before, but specify both a custom window length and a custom Overlap
%   eeg = flt_fourier(eeg,{'welch','WindowLength',64,'Overlapping,48})
%
%   % use the Multi-taper method to get highest-quality outputs, and do a log-transform
%   eeg = flt_fourier(eeg,'multitaper',true,true)
%
%   % as before, but specify a fairly large degree of smoothing
%   eeg = flt_fourier(eeg,{'multitaper','TimeBandwidth',8},true,true)
%
% References:
%  [1] Welch, PD "The Use of Fast Fourier Transform for the Estimation of Power Spectra: A Method Based on Time Averaging Over Short, Modified Periodograms", 
%      IEEE Transactions on Audio Electroacoustics, Volume AU-15 (June 1967), pages 70?73.
%  [2] Thomson, D. J. "Spectrum estimation and harmonic analysis." 
%      In Proceedings of the IEEE, Volume 70 (1982), 1055?1096.
%  [3] Slepian, D. "Prolate spheroidal wave functions, Fourier analysis, and uncertainty ? V: The discrete case." 
%      Bell System Technical Journal, Volume 57 (1978), 1371?430.
%
% See also:
%   fft, mtspectrumc, pwelch
%
% TODO:
%   Try to speed up the computation (GPU?)
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-03-28

% flt_fourier_version<1.0> -- for the cache

if ~exp_beginfun('filter') return; end

% requires epoched data, works best on spatially filtered data
declare_properties('name','SpectralTransform', 'depends','set_makepos', 'follows',{'flt_reconstruct','flt_project','flt_window'}, 'independent_channels',true, 'independent_trials',true);

arg_define(varargin,...
    arg_norep({'signal','Signal'}), ...
    arg_subswitch({'rep','Representation'},'power', ...
        {'complex',{},'power',{},'logpower',{},'amplitude',{},'phase',{},'polar',{}, ...
         'welch',{arg({'wndlen','WindowLength'},0.25,[0 0.1 5 Inf],'Window length. Spectral density is estimated in overlapped windows of this length, and then averaged. If this is a fraction, it is implicitly multiplied by the epoch length.'),...
                  arg({'noverlap','Overlapping'},[],[],'Overlap samples. Number of samples by which windows overlap. By default, this is 50% of the window length.','guru',true), ...
                  arg({'nfft','FFTLength'},[],[],'Number of FFT bins. Larger numbers give a smoother spectrum.','guru',true)}, ...
         'multitaper',{arg({'bandwidth','TimeBandwidth'},3,[],'Spectral smoothing. Controls the bias vs. variance of the spectral estimation. Reasonable values are 1 to 3 (1 being fairly noisy, and 3 being fairly smooth but 5x slower)'), ...
                       arg({'tapers','Tapers'},[],[],'Number of tapers. Should be an integer smaller than 2*TimeBandwith; default 2*TimeBandwidth-1','guru',true), ...
                       arg({'padding','Padding'},0,[],'FFT padding factor. Controls the oversampling of the spectrum; 0 is the next largest power of two, 1 is 2x as much, etc.','guru',true)} ...
        }, 'Spectral representation of the signal. Complex and polar are complete representations, the others extract some aspect of the signal. Complex is the only one that usually requires a strong non-linear classifier. Welch and Multitaper are higher-quality spectral estimation techniques. Note: Welch is very slow, be sure to test online performance before relying on it.'),...
    arg({'normalized','Normalized'}, true, [], 'Normalize the spectrum by 1/f. Doing this has benefits for classifiers that work best with naturally normalized features (e.g. some regularized classifiers).'), ...
    arg({'logtransform','LogTransform'}, false, [], 'Log-Transform. Whether to transform the resulting spectral data into the log domain; can facilitates the use of simple (linear) classifiers.'), ...    
    arg({'logspacing','LogSpacing'}, 0, [], 'Log-Spacing. Whether to sub-sample the data in the log domain. If this is a number (>1) it determines the number of samples taken. If this is a fractional number < 1, it is a fraction of the number of trials.'));

[C,S,T] = size(signal.data);    
if strcmp(rep.arg_selection,'welch')
    % compute the PSD using the Welch method
    if ~isempty(rep.wndlen) && rep.wndlen <=1
        rep.wndlen = floor(rep.wndlen*S); end
    [tmp,F] = pwelch(signal.data(1,:,1),rep.wndlen,rep.noverlap,rep.nfft,signal.srate);
    tmp = zeros(C,length(tmp),T);
    for c = 1:C
        for t = 1:T
            [tmp(c,:,t),F] = pwelch(signal.data(c,:,t),rep.wndlen,rep.noverlap,rep.nfft,signal.srate); end
    end
    signal.data = tmp;
    signal.pnts = size(signal.data,2);
    % normalize the data
    if normalized
        signal.data = bsxfun(@times,signal.data,F([2 2:end])'); end
elseif strcmp(rep.arg_selection,'multitaper')
    try
        % compute the PSD using Slepian tapers
        if isempty(rep.tapers)
            rep.tapers = 2*rep.bandwidth-1; end
        if T == 1
            % one-shot calculation
            [signal.data,F] = mtspectrumc(signal.data',struct('tapers',[2*rep.bandwidth rep.tapers],'pad',rep.padding,'Fs',signal.srate));
            signal.data = signal.data';
        else
            % process data over trials
            [first,F] = mtspectrumc(signal.data(:,:,1)',struct('tapers',[2*rep.bandwidth rep.tapers],'pad',rep.padding,'Fs',signal.srate));
            first = first';
            tmp = zeros(C,size(first,2),T);
            tmp(:,:,1) = first;
            for t = 2:size(tmp,3)
                tmp(:,:,t) = mtspectrumc(signal.data(:,:,t)',struct('tapers',[2*rep.bandwidth rep.tapers],'pad',rep.padding,'Fs',signal.srate))'; end
            signal.data = tmp;
        end
    catch e
        if strcmp(e.identifier,'MATLAB:UndefinedFunction')
            error('BCILAB:flt_iir:no_license','Apparently you don''t have a Signal Processing Toolbox license, so you cannot use the ''multitaper'' mode of the SpectralTransform filter.\nYou can switch to a mode that runs under vanilla MATLAB in the "Review/Edit approach" dialog by selecting ''power'' in the SpectralTransform>Representation item instead (which is somewhat lower-quality). If you are using a standard paradigm, you may also look for an equivalent of it that does not require the SigProc toolbox (these are at the bottom of the list under "New Approach").');
        else
            rethrow(e);
        end
    end
    signal.pnts = size(signal.data,2);
    % normalize the data
    if normalized
        signal.data = bsxfun(@times,signal.data,([1 1:size(signal.data,2)-1])/size(signal.data,2)); end
else
    % first transform into the spectrum
    F = (0:ceil(S/2))*(2*signal.srate/S);
    signal.data = fft(signal.data,[],2);
    if normalized
        signal.data = bsxfun(@times,signal.data,([1 1:size(signal.data,2)-1])/size(signal.data,2)); end
    
    % then map to the desired representation
    switch rep.arg_selection
        case 'power'
            signal.data = abs(signal.data);
        case 'amplitude'
            signal.data = real(sqrt(abs(signal.data)));
        case 'phase'
            signal.data = real(angle(signal.data));
        case 'polar'
            signal.data = real([sqrt(abs(signal.data)); angle(signal.data)]);
            % replicate & rename channels
            labels = [cellfun(@(s) [s '_a'], {signal.chanlocs.labels},'UniformOutput',false) cellfun(@(s) [s '_p'], {signal.chanlocs.labels},'UniformOutput',false)];
            signal.chanlocs = signal.chanlocs([1:end 1:end]);
            [signal.chanlocs.labels] = labels{:};
        case 'complex'
            % nothing to do
        otherwise
            error('unsupported representation selected.');
    end
    
    % restrict to the first half of the spectrum
    signal.data = signal.data(:,1:size(signal.data,2)/2,:);
end

if logtransform
    signal.data = log(signal.data); 
    signal.data(~isfinite(signal.data(:))) = 0;
end

if logspacing
    if logspacing <= 1
        logspacing = round(logspacing*size(signal.data,2)); end
    idx = unique(round(logspace(log10(3),log10(size(signal.data,2)),logspacing)));
    signal.data = signal.data(:,idx,:);
    F = F(idx);
    signal.pnts = numel(idx);    
end

signal.freqs = F;

exp_endfun;
