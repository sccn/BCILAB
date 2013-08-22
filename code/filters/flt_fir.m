function [signal,state] = flt_fir(varargin)
% Filter a continuous data set by a digital FIR filter.
% [Signal,State] = flt_fir(Signal, Frequencies, Mode, Type, PassbandRipple, StopbandRipple, State)
%
% Digital FIR filters [1] are computationally less efficient than IIR filters, but allow for
% somewhat more control. Specifically, FIR filters can not "blow up" (diverge), even if extremely
% demanding frequency responses are implemented (e.g., drift removal). The computational cost of
% very low-frequency filters during online processing may be prohibitive, though. FIR filters can be
% designed with different phase (delay/distortion) behavior, depending on the desired application.
% Linear phase filters are the most commonly used ones, as they do not distort the data (which makes
% interpretation easier) but only delay it, and because they are causal (i.e. can be used online).
% The delay can, however, easily be too large for certain time-sensitive online tasks (it is a
% function of the lower transition edge). Zero-phase filters are mostly interesting for
% visualization, as they neither delay nor distort the signal, but cannot be used for online
% applications, or within data analyses that estimate online application behavior. Minimum-phase
% filters can be used online, have very low latency, and can implement extreme frequency responses,
% but distort the signal. In that case, some assumptions about signal shape may turn invalid, and
% have to be revised given the filtered data.
% 
% In:
%   Signal        :   continuous data set to be filtered
%
%   Frequencies  :   frequency specification:
%                    * for a low/high-pass filter, this is: [transition-start, transition-end],in Hz
%                    * for a band-pass/stop filter, this is: [low-transition-start,
%                      low-transition-end, hi-transition-start, hi-transition-end], in Hz
%                    * for a free-form filter, this is a cell array of {[frequency, frequency, ...],
%                      [amplitude, amplitude, ...]} (where the amplitudes specify piecewise constant 
%                      regions in the desired filter response, usually between 0 and 1, and the 
%                      frequencies are the lower and upper frequency edge of each of the bands, 
%                      omitting the lower edge of the first band and upper edge of the last band, 
%                      which are assumed to be 0Hz and the Nyquist frequency, respectively) 
%
%                      Alternatively, it can also be a 3xN array of the form;
%                      [freq_lo,freq_hi,amp; freq_lo,freq_hi,amp; freq_lo,freq_hi,amp; ...]
%                      specifying the lower edge, upper edge and amplitude of each constant segment.
%                      The lower edge of the first segment and upper edge of the last segment are 
%                      ignored and assumed as explained above.
%
%   Mode         :   filter mode, 'bp' for bandpass, 'hp' for highpass, 'lp' for lowpass, 'bs' for
%                    bandstop, 'ff' for free-form (default: 'bp')
%
%   Type         :   * 'minimum-phase' minimum-hase filter -- pro: introduces minimal signal delay;
%                       con: distorts the signal (default)
%                    * 'linear-phase' linear-phase filter -- pro: no signal distortion; con: delays
%                       the signal
%                    * 'zero-phase' zero-phase filter -- pro: no signal delay or distortion; con:
%                       can not be used for online purposes
%
%   PassbandRipple  :   Maximum relative ripple amplitude in pass-band. Relative to nominal 
%                       pass-band gain. Assumed to be in db if negative, otherwise taken as a ratio.
%                       (default: -20)
%
%   StopbandRipple  :   Maximum relative ripple amplitude in stop-band. Relative to nominal
%                       pass-band gain. Assumed to be in db if negative, otherwise taken as a ratio.
%                       (default: -30)
%
%   DesignRule : Filter design rule. Parks-McClellan minimizes the maximum error, the Window Method
%                minimizes the square error, and Frequency Sampling constructs the filter via the 
%                Fourier transform without tuning. (default: 'Frequency Sampling')
%
%   ChunkLength : Maximum chunk length. Process the data in chunks of no larger than this (to avoid
%                 memory limitations). (default: 50000)
%
%   NormalizeAmplitude : Normalize amplitude. Normalizes the amplitude such that the maximum gain is
%                        1. This helps with the occasional erratic filter design result. (default: true)
%
%   State        :   previous filter state, as obtained by a previous execution of flt_fir on an
%                    immediately preceding data set (default: [])
%
% Out: 
%   Signal       :  filtered, continuous data set
%   State        :  state of the filter, after it got applied
%
% Examples:
%   % use a 7-30 Hz bandpass filter, with transition regions that are 2 Hz wide
%   eeg = flt_fir(eeg,[6 8 29 31])
%
%   % use a 1Hz highpass filter (with a transition between 0.9 and 1.1 Hz)
%   eeg = flt_fir(eeg,[0.9 1.1],'highpass')
%
%   % use a 1Hz highpass filter (with very generous transition region) that is linear phase (i.e. 
%   % does not distort the signal)
%   eeg = flt_fir(eeg,[0.5 1.5],'highpass','linear-phase')
%
%   % use a 7.5-30 Hz bandpass filter, with transition regions that are 5 Hz wide, and a particular
%   % specification of pass-band and stop-band rippling constraints, passing all arguments by name
%   eeg = flt_fir('Signal',eeg,'Frequencies',[5 10 27.5 32.5],'PassbandRipple',-10,'StopbandRipple',-50);
%
%   % as previous, but using the short argument names
%   eeg = flt_fir('signal',eeg,'fspec',[5 10 27.5 32.5],'passripple',-10,'stopripple',-50);
%
%   % implement a free-form FIR filter with peaks within 12-15 Hz and within 26-35 Hz
%   eeg = flt_fir(eeg,[0 11 0; 12 15 1; 16 25 0; 26 35 1; 36 100 0],'freeform')
%
% Notes:
%   The Parks-McClellan and Window Method design rules require the Signal Processing toolbox.
%
% References:
%   [1] A.V. Oppenheim and R.W. Schafer, "Digital Signal Processing",
%       Prentice-Hall, 1975.
%
% See also:
%   firpm, fir1, design_fir, filter
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-17  

% flt_fir_version<1.03> -- for the cache

if ~exp_beginfun('filter') return; end

% does not make sense on epoched data
declare_properties('name','FIRFilter', 'follows','flt_iir', 'cannot_follow','set_makepos', 'independent_channels',true, 'independent_trials',true);

arg_define(varargin, ... 
    arg_norep({'signal','Signal'}), ...
    arg({'fspec','Frequencies'}, [], [], ['Frequency specification of the filter. For a low/high-pass filter, this is: [transition-start, transition-end], in Hz and for a band-pass/stop filter, this is: [low-transition-start, low-transition-end, hi-transition-start, hi-transition-end], in Hz. For a free-form filter, this is a 2d matrix of the form [frequency,frequency,frequency, ...; amplitude, amplitude, amplitude, ...] or [frequency,frequency,frequency, ...; amplitude, amplitude, amplitude, ...; ripple, ripple, ripple, ...]']), ...
    arg({'fmode','Mode'}, 'bandpass', {'bandpass','highpass','lowpass','bandstop','freeform'}, 'Filtering mode. Determines how the Frequencies parameter is interpreted.'), ...
    arg({'ftype','Type'},'minimum-phase', {'minimum-phase','linear-phase','zero-phase'}, 'Filter type. Minimum-phase introduces only minimal signal delay but distorts the signal. Linear-phase has no signal distortion but delays the signal. Zero-phase has neither signal delay nor distortion but can not be used for online purposes.'), ...
    arg({'passripple','PassbandRipple'}, -20, [-180 1], 'Maximum relative ripple amplitude in pass-band. Relative to nominal pass-band gain. Affects the filter length (i.e., delay). Assumed to be in db if negative, otherwise taken as a ratio.'), ...
    arg({'stopripple','StopbandRipple'}, -40, [-180 1], 'Maximum relative ripple amplitude in stop-band. Relative to nominal pass-band gain. Affects the filter length (i.e., delay). Assumed to be in db if negative, otherwise taken as a ratio.'), ...
    arg({'designrule','DesignRule'}, 'Frequency Sampling', {'Parks-McClellan','Window Method','Frequency Sampling'}, 'Filter design rule. Parks-McClellan minimizes the maximum error, the Window Method minimizes the square error, and Frequency Sampling constructs the filter via the Fourier transform without tuning (the latter requires no sigproc toolbox).'), ...
    arg({'chunk_length','ChunkLength'},50000,[], 'Maximum chunk length. Process the data in chunks of no larger than this (to avoid memory limitations).','guru',true), ...
    arg({'normalize_amplitude','NormalizeAmplitude'},true,[], 'Normalize amplitude. Normalizes the amplitude such that the maximum gain is as desired. This helps with the occasional erratic filter design result.','guru',true), ...
    arg_nogui({'state','State'}));

if isempty(state)
    % design filter kernel
    if passripple < 0
        passripple = 10^(passripple/10); end
    if stopripple < 0
        stopripple = 10^(stopripple/10); end
    
    if ~(iscell(fspec) && isscalar(fspec))
        
        % create a filter specification accepted by firpmord or kaiserord
        switch fmode
            case {'freeform','ff'}
                % a free-form frequency spec is given
                if iscell(fspec)
                    % given as a cell array of {bandfreqs,amps} or {bandfreqs, amps, ripple}
                    if length(fspec{1}) == 2*length(fspec{2})
                        error('When specifying the bands for each constant-amplitude region of the filter response, the first band is assumed to begin at 0Hz and the last band is assumed to end at the Nyquist frequency -- thus, these 2 numbers in the band specification should be omitted.'); 
                    elseif length(fspec{1}) ~= 2*length(fspec{2})-2
                        error('The specification of band edges does not match the specification of band amplitudes; for each band, a lower and an upper edge frequency must be listed, and both the lower edge of the first band and upper edge of the last band must be omitted (they equal 0Hz and the Nyquist frequency, respectively).');
                    end
                elseif ~isvector(fspec)
                    if size(fspec,2) == 3
                        bands = fspec(:,1:2)'; 
                        fspec = {bands(2:end-1),fspec(:,3)'};
                    elseif size(fspec,1) == 4
                        bands = fspec(:,1:2)'; 
                        fspec = {bands(2:end-1),fspec(:,3)',fspec(:,4)'};
                    else
                        error('When specifying the piecewise-constant filter design in matrix form, a 3xB or 4xB matrix (B = number of bands) of the form [freq_lo,freq_hi,amp; freq_lo,freq_hi,amp; freq_lo,freq_hi,amp; ...] or [freq_lo,freq_hi,amp,ripple; freq_lo,freq_hi,amp,ripple; ...] must be given.');
                    end
                else
                    error('When specifying the piecewise-constant filter design in matrix form, a 3xB or 4xB matrix (B = number of bands) of the form [freq_lo,freq_hi,amp; freq_lo,freq_hi,amp; freq_lo,freq_hi,amp; ...] or [freq_lo,freq_hi,amp,ripple; freq_lo,freq_hi,amp,ripple; ...] must be given.');
                end
            case {'bandpass','bp'}
                fspec = {fspec,[0 1 0]};
            case {'bandstop','bs'}
                fspec = {fspec,[1 0 1]};
            case {'lowpass','lp'}
                fspec = {fspec,[1 0]};
            case {'highpass','hp'}
                fspec = {fspec,[0 1]};
            otherwise
                error(['Unrecognized filter mode: ' hlp_tostring(fmode)]);
        end
        
        % is the filter being applied twice? correct for that.
        if strcmp(ftype,'zero-phase')
            fspec{2} = sqrt(fspec{2}); end
        % derive the rippling specification from the amplitudes and passripple/stopripple
        if length(fspec) < 3        
            fspec{3} = stopripple + fspec{2}*(passripple-stopripple); end    

        try
            % design the filter
            switch designrule
                case {'Parks-McClellan','pm'}
                    args = hlp_diskcache('filterdesign',@firpmord,fspec{:},signal.srate,'cell');
                    b = hlp_diskcache('filterdesign',@firpm,max(3,args{1}),args{2:end});
                case {'Window Method','wm'}
                    args = hlp_diskcache('filterdesign',@kaiserord,fspec{:},signal.srate,'cell');
                    b = hlp_diskcache('filterdesign',@fir1,max(3,args{1}),args{2:end});
                case {'Frequency Sampling','fs'}
                    % get frequencies and amplitudes
                    freqs = [0 fspec{1}*2/signal.srate 1];
                    amps = vec([fspec{2}; fspec{2}]);
                    % design Kaiser window for smallest transition region
                    [dummy,pos] = min(diff(freqs)); %#ok<ASGLU>
                    wnd = design_kaiser(freqs(pos),freqs(pos+1),-20*log10(stopripple),amps(end)~=0);
                    % design FIR filter with that window
                    b = design_fir(length(wnd)-1,freqs,amps,[],wnd);
                otherwise
                    error(['Unrecognized filter design rule: ' hlp_tostring(designrule)]);
            end
        catch e
            if strcmp(e.identifier,'MATLAB:UndefinedFunction')
                error('A function was not found (likely from the signal processing toolbox). Consider falling back to the Frequency Sampling method.');
            else
                rethrow(e);
            end
        end
        state.b = b;
        n = length(state.b);
        
        % use cepstral windowing to calculate a minimum-phase filter (note: the min-phase change applies
        if strcmp(ftype,'minimum-phase') && ~any(strcmp(fmode,{'highpass','hp'}))
            wnd = [1 2*ones(1,(n+mod(n,2))/2-1) ones(1,1-mod(n,2)) zeros(1,(n+mod(n,2))/2-1)];
            state.b = real(ifft(exp(fft(wnd.*real(ifft(log(abs(fft(state.b))+stopripple)))))));
        end
        
        % normalize the magnitude
        if normalize_amplitude
            maxamp = max(abs(fft([state.b(:); zeros(1000,1)])));
            state.b = state.b*max(fspec{2})/maxamp*(1+passripple);
        end
    else
        % precomputed filter coefficients
        state.b = fspec{1};
    end 
end

[b,n] = deal(state.b,length(state.b));
% process each known time series field
for fld = {'data','srcpot','icaact'}
    field = fld{1};
    if isfield(signal,field) && ~isempty(signal.(field)) && ~isequal(signal.(field),1)
        if ~isfield(state,field)
            state.(field) = []; end

        % phase 2: filter the data
        sig = double(signal.(field))';
        if isempty(state.(field))
            % no prior state: prepend the signal with a mirror section of itself, to minimize
            % start-up transients (and if the signal is too short, we repeat it as much as we need)
            sig = [repmat(2*sig(1,:),n,1) - sig(1+mod(((n+1):-1:2)-1,size(sig,1)),:); sig];
            if strcmp(ftype,'zero-phase')
                % to get a zero-phase filter, we run the filter backwards first
                % reverse the signal and prepend it with a mirror section (to minimize start-up transients)
                sig = sig(end:-1:1,:); sig = [repmat(2*sig(1,:),n,1) - sig((n+1):-1:2,:); sig];
                % run the filter
                sig = filter_fast(b,1,sig);
                % reverse and cut startup segment again
                sig = sig(end:-1:(n+1),:);
            end
            prepended = true;
        else
            % online case: check for misuses
            if strcmp(ftype,'zero-phase')
                error('zero-phase filters are non-causal and cannot be run online (or on continued data); use linear-phase or minimum-phase filters, or flt_iir.'); end
            prepended = false;
        end
        
        % apply the filter
        S = size(sig,1);
        numsplits = ceil(S/chunk_length);
        for i=0:numsplits-1
            range = 1+floor(i*S/numsplits) : min(S,floor((i+1)*S/numsplits));
            [sig(range,:),state.(field)] = filter_fast(b,1,sig(range,:),state.(field),1);
        end
        
        % cut off the data segment that we previously prepended
        if prepended
            sig(1:n,:) = []; end
        
        % write the data back
        signal.(field) = sig';
    end
end

if ~isfield(signal.etc,'filter_delay')
    signal.etc.filter_delay = 0; end

if strcmp(ftype,'linear-phase')
    signal.etc.filter_delay = signal.etc.filter_delay + length(b)/2/signal.srate;
elseif strcmp(ftype,'minimum-phase')
    signal.etc.filter_delay = signal.etc.filter_delay + hlp_getresult(2,@max,b)/signal.srate;
end    

exp_endfun;
