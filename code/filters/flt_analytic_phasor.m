function [signal,state] = flt_analytic_phasor(varargin)
% Calculate analytic amplitude and phase of a signal.
%
% The process implemented by this filter is as follows:
% 1) filter the original data at some frequency to get the cosine (real)
%    part of the signal (Re).
% 2) construct a dif filter B=firls(N,F,A,W,'differentiator')
% 3) filter (1) with B to get sine (imaginary) part of signal (Im)
% 4) compute analytic amp as A = sqrt(Re.^2 + Im.^2)
% 5) compute phase as Phi = atan2(Im,Re)
%
% This filter applies to all time-series fields of the given signal.
%
% In:
%   Signal : EEGLAB data set structure
%
%   DiffFilter : type of differencing filter to use and sub-options; default: 'hilbert'
%                'hilbert' : Hilbert differentiator; suboptions are:
%                   FilterLength : length of the filter, auto-determined if [] (default: [])
%                   FrequencyBand : Frequency band to use. The 4 numbers characterize the onset,
%                                   pass-band limits, and offset frequency of the filter. 
%                                   (default: [7 8 12 13])
%                   PassbandRipple : Maximum relative rippling in pass-band. Assumed to be in db if
%                                    negative, otherwise taken as a ratio. (default: -20)
%                   StopbandRipple : Maximum relative rippling in stop-band. Assumed to be in db if
%                                    negative, otherwise taken as a ratio. (default: -30)
%                   DesignRule : Filter design rule. Can be either Least-Squares filter design or
%                                Parks-McClellan filter design (default: 'Parks-McClellan')
%                'differentiator' : odd-symmetry differentiator with special weighting; suboptions are:
%                   FilterLength : length of the filter, auto-determined if [] (default: [])
%                   FrequencyBand : Frequency band to use. The 4 numbers characterize the onset,
%                                   pass-band limits, and offset frequency of the filter. 
%                                   (default: [7 8 12 13])
%                   PassbandRipple : Maximum relative rippling in pass-band. Assumed to be in db if
%                                    negative, otherwise taken as a ratio. (default: -20)
%                   StopbandRipple : Maximum relative rippling in stop-band. Assumed to be in db if
%                                    negative, otherwise taken as a ratio. (default: -30)
%                   DesignRule : Filter design rule. Can be either Least-Squares filter design or
%                                Parks-McClellan filter design (default: 'Parks-McClellan')
%                'smooth_diff' : rectangular low-pass differentiation filter; suboptions are:
%                   FilterLength : length of the filter (default: 10)
%
%   OverrideOriginal : Override original data. If checked, the original signals will be replaced by
%                      their analytic amplitudes. (default: true)
%
%   IncludeAnalyticAmplitude : Include analytic amplitude. If checked, extra fields ending in _aamp
%                              are included for each time series field in the signal. (default: false)
%
%   IncludeAnalyticPhase : Include analytic phase. If checked, extra fields ending in _aphase are
%                          included for each time series field in the signal. (default: false)
%
%   State : optionally the initial filter state from a previous invocation
%
% Out:
%   Signal : filtered signal; the original time series fields now contain the analytic amplitude
%            while new time series fields ending in _aphase are added that contain the analytic 
%            phases.
%
%   State : final filter state
% 
%
% Notes:
%	Requires the Signal Processing toolbox for filter design.
%
%                                Tim Mullen, Swartz Center for Computational Neuroscience, UCSD
%                                2012-03-27

if ~exp_beginfun('filter') return; end

% follows IIR/FIR, as it should operate on a clean signal (rather than depend on HF noise, etc.)
declare_properties('name','AnalyticPhasor', 'cannot_follow','set_makepos', 'follows',{'flt_iir','flt_fir'},'independent_channels',true, 'independent_trials',false);

arg_define(varargin,...
    arg_norep({'signal','Signal'}), ...
    arg_subswitch({'flt','DiffFilter','diffFilter'},'hilbert', ...
        {'hilbert' ...
            {arg({'freqband','FrequencyBand'}, [7 8 12 13], [], 'Frequency band to use. The 4 numbers characterize the onset, pass-band limits, and offset frequency of the filter.','type','expression'), ...
             arg({'fltlen','FilterLength'}, [], [], 'Filter length. This determines both the quality and the delay of the differentiator. If left empty, the optimal order will be estimated based on tolerated pass/stopband ripple.','shape','scalar'), ...
             arg({'passripple','PassbandRipple'}, -20, [-180 1], 'Maximum relative ripple amplitude in pass-band. Relative to nominal pass-band gain. Assumed to be in db if negative, otherwise taken as a ratio.'), ...
             arg({'stopripple','StopbandRipple'}, -30, [-180 1], 'Maximum relative ripple amplitude in stop-band. Relative to nominal pass-band gain. Assumed to be in db if negative, otherwise taken as a ratio.'), ...
             arg({'designrule','DesignRule'}, 'Parks-McClellan', {'Least-Squares','Parks-McClellan'}, 'Filter design rule. Can be either least-squares filter design or Parks-McClellan filter design.')}, ...
         'differentiator' ...
            {arg({'freqband','FrequencyBand'}, [7 8 12 13], [], 'Frequency band to use. The 4 numbers characterize the onset, pass-band limits, and offset frequency of the filter.','type','expression'), ...
             arg({'fltlen','FilterLength'}, [], [], 'Filter length. This determines both the quality and the delay of the differentiator. If left empty, the optimal order will be estimated based on tolerated pass/stopband ripple.','shape','scalar'), ...
             arg({'passripple','PassbandRipple'}, -20, [-180 1], 'Maximum relative ripple amplitude in pass-band. Relative to nominal pass-band gain. Assumed to be in db if negative, otherwise taken as a ratio.'), ...
             arg({'stopripple','StopbandRipple'}, -30, [-180 1], 'Maximum relative ripple amplitude in stop-band. Relative to nominal pass-band gain. Assumed to be in db if negative, otherwise taken as a ratio.'), ...
             arg({'designrule','DesignRule'}, 'Parks-McClellan', {'Least-Squares','Parks-McClellan'}, 'Filter design rule. Can be either least-squares filter design or Parks-McClellan filter design.')}, ...
         'smooth_diff' ...
            {arg({'fltlen','FilterLength'}, 10, [], 'Filter length. This determines both the quality and the delay of the differentiator. The default should be fine.','shape','scalar')} ...
        },'Differencing filter. The smooth_diff implementation is currently experimental; hilbert and differentiator should be fine.'), ...
    arg({'override_original','OverrideOriginal'},true,[],'Override original data. If checked, the original signals will be replaced by their analytic amplitudes.'), ...
    arg({'include_aamp','IncludeAnalyticAmplitude'},false,[],'Include analytic amplitude. If checked, extra fields ending in _aamp are included for each time series field in the signal.'), ...
    arg({'include_aphase','IncludeAnalyticPhase'},false,[],'Include analytic phase. If checked, extra fields ending in _aphase are included for each time series field in the signal.'), ...
    arg({'include_afreq','IncludeAnalyticFrequency'},false,[],'Include analytic frequency. If checked, extra fields ending in _afreq are included for each time series field in the signal.'), ...
    arg_nogui({'state','State'},[]));
    
% make up prior state if necessary
if isempty(state)
    % construct differentiator
    switch flt.arg_selection
        case 'smooth_diff'
            state.b_im = hlp_microcache('fdesign',@smooth_diff,flt.fltlen);
            state.b_re = abs(state.b_im);
        case {'hilbert','differentiator'}
            if ~iscell(flt.freqband)
                flt.freqband = {flt.freqband}; end
            nFreqs = length(flt.freqband);
            % build firpm freq spec
            if flt.passripple < 0
                flt.passripple = 10^(flt.passripple/10); end
            if flt.stopripple < 0
                flt.stopripple = 10^(flt.stopripple/10); end            
            if isempty(flt.fltlen)
                % optionally estimate filters order for each band
                for f=nFreqs:-1:1
                    pmspec{f} = hlp_diskcache('filterdesign',@firpmord,flt.freqband{f},[0 1 0],flt.stopripple + [0 1 0]*(flt.passripple-flt.stopripple),signal.srate,'cell'); end
                % make sure that each filter uses the same order (=delay)
                maxlen = max([3,cellfun(@(x)x{1},pmspec)]);
                for f=nFreqs:-1:1
                    pmspec{f}{1} = maxlen; end
            else
                for f=nFreqs:-1:1
                    pmspec{f} = {flt.fltlen,[0 flt.freqband{f}*2/signal.srate 1],[0 0 1 1 0 0]}; end
            end
            % design filters
            for f=nFreqs:-1:1
                switch flt.designrule
                    case {'Parks-McClellan','pm'}                    
                        state.b_im{f} = hlp_diskcache('filterdesign',@firpm,pmspec{f}{:},flt.arg_selection); 
                        state.b_re{f} = hlp_diskcache('filterdesign',@firpm,pmspec{f}{:}); 
                    case {'Least-Squares','ls'}
                        state.b_im{f} = hlp_diskcache('filterdesign',@firls,pmspec{f}{:},flt.arg_selection);
                        state.b_re{f} = hlp_diskcache('filterdesign',@firls,pmspec{f}{:});
                    otherwise 
                        error(['Unrecognized filter design rule:' hlp_tostring(flt.designrule)]);
                end
            end
        otherwise
            error(['Unrecognized differentiator type selected: ' hlp_tostring(flt.arg_selection)])
    end
    % set up initial state
    state.cached_chanlocs = update_chanlocs(signal.chanlocs,nFreqs);
    for fld = utl_timeseries_fields(signal)
        state.(fld{1}) = struct('zi_re',{cell(1,length(flt.freqband))},'zi_im',{cell(1,length(flt.freqband))}); end
end

% filter bandpassed signal with differentiator to get sine part
nFreqs = length(state.b_im);
for fld = utl_timeseries_fields(signal)
    field = fld{1};
    if ~isempty(signal.(field))
        % get rid of NaN's and Inf's
        signal.(field)(~isfinite(signal.(field)(:))) = 0;
        % apply differentiator to get real (cosine) and imaginary (sine) part of signal
        [impart,repart] = deal(cell(1,nFreqs));
        for f=1:nFreqs
            [impart{f},state.(field).zi_im{f}] = filter_fast(state.b_im{f},1,signal.(field),state.(field).zi_im{f},2);
            [repart{f},state.(field).zi_re{f}] = filter_fast(state.b_re{f},1,signal.(field),state.(field).zi_re{f},2);
        end
        % concatenate across all frequencies
        impart = vertcat(impart{:});
        repart = vertcat(repart{:});
        % compute analytic phase and amplitude
        aamp = sqrt(repart.^2 + impart.^2);
        if include_aphase || include_afreq
            aphase = atan2(impart,repart); end
        if include_aphase
            signal = utl_register_field(signal,'timeseries',[field '_aphase'],aphase); end
        if include_aamp
            signal = utl_register_field(signal,'timeseries',[field '_aamp'],aamp); end        
        if include_afreq
            % compute instantaneous frequency
            afreq = diff(aphase,[],2);
            flips = afreq<-pi;
            afreq(flips) = afreq(flips)+2*pi;
            signal = utl_register_field(signal,'timeseries',[field '_afreq'],afreq/(2*pi)*signal.srate);
        end
        if override_original
            signal.(field) = aamp; 
            % also override channel labels
            signal.chanlocs = state.cached_chanlocs;
            signal.nbchan = length(signal.chanlocs);
        end
    end
end

exp_endfun;


function locs = update_chanlocs(locs,nFreqs)
nChans = length(locs);
locs = repmat(locs,nFreqs,1);
[locs.labels] = celldeal(cellfun(@strcat,{locs.labels},vec(repmat(cellfun(@(f)sprintf('_%i',f),num2cell(1:nFreqs),'UniformOutput',false),nChans,1))','UniformOutput',false));
