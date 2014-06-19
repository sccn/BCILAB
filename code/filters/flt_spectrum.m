function signal = flt_spectrum(varargin)
% Select a frequency portion of the data in an epoched data set.
% Signal = flt_spectrum(Signal, FrequencySpec)
%
% flt_select can be used to implement spectral filtering of the data. This is done by specifying one
% or more frequency windows (or a function of frequency), to which the data should be restricted.
% Spectral filtering is critical, on one hand, to separate oscillatory processes of interest (e.g.
% thalamo-cortical loops) from other signal contents (e.g., muscle artifacts, slow cortical
% potentials) and on the other hand to isolate slow cortical potentials (e.g. some event-related
% negativity) from higher-frequency signal content.
%
% In:
%   Signal  : epoched data set
%
%   FrequencySpec: frequency-domain selection, can be specified as one of the following
%                   * [low high] or [low high; low high; low high; ...] for a brickwall filter or
%                     sum of brickwall filters (numbers in Hz)
%                   * [low_begin low_end high_begin high_end; ...] for a linearly sloped filter or
%                     sum of sloped filters (numbers in Hz)
%                   * function_handle for an arbitrary function of frequency (in Hz)
%                   low and high frequencies are generally inclusive
%
% Out:
%   Signal      :   subset of the data set
%
% Examples:
%   % apply a 10-15Hz band-pass filter
%   eeg = flt_spectrum(eeg,[10 15])
%
%   % apply a filter which retains 10-15 Hz and 20-25 Hz ranges
%   eeg = flt_spectrum(eeg,[10 15; 20 25]);
%
%   % apply a filter which retains 10-15 Hz but with linear falloffs at both edges
%   eeg = flt_spectrum(eeg,[7 10 15 17]);
%
%   % apply a filter which retains 10-15 Hz and 25-30Hz but with linear falloffs at the edges
%   eeg = flt_spectrum(eeg,[7 10 15 17; 22 25 30 33]);
%
%   % apply a filter whose frequency response is specified by a function
%   % (here: an unnamed function implementing 7-30 Hz brickwall)
%   eeg = flt_spectrum(eeg,@(f)f>7&f<30)
%
%   % apply a filter whose frequency response is specified by a function
%   % (here: a smooth sine-based function with peaks at 12 and 25 Hz)
%   func = @(f)(f>7&f<30).*(1-cos((f-(7+30)/2)/(7-30)*pi*4));
%   eeg = flt_spectrum(eeg,func)
%
%   % retain the original data (i.e., do nothing)
%   eeg = flt_spectrum(eeg,[])
%
% References:
%   [1] Oppenheim, A.V., and R.W. Schafer, "Discrete-Time Signal Processing", 
%       Prentice-Hall, 1989.
%
% See also:
%   fft, ifft
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-06-28

% flt_spectrum_version<1.0> -- for the cache

if ~exp_beginfun('filter') return; end

declare_properties('name','SpectralSelection', 'depends','set_makepos', 'follows',{'flt_window','flt_rmbase'}, 'independent_channels',true, 'independent_trials',true);

arg_define(varargin, ... 
    arg_norep({'signal','Signal'}), ...
    arg({'freq','FrequencySpecification','Frequencies'}, [], [], ['Frequency-domain selection. Can be specified either as [low high] or [low high; low high; low high; ...] for a brickwall filter or sum of brickwall filters, ' ...
                                                    'or alternatively as [low_begin low_end high_begin high_end; ...] for a linearly sloped filter or sum of sloped filters (frequencies in Hz).']));

if ~isempty(freq) %#ok<*NODEF>    
    % apply frequency-domain selection to all known time-series fields in
    % the data
    for fld = utl_timeseries_fields(signal)
        field = fld{1};
        if ~isempty(signal.(field))            
            T = size(signal.(field),2);
            % almost correct frequency vector (fixed further below)
            frqs = (0:T-1)*signal.srate/T;
            if isnumeric(freq) 
                if size(freq,2) == 2
                    % brickwall filters
                    flt = zeros(1,T);
                    for k=1:size(freq,1)
                        f = freq(k,:);
                        flt = max(flt,(frqs>=f(1) & frqs<=f(2))); end
                elseif size(freq,2) == 4
                    % sloped filters
                    flt = zeros(1,T);
                    for k=1:size(freq,1)
                        f = freq(k,:);
                        flt = max(flt,min(max(0,min(1,(frqs-f(1)) / (f(2)-f(1)))), max(0,min(1, 1-(frqs-f(3)) / (f(4)-f(3))))));
                    end
                else
                    error('Unsupported format for the FrequencySpecification: %s',hlp_tostring(freq,10000));
                end
            else
                % function handle filters
                flt = freq(frqs);
            end
            % fix the filter beyond the Nyquist freq
            flt = (flt(:) .* (frqs(:)<(signal.srate/2)))';
            % filter the data
            signal.(field) = 2*real(ifft(bsxfun(@times,fft(signal.(field),[],2),flt),[],2));
        end
    end
end

exp_endfun;
