function signal = flt_fft(varargin)
% Apply an FFT to each epoch of an epoched signal (Example).
% Signal = flt_fft(Signal, LogPower)
%
% This is example code to transform a signal into the power domain, or log-power domain. A 
% fully-featured version of this is flt_fourier.
%
% In:
%   Signal :   Epoched data set to be processed
%
%   LogPower : whether to take the logarithm of the power (instead of the raw power) (default: false)
%
% Out:
%   Signal  :   processed data set
%
% Example:
%   % use default settings
%   eeg = flt_fft(eeg,true);
%
%   % use it with log-power enabled
%   eeg = flt_fft(eeg,true);
% 
% See also:
%   fft
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2011-01-19

% flt_fft_version<1.0> -- for the cache

if ~exp_beginfun('filter') return; end

% requires epoched data, works best on spatially filtered data
declare_properties('name','EpochedFFT', 'depends','set_makepos', 'follows',{'flt_project','flt_window'}, 'independent_channels',true, 'independent_trials',true);

% declare arguments
arg_define(varargin,...
    arg_norep({'signal','Signal'}), ...
    arg({'do_logpower','LogPower'}, false, [], 'Compute log-power. Taking the logarithm of the power in each frequency band is easier to handle for simple statistical classifiers, such as LDA.'));

% apply FFT and cut mirror half of the resulting samples
tmp = fft(signal.data,[],2);
tmp = tmp(:,1:signal.pnts/2,:);

% take signal power or log(power)
if do_logpower 
    signal.data = log(abs(tmp));
else
    signal.data = abs(tmp);
end

exp_endfun;