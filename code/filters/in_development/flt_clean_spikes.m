function signal = flt_clean_spikes(varargin)
% Set outlier samples in the data to zero.
% Signal = flt_clean_spikes(Signal,Quantile)
%
% Note: This basically needs to be reimplemented with robust moving-window statistics.
% It appears that this requires a mex file that does partial sorting...
%
% In:
%   Signal   : a continous data set
%
%   Quantile : quantile of the data that should be retained (data higher than that is set to zero).
%              (default: 0.9)
% 
% Out:
%   Signal   : data set with outliers set to zero
% 
% Examples:
%   % use the defaults
%   eeg = flt_clean_spikes(eeg)
%
%   % use a more aggressive threshold
%   eeg = flt_clean_spikes(eeg,0.8);
% 
% See also:
%   flt_clean_peaks
%
% TODO:
%   Replace by fast running median/mad based spike detection
%
%                                  Laura Frolich, Swartz Center for Computational Neuroscience, UCSD
%                                  2010-09-22

% flt_clean_spikes_version<0.1> -- for the cache

if ~exp_beginfun('offline') return; end

arg_define(varargin,...
    arg_norep({'signal','Signal'}), ...
    arg({'qq','Quantile'}, 0.995, [0 1], 'Quantile of the data to retain. Data higher than that is set to zero.'));

for i = size(signal.data,1) %#ok<*NODEF>
    signal.data(i,abs(signal.data(i,:)) > quantile(abs(signal.data(i,:)),qq)) = 0; end

exp_endfun;
