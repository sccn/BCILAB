function signal = set_picktimes(varargin)
% Pick time intervals from an epoched data set.
% Signal = set_picktimes(Signal, TimeRanges)
%
% In:
%   Signal     : epoched data set for which time ranges should be averaged
%
%   TimeRanges : 2d matrix specifying time ranges which should be selected, averaged,
%                 and emitted as new samples into the resulting data set
%                 specified as [Range1-Begin Range1-End; Range2-Begin Range2-End, ...]
%
% Out:
%   Signal : newly created epoched data set, whose samples are averages over the supplied time ranges
%
% Examples:
%   % for an epoched data set, average each signal segment in the time intervals (relative to the time-
%   % locking event) within 0.25s to 0.5s, 0.5s to 0.6s and 0.8s to 1.5s, and output a new epoched
%   % data set which has only 3 samples per epoch and channel, containing the averaged values
%   eeg = set_picktimes(eeg,[0.25 0.5; 0.5 0.6; 0.7 1.5])
%
%   % as before, but passing all arguments by name
%   eeg = set_picktimes('Signal',eeg, 'TimeRanges',[0.25 0.5; 0.5 0.6; 0.7 1.5])
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-01

% set_picktimes_version<1.0> -- for the cache

if ~exp_beginfun('editing') return; end

declare_properties('independent_channels',true,'independent_trials',true);

arg_define(varargin, ...
    arg_norep({'signal','Signal'}), ...
    arg({'ranges','TimeRanges'},[],[],'Selected time ranges. Array of the form [start end; start end; start end; ...].'));
    
signal.data = utl_picktimes(signal.data, (ranges-signal.xmin) * signal.srate);

exp_endfun;
