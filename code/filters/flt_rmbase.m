function signal = flt_rmbase(varargin)
% Subtract a baseline from an data set, computed over the given baseline window.
% Signal = flt_rmbase(Signal, Window)
%
% Baseline correction is a commonly used method in the analysis of Event-Related Potentials, to
% factor out irrelevant variations in a signal's baseline. It is applied to an epoched data set, by
% specifying a sub-window that (typically) lies before the phenomena of interest, where the mean
% signal value in that window is subtracted from the entire signal (per channel or component).
% Baseline considerations are relevant for BCIs operating on slow cortical potentials, and can be
% implemented in a variety of other ways, too. One is to use a highpass filter (flt_iir, flt_fir,
% flt_select) to subtract low-frequency drifts in the signal. Another one is to add features from
% the signal which measure the baseline, for example one or a collection of windows of various
% lengths prior to the phenomenon of interest (in this case, the baseline correction is done by the
% machine learning algorithm that operates on these features). Paradigms which assign different
% weights for every time point in an epoch typically do baseline correction implicitly.
%
% In:
%   Signal          :   epoched data set to be processed
%
%   BaselineWindow  :   baseline window in seconds, e.g. [-0.5 -0.3]
%
% Out: 
%   Signal  :   baseline-corrected data set
%
% Examples:
%   % remove the baseline of a continuous or epoched data set
%   eeg = flt_rmbase(eeg)
%
%   % in an epoched signal, subtract the average of the signal from 250ms before the time-locking 
%   % event to 100ms after the time-locking event
%   eeg = flt_rmbase(eeg,[-0.25 0.1])
%
%   % pass the arguments by name
%   eeg = flt_rmbase('Signal'eeg,'BaselineWindow',[-0.25 0.1])
%
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-03-28

% flt_rmbase_version<1.0> -- for the cache

if ~exp_beginfun('filter') return; end

% only useful for epoch baselines
declare_properties('name','BaselineRemoval', 'depends','set_makepos', 'precedes','flt_window', 'independent_channels',true, 'independent_trials',true);

arg_define(varargin, ...
    arg_norep({'signal','Signal'}), ...
    arg({'wnd','BaselineWindow'}, [], [], 'Baseline window in seconds.','shape','row'));

if isempty(wnd)  %#ok<*NODEF>
    wnd = [-Inf Inf]; end

for f = utl_timeseries_fields(signal)
    signal.(f{1}) = double(signal.(f{1}));
    wnd = round(max(1,min(size(signal.(f{1}),2),(wnd-signal.xmin)*signal.srate+1)));
    signal.(f{1}) = signal.(f{1}) - repmat(mean(signal.(f{1})(:,wnd(1):wnd(2),:),2),[1,size(signal.(f{1}),2),1]);
end

exp_endfun;
