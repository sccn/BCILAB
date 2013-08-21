function signal = set_selinterval(varargin)
% Selects a time interval from a data set.
% Signal = set_selinterval(Signal,Intervals,Unit)
%
% In:
%   Signal    :   continuous EEGLAB data set from which an interval should be selected
%
%   Intervals :   selection interval(s), [start end; start end; start end; ...]
%
%   Unit      :   Unit of measurement for the interval. Either 'seconds', 'samples', or 'fraction' 
%                 (default: seconds)
%
% Out:
%   Signal    :   data set restricted to the selected range
%
% Examples:
%   % for a continuous data set, retain only the data within 50s and 200s, as well as 1200s and 1500s
%   eeg = set_selinterval(eeg,[50 200; 1200 1500])
%
%   % for a continuous data set, retain only the data within the 1000's and the 10000's sample
%   eeg = set_selinterval(eeg,[1000 10000],'samples')
%
%   % for a continuous data set, retain only the last half of the data
%   eeg = set_selinterval(eeg,[0.5 1],'fraction')
%
%   % as before, but pass the arguments by name
%   eeg = set_selinterval('Signal',eeg, 'Intervals',[0.5 1], 'Unit','fraction')
%
% See also:
%   set_selepos, set_partition
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-01

% set_selinterval_version<1.0> -- for the cache

if ~exp_beginfun('editing') return; end

declare_properties('independent_channels',true,'independent_trials',false);

arg_define(varargin, ...
    arg_norep({'signal','Signal'}), ...
    arg({'intervals','Intervals'},[],[], 'Selection intervals. Array of the form [start end; start end; start end; ...].'), ...
    arg({'selunit','Unit','unit'},'seconds',{'seconds','samples','fraction'}, 'Unit of measurement.'));

if strcmp(selunit,'fraction')
    signal = pop_select(signal,'point',intervals*signal.pnts);
else
    signal = pop_select(signal,hlp_rewrite(selunit,'seconds','time','samples','point'),intervals);
end

exp_endfun;
