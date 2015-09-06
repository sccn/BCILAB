function params = hlp_getDefaultArglist(prefix)
% params = hlp_getDefaultArglist(prefix)
% Returns a cell array of default arguments for finputcheck. The set of
% parameters depends on the specific prefix ('est','hlp','viz','pop')
%
% This function will be deprecated in the beta release


if ~ischar(prefix)
    help(hlp_getDefaultParams);
    error('input argument must be a prefix string (est,hlp,viz,pop)');
end

switch lower(prefix)
    case 'pre'
        params = {
            'resample',         'real'         []                0;              ...                     % new sampling rate
            'components',       'integer'      []               [];              ...
            'newtlims',         'real'         []               [];              ...                     % seconds
            'filter',           'cell'         {}               {};              ...                     % {low-edge hi-edge 'eegfilt'}
            'detrend',          'boolean'      []                0;              ...
            'normalize',        ''             {'ensemble','time',{'time','ensemble'},''}   '';    ...   % 'time' | 'ensemble' | 'both'
            'center',           'boolean'      []                0;              ...                     % remove the mean of each series
            'diff',             'integer'      [0 10]            0;              ...                     % number of times to difference
            'verb',             'integer'      [0 2]             1;              ...                     % verbosity level
            'DO_JACKET',        'boolean'      []                0;              ...                     % use jacket (GPU parallel)
            'newtrials',        'integer'      []               [];              ...                     % list of trials to include
            'badsegments',      'real'         []               [];              ...                     % K x 2 matrix of [lo hi] intervals (seconds) of data within ea. trial to set to nan
            'equalizetrials',   'boolean'      {0 1}             0;              ...                     % equalize the number of trials between two conditions
            };
    case 'est'
        params = {                                                               ...
            'algorithm',          ''           {'vieira-morf','arfit'}     'vieira-morf'; ...            % which algorithm to use for model fitting
            'winStartIdx',        'real'       []               [];              ...                     % vector of sample points (start of windows) at which to estimate windowed VAR model
            'morder',             'real'       []               [];              ...                     % VAR model order
            'winlen',             'real'       []               0.5;             ...                     % window length (sec)
            'winstep',            'real'       []               0.03;            ...                     % window step size (sec)
            'epochTimeLims',      'real'       []               [];              ...                     % time range to analyze (sec) where 0 = event time
            'prctWinToSample',    'real'       [0 100]          100;             ...                     % percent of time windows to randomly select
            'verb',               'real'       [0 2]            2;               ...                     % verbosity level (0=no output, 1=text, 2=gui)
            'timer'               'boolean'    []               0;               ...                     % estimate time required to fit model
            'normalize',          ''           {'ensemble','time',{'time','ensemble'},''}   '';    ...   % normalization 'time' | 'ensemble' | 'both'
            'icselector',         ''           {'sbc','aic','fpe','hq'}   {'sbc','aic'}; ...             % information criteria to use for selecting model order
            'detrend',            ''           {} [] ...
            };
    case 'hlp'
        params = {};
    case 'stat'
        params = {};
    case 'vis'
        params = {};
end
