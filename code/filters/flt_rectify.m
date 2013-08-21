function signal = flt_rectify(varargin)
% Rectify a signal. Result is the magnitude (absval) of the signal
% If squareData=true, signal is squared after absvalue.
% Tim Mullen, SCCN/INC

if ~exp_beginfun('filter') return; end

declare_properties('name', 'Rectify', 'precedes','set_makepos', 'independent_trials',true, 'independent_channels',true);

arg_define(varargin, ...
        arg_norep({'signal','Signal'}), ...
        arg({'squareData','SquareData'},false,[],'Return abs(x).^2. Otherwise only abs(x) is returned.'));
     
if ~isreal(signal.data) || ~squareData
    signal.data = abs(signal.data); end

if squareData
    signal.data = signal.data.^2; end

exp_endfun;