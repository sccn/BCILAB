function signal = flt_divideBySum(varargin)
% Normalize data by the sum across channels.
% Tim Mullen, SCCN/INC

if ~exp_beginfun('filter') return; end

declare_properties('name', 'DivideBySum', 'cannot_follow','set_makepos', 'independent_trials',true, 'independent_channels',true);

g = arg_define(varargin, ...
        arg_norep({'signal','Signal'}));

signal = g.signal;

% normalization by sum
signal.data = signal.data./sum(signal.data);

exp_endfun;