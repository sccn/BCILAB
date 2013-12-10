function signal = flt_sigmoid(varargin)
% Sigmoid transform of data
% Tim Mullen, SCCN/INC

if ~exp_beginfun('filter') return; end

declare_properties('name', 'SigmoidTransform','cannot_follow','set_makepos', 'independent_trials',true, 'independent_channels',true);

g = arg_define(varargin, ...
        arg_norep({'signal','Signal'}), ...
        arg({'amp','Amplitude'},1,[],'Amplitude'), ...
        arg({'phase','Phase'},1,[],'Phase. X-shifting'), ...
        arg({'slope','Slope'},0.7,[],'Slope. Larger means narrower range for squashing.'), ...
        arg({'x_shift'},0,[],'X translation'), ...
        arg({'y_shift'},0,[],'Y translation'));

signal = g.signal;
% sigmoid transform parameters: [amp, phase, slope, x_shift, y_shift]
for f=utl_timeseries_fields(signal)
    if ~isempty(signal.(f{1}))
        signal.(f{1}) = sigmoid([g.amp,g.phase,g.slope,g.x_shift,g.y_shift],signal.(f{1})); end
end

exp_endfun;