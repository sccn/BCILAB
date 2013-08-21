function defaults = exp_settings(setting)
% Define custom settings (default options) for exp_beginfun.
% Defaults = exp_settings(Setting)
%
% In:
%   Setting : string that identifies the setting; first argument to exp_beginfun
%
% Out:
%   Defaults : cell array of name-value pairs, denoting custom defaults to exp_beginfun's predefined
%              or custom attributes
%
% See also:
%   exp_beginfun
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-14


% define custom settings for use in exp_begindef
switch(setting)
    case 'filter'
        defaults = {'argsteps',{@utl_check_dataset}, 'poststeps',{@utl_add_online}, 'set_online','reproduce'};
    case 'editing'
        defaults = {'argsteps',{@utl_check_dataset}, 'poststeps',{@utl_add_online}, 'set_online','passthrough'};
    case 'offline'
        defaults = {'argsteps',{@utl_check_dataset}, 'poststeps',{@utl_add_online}, 'set_online','inapplicable'};
    otherwise
        error('Unsupported setting.');
end
