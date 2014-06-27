function [signal,state] = flt_movavg(varargin)
% Filter a continuous data set by a moving-average filter.
% [Signal,State] = flt_fir(Signal, Length, State)
%
% In:
%   Signal        :   continuous data set to be filtered
%
%   Length        :   Filter length in samples.
%
%   State        :   previous filter state, as obtained by a previous execution of flt_fir on an
%                    immediately preceding data set (default: [])
%
% Out: 
%   Signal       :  filtered, continuous data set
%
%   State        :  state of the filter, after it got applied
%
% See also:
%   filter
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2012-12-11  

% flt_movavg_version<1.0> -- for the cache

if ~exp_beginfun('filter') return; end

% does not make sense on epoched data
declare_properties('name','MovingAverage', 'precedes','flt_resample', 'cannot_follow','set_makepos', 'independent_channels',false, 'independent_trials',true);

arg_define(varargin, ... 
    arg_norep({'signal','Signal'}), ...
    arg({'fltlength','Length'}, 4, uint32([1 1000000]), 'Length of the filter kernel.'), ...
    arg_nogui({'state','State'}));

% handle inputs
b = ones(1,fltlength)/fltlength;

for f = utl_timeseries_fields(signal) %#ok<NODEF>
    if ~isempty(signal.(f{1}))        
        % get the data (transposed)
        [X,dims] = spatialize_transpose(double(signal.(f{1})));        
        % if necessary prepend the signal with a mirror section of itself, to minimize start-up transients
        prepend = ~isfield(state,f{1});
        if prepend
            X = [repmat(2*X(1,:),length(b),1) - X(1+mod(((length(b)+1):-1:2)-1,size(X,1)),:); X]; %#ok<AGROW>
            state.(f{1}) = [];
        end
        % apply the filter
        [X,state.(f{1})] = filter(b,1,X,state.(f{1}),1);
        % check if we need to cut off a data segment that we previously prepended
        if prepend
            X(1:length(b),:) = []; end
        % write the data back
        signal.(f{1}) = unspatialize_transpose(X,dims);        
    end
end

if ~isfield(signal.etc,'filter_delay')
    signal.etc.filter_delay = 0; end
signal.etc.filter_delay = signal.etc.filter_delay + fltlength/(2*signal.srate);

exp_endfun;
