function [signal,newstate] = flt_movavg(varargin)
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
    arg({'fltlength','Length'}, 4, [], 'Length of the filter kernel.'), ...
    arg_nogui({'state','State'}));

% handle inputs
b = ones(1,fltlength)/fltlength;

for fld = utl_timeseries_fields(signal)
    field = fld{1};
    if ~isempty(signal.(field))        
        % get the data (transposed)
        sig = double(signal.(field))';
        
        % if necessary prepend the signal with a mirror section of itself, to minimize start-up transients
        if ~isfield(state,field) || isempty(state.(field))
            sig = [repmat(2*sig(1,:),length(b),1) - sig((length(b)+1):-1:2,:); sig];
            state.(field) = [];
        end
        
        % apply the filter
        [sig,newstate.(field)] = filter(b,1,sig,state.(field),1);
        
        % check if we need to cut off a data segment that we previously prepended
        if isempty(state.(field))
            sig(1:length(b),:) = []; end
        
        % write the data back
        signal.(field) = sig';
    end
end

exp_endfun;
