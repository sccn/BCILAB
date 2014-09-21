function fields = utl_timeseries_fields(signal, include_timeaxis)
% Get the time-series fields of te given signal.
% function Fields = utl_timeseries_fields(Signal, IncludeTimeaxis)
%
% This function returns the field names in the given signal that are 
% carrying time-series information (and which therefore should be filtered,
% buffered, etc.). The result is a combination of a set of hard-coded fields
% and the names listed in the field .tracking.timeseries_fields. Only fields
% that are present in the signal are returned.
%
% In:
%   Signal : a signal for which time-series field names shall be looked up
%
%   IncludeTimeaxis : whether to also return the time axis fields (which shall not be filtered but
%                     buffered) (e.g., .times)
%
% Out:
%   Fields : Cell array of field names (row vector). It is assumed that the second dimension of
%            these fields is the time axis along which buffering and/or filtering should happen; any
%            number of other dimensions may be present for any field returned by this function.
%
% See also:
%   utl_register_field
%
%                                       Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                       2013-08-16

if include_timeaxis
    fields = utl_registered_fields(signal,{'timeseries','timeaxis'});
else
    fields = utl_registered_fields(signal,'timeseries');
end
