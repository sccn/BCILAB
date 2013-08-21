function signal = utl_register_field(signal,fieldtype,fieldname,fieldvalue)
% Register a custom field in a signal.
% Signal = utl_register_field(Signal,Fieldtype,Fieldname,Value)
%
% Some processing functions are able to treat custom fields in an appropriate manner (for example,
% extra time-series fields in a signal, which would have to be buffered during online processing).
% This function allows to register such fields as having a certain type.
%
% In:
%   Signal : EEGLAB data set struct
%
%   Fieldtype : type of the field to register; possible options include:
%               * 'timeseries': the field is a time-series field, and will be filtered by certain
%                               filters, and will be buffered during online processing
%
%   Fieldname : name of a field that shall be registered in the struct
%               for use by certain processing functions
%
%   FieldValue : optionally the initial value to assign to the field
%
% Out:
%   Signal : updated signal
%
% See also:
%   utl_timeseries_fields
%
%                                       Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                       2013-08-16

if nargin < 2
    error('You need to pass in the field type.'); end
if nargin < 3
    error('You need to pass in the field name.'); end

switch fieldtype
    case 'timeseries'
        if ~isfield(signal,'tracking') || ~isfield(signal.tracking,'timeseries_fields')
            signal.tracking.timeseries_fields = {fieldname};
        else
            signal.tracking.timeseries_fields = [signal.tracking.timeseries_fields {fieldname}];
        end
    otherwise
        error(['Unrecognized field type passed in: ' hlp_tostring(fieldtype)]);
end

if nargin >= 4
    signal.(fieldname) = fieldvalue; end
