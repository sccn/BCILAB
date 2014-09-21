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
%               * 'timeaxis': the field is a time-axis field, which will be buffered during online 
%                             processing but will not be filtered
%               * 'samplingrate': the field indicates a sampling rate
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

% this is the list of fields that can be registered
permitted_fields = {'timeseries', 'timeaxis', 'samplingrate'};

if nargin < 2
    error('You need to pass in the field type.'); end
if nargin < 3
    error('You need to pass in the field name.'); end

if ischar(fieldtype) && any(strcmp(fieldtype,permitted_fields))
    if ~isfield(signal,'tracking') || ~isfield(signal.tracking,[fieldtype '_fields'])
        signal.tracking.([fieldtype '_fields']) = {fieldname};
    else
        signal.tracking.([fieldtype '_fields']) = [signal.tracking.([fieldtype '_fields']) {fieldname}];
    end
elseif ischar(fieldtype)
    error('Unrecognized field type passed in: %s',fieldtype);
else
    error('The second input needs to be type of field to append (e.g.,''timeseries''), but was: %s',hlp_tostring(fieldtype));
end

if nargin >= 4
    signal.(fieldname) = fieldvalue;
elseif ~isvarname(fieldname)
    error('The given Fieldname argument is not a valid field name: %s',hlp_tostring(fieldname));
end
