function utl_check_fields(x,required_fields,argname,typename)
% Check whether a given value is a scalar struct that has the required fields.
%
% This function ensures that the given value is a 1x1 struct that has the required fields and
% produces and informative error message if not.
% 
% In:
%   Value : value to check
%
%   RequiredFields : fields that need to be present
%
%   ParameterName : name of the parameter in question
%
%   TypeName : type that the parameter should have (e.g., "signal")
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2014-03-02

if ~isstruct(x)
    error('The given %s parameter must be a struct, but was: %s',argname,hlp_tostring(x,10000)); end
if ~isscalar(x)
    error('The given %s parameter cannot be a struct array (must be a single struct), but was: %s',argname,hlp_tostring(x,10000)); end
if ~all(isfield(x,required_fields))
    error('The given %s parameter is lacking the required fields %s; it is not a valid %s (its fields were: %s).',argname,hlp_tostring(setdiff(required_fields,fieldnames(x))),typename,hlp_tostring(fieldnames(x))); end
