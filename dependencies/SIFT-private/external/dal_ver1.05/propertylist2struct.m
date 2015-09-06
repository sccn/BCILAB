function opt = propertylist2struct(varargin)
% PROPERTYLIST2STRUCT - Make options structure from parameter/value list
%
%   OPT = PROPERTYLIST2STRUCT('param1', VALUE1, 'param2', VALUE2, ...)
%   Generate a structure OPT with fields 'param1' set to value VALUE1, field
%   'param2' set to value VALUE2, and so forth.
%
%   OPT has an additional field 'isPropertyStruct' that is meant to identify
%   OPT is a structure containing options data. Only in the case of missing
%   input arguments, no such identification field is written, that is,
%   PROPERTYLIST2STRUCT() returns [].
%
%   OPT2 = PROPERTYLIST2STRUCT(OPT1, 'param', VALUE, ...) takes the options
%   structure OPT1 and adds new fields 'param' with according VALUE.
%
%   See also SET_DEFAULTS
%

% Copyright Fraunhofer FIRST.IDA (2004)

if nargin==0,
  % Return an empty struct without identification tag
  opt= [];
  return;
end

if isstruct(varargin{1}) | isempty(varargin{1}),
  % First input argument is already a structure: Start with that, write
  % the additional fields
  opt= varargin{1};
  iListOffset= 1;
else
  % First argument is not a structure: Assume this is the start of the
  % parameter/value list
  opt = [];
  iListOffset = 0;
end
% Write the identification field. ID field contains a 'version number' of
% how parameters are passed.
opt.isPropertyStruct = 1;

nFields= (nargin-iListOffset)/2;
if nFields~=round(nFields),
  error('Invalid parameter/value list');
end

for ff= 1:nFields,
  fld = varargin{iListOffset+2*ff-1};
  if ~ischar(fld),
    error(sprintf('String required on position %i of the parameter/value list', ...
                  iListOffset+2*ff-1));
  end
  prp= varargin{iListOffset+2*ff};
  opt= setfield(opt, fld, prp);
end




