function options = adigatorOptions(varargin)
% function options = adigatorOptions(field1,value1,field2,value2,...)
% OPTIONS:
%    auxdata:  1 - auxiliary inputs will always have the same sparsity
%                  pattern but may change numeric values 
%              0 - auxiliary inputs will always have the same numeric
%                  values (required if inputs are a size to be looped on or
%                  a reference index, etc.) (default)
%       echo:  1 - echo to screen the transformation progress (default)
%              0 - dont echo
%     unroll:  1 - unroll loops and sub functions in derivative program
%              0 - keep loops and sub functions rolled in derivative 
%                  (default)
%   comments:  1 - print comments to derivative file giving the lines of
%                  user code which correspond to the printed statements
%                  (default)
%              0 - dont print comments
%  overwrite:  1 - if the user supplies a derivative file name which
%                  corresponds to an already existing file, then setting
%                  this option will overwrite the existing file.
%              0 - if a file already exists with the same name as the given
%                  derivative file name, then it will not overwrite and
%                  error out instead (default)
%
% If desired, the defaults of each option may be changed by editing this
% file.
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
%
% See also adigator

% Set Defaults
options.auxdata   = 0;
options.echo      = 1;
options.unroll    = 0;
options.comments  = 1;
options.overwrite = 0;
if nargin/2 ~= floor(nargin/2)
  error('Inputs to adigatorOptions must come in field/value pairs')
end

% Set user wanted options
for i = 1:nargin/2
  field = varargin{2*(i-1)+1};
  value = varargin{2*i};
  if value ~= 0 && value ~= 1
    error('Value must be binary 1 or 0')
  end
  switch lower(field)
    case {'auxdata','echo','unroll','comments','overwrite'}
      options.(field) = value;
    otherwise
      error(['Invalid option field: ',field])
  end
end