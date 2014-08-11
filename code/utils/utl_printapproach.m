function string = utl_printapproach(app,strip_direct,indent,indent_incr,subset)
% Convert an approach to a string representation
% String = utl_printapproach(Approach)
%
% In:
%   Approach : a BCI approach, either designed in the GUI or constructed in a script
%   
%   StripDirect : strip arg_direct flags (default: true)
%
%   Indent : initial indent (default: 0)
%
%   IndentIncrement : indentation increment (default: 4)
%
%   Subset : subset of parameters to display; if 'all' then all parameters are printed, 
%            and if 'diff' then only those parameters are printed that differ from the defaults
%            of the respective paradigm (default: 'diff')
%
% Out:
%   String : a string representation of the approach, for use in scripts
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2013-10-22

% check inputs
if nargin < 2 || isempty(strip_direct)
    strip_direct = true; end
if nargin < 3 || isempty(indent)
    indent = 0; end
if nargin < 4 || isempty(indent_incr)
    indent_incr = 4; end
if nargin < 5 || isempty(subset)
    subset = 'diff'; end
if ~isa(strip_direct,'logical')
    error('The StripDirect argument must be a logical/boolean value.'); end

% get required approach properties
if ischar(app)
    paradigm = ['Paradigm' app];
    parameters = {};
elseif iscell(app)
    paradigm = ['Paradigm' app{1}];
    parameters = app(2:end);
elseif all(isfield(app,{'paradigm','parameters'}))
    paradigm = app.paradigm;
    parameters = app.parameters;
else
    error('The given data structure is not an approach.');
end

% get a handle to the paradigm's calibrate() function
try
    instance = eval(paradigm);
catch e
    if ~exist(char(paradigm),'file')
        error('A BCI paradigm with name %s does not exist.',char(paradigm));
    else
        error('Failed to instantiate the paradigm named %s with error: %s; this is likely an error in the Paradigm''s code.',char(paradigm),e.message);
    end
end
func = @instance.calibrate;

% report the defaults of the paradigm
defaults = clean_fields(arg_report('rich',func));

% report the current settings of the paradigm in form of an argument specification
settings = clean_fields(arg_report('lean',func,parameters));

% get the difference between the defaults and settings
if strcmp(subset,'diff')
    specdiff = arg_diff(defaults,settings);
elseif strcmp(subset,'all')
    specdiff = settings;
else
    error('Unsupported subset: %s',hlp_tostring(subset,100));
end

% convert to nested cell arrays of human-readable name-value pairs
difference = arg_tovals(specdiff,[],'HumanReadableCell',false);

% pre-pend the paradigm choice 
paradigm_name = char(paradigm);
difference = [{'arg_selection',paradigm_name(9:end)} difference];

% and convert to string
string = arg_tostring(difference,strip_direct,indent,indent_incr);


% clean fields of x, by removing all arg_direct instances and 
% all skippable fields
function x = clean_fields(x)
x(strcmp({x.first_name},'arg_direct') | [x.skippable]) = [];
try
    children = {x.children};
    empty_children = cellfun('isempty',children);
    [x(~empty_children).children] = celldeal(cellfun(@clean_fields,children(~empty_children),'UniformOutput',false));    
catch %#ok<CTCH>
end
