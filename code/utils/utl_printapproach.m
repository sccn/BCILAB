function string = utl_printapproach(app,strip_direct,indent,indent_incr)
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
% Out:
%   String : a string representation of the approach, for use in scripts
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2013-10-22

% check inputs
if nargin < 2
    strip_direct = true; end
if nargin < 3
    indent = 0; end
if nargin < 4
    indent_incr = 4; end

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
    error('Failed to instantiate the paradigm %s with error: %s.',char(paradigm),e.message);
end
func = @instance.calibrate;

% report both the defaults of the paradigm and 
% the current settings in form of argument specifications
try
    defaults = clean_fields(arg_report('rich',func));
catch e
    hlp_handleerror(e);
    error('Failed to report default arguments of the given paradigm''s calibrate() method with error: %s',e.message);
end

try
    settings = clean_fields(arg_report('lean',func,parameters));
catch e
    hlp_handleerror(e);
    error('Failed to process parameters of the given paradigm''s calibrate() method with error: %s',e.message);
end

% get the difference as cell array of human-readable name-value pairs
try
    specdiff = arg_diff(defaults,settings);
catch e
    hlp_handleerror(e);
    error('Failed to calculate difference between default parameters and overridden parameters with error: %s',e.message);
end

try
    difference = arg_tovals(specdiff,[],'HumanReadableCell',false);
catch e
    hlp_handleerror(e);
    error('Failed to convert argument difference to nested cell-array form with error: %s',e.message);    
end

% pre-pend the paradigm choice 
paradigm_name = char(paradigm);
difference = [{'arg_selection',paradigm_name(9:end)} difference];

% and convert to string
try
    string = arg_tostring(difference,strip_direct,indent,indent_incr);
catch e
    hlp_handleerror(e);
    error('Failed to convert argument difference to string with error: %s',e.message);        
end


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
