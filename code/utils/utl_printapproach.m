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
else
    paradigm = app.paradigm;
    parameters = app.parameters;
end

% get a handle to the paradigm's calibrate() function
instance = eval(paradigm);
func = @instance.calibrate;

% report both the defaults of the paradigm and 
% the current settings in form of argument specifications
defaults = clean_fields(arg_report('rich',func));
settings = clean_fields(arg_report('lean',func,parameters));

% get the difference as cell array of human-readable name-value pairs
specdiff = arg_diff(defaults,settings);
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
catch
end
