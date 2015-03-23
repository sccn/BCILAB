function result = arg_report(type,func,args)
% Report information of a certain Type about the given Function.
% Result = arg_report(Type,Function,Arguments)
%
% Functions that declare their arguments via arg_define() make their parameter specification
% accessible to outside functions. This can be used to display auto-generated settings dialogs, to
% record function calls, and so on. Some other functions can enable reports, as well (like
% declare_properties() and expose_handles()).
%
% In:
%   Type : Type of information to report. The currently defined report types are as follows:
%
%          'vals' : Report the values of the function's arguments as a struct, possibly with
%                   sub-structs, after assignment of Arguments. (provided by arg_define)
%          'nvps' : Report the values of the function's arguments as cell array of name-value pairs.
%                   (provided by arg_define)
%          'lean' : Report a lean specification of the function's arguments as a struct array, with
%                   fields as in arg_specifier, including help text, etc. This is after assignment
%                   of Arguments. (provided by arg_define)
%          'rich' : Report a rich specification of the function's arguments as a struct array, with
%                   fields as in arg_specifier. In addition to lean this includes information about
%                   alternative (non-default) settings, for use in GUI generation. (provided by
%                   arg_define)
%
%          'properties' : Report declared properties of the function, if any. Arguments must be empty.
%                         (provided by declare_properties)
%
%          'handle': Report function handles to scoped and nested functions within the Function. The
%                    names of those functions are given as a cell array in Arguments, and the result
%                    is a cell array of associated function handles (or one handle if only one
%                    function was requested).
%
%   Function : Handle to a function that supports the given report type.
%
%   Arguments : Cell array of inputs to be passed to the function.
%
% Out:
%   Result : the reported data.
%
% Examples:
%   % for a function call with some arguments assigned, obtain a struct with all parameter
%   % names and values, including defaults
%   params = arg_report('vals',@myfunction,{4,10,true,'option1','xxx','option5',10})
%
%   % obtain a specification of all function arguments, with defaults, help text, type, shape, and other
%   % meta-data (with a subset of settings customized according to arguments)
%   spec = arg_report('rich',@myfunction,myarguments)
%
%   % obtain a report of properties of the function (declared via declared_properties() within the
%   % function)
%   props = arg_report('properties',@myfunction)
%
% See also:
%   arg_define, declare_properties, expose_handles
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-09-24

% Copyright (C) Christian Kothe, SCCN, 2010, christian@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License as published by the Free Software Foundation; either version 2 of the
% License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
% even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License along with this program; if not,
% write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
% USA

% check inputs
if ischar(func)
    func = str2func(func); end
if ~isa(func,'function_handle')
    error('The given Function argument must be a function handle.'); end
if ~ischar(type)
    error('The given Type argument must be a string.'); end
if nargin < 3
    if strcmp(type,'properties')
        % for the properties report we implicitly pad the arguments with blanks to allow functions
        % to define up to this number of traditional non-varargin arguments prior to the varargin
        args(1:15) = {[]};
    else
        args = {};
    end
elseif isequal(args,[])
    args = {};
elseif ~iscell(args)
    error('Arguments must be a cell array, if given.'); 
end

% disable the expression features of functions called by arg_report
persistent handle_expressions;
if isempty(handle_expressions)
    handle_expressions = exist('exp_eval','file'); end

result = {};
try
    % evaluate the funtion with two special arguments appended
    if handle_expressions
        hlp_scope({'disable_expressions',1},func,args{:},'__arg_report__',type);
    else
        func(args{:},'__arg_report__',type);
    end
catch report
    if strcmp(report.identifier,'BCILAB:arg:report_args')
        global tracking; %#ok<TLEV>
        % get the ticket of the report
        ticket = sscanf(report.message(end-4:end),'%u');
        % read out the payload and return the ticket
        result = tracking.arg_sys.reports{ticket};
        tracking.arg_sys.tickets.addLast(ticket);
    elseif ~strcmp(type,'properties')        
        % got a genuine error
        settings = dbstatus;
        if any(strcmp({settings.cond},'error')) && ~hlp_resolve('disable_dbstop_if_error_msg',false)
            % rerun the code to propagate it properly (we're not using rethrow here to not confuse dbstop if error)
            disp_once('arg_report: caught error while in "dbstop if error" mode; re-running offending code to get break point...\n  error traceback was: %s\n',hlp_handleerror(report));
            if handle_expressions
                hlp_scope({'disable_expressions',1},func,args{:},'__arg_report__',type);
            else
                func(args{:},'__arg_report__',type);
            end
        else
            rethrow(report);
        end
    end
end
