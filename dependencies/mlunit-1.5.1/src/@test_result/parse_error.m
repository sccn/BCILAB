function [error, stacktrace] = parse_error(self, error, stacktrace) %#ok<INUSL>
%test_result/parse_error parses special errors to extract further
%information for the stacktrace.
%
%  Example: See test_result/add_error_with_stack.m.
%
%  See also TEST_RESULT/ADD_ERROR_WITH_STACK.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: parse_error.m 271 2007-04-06 13:57:40Z thomi $

if (~isempty(strfind(error, 'Unbalanced or misused parentheses or brackets.')) || ...
    ~isempty(strfind(error, 'Unbalanced or unexpected parenthesis or bracket.')))
    [tokens] = regexp(error, 'Error:.*File:\ ([\w\ \.,$&/\\:@]*.m)\ Line: (\w*)\ Column: (\w*).*', 'tokens', 'once');
    if (length(tokens) == 3)
        fullname = which(char(tokens(1)));
        if (~isempty(fullname))
            stacktrace = sprintf('\n  In %s at line %s%s', ...
                fullname, char(tokens(2)), ...
                stacktrace);
        else
            stacktrace = sprintf('\n  In %s at line %s%s', ...
                char(tokens(1)), char(tokens(2)), ...
                stacktrace);
        end;
        error = 'Unbalanced or misused parentheses or brackets.';
    end;
else
    [tokens] = regexp(error, 'Error using ==> <a href.*>(.*)</a>\n(.*)', 'tokens', 'once');
    if (length(tokens) == 2)
        error = char(tokens(2));
    end;
end;