function shortened_text = shorten_error_text(self, error_text)
%gui_test_runner/set_shorten shortens the directory paths of
%an error message. This is helpful for a small width of the gui window.
%
%  Example
%  =======
%  The method is internal to the mlUnit framework and should not be called
%  directly.
%
%  See also GUI_TEST_RUNNER, GUI_TEST_RUNNER/GUI,
%           GUI_TEST_RUNNER/SET_SHORTEN. 

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: shorten_error_text.m 252 2007-01-27 21:02:36Z thomi $

if (self.shorten == 0)
    shortened_text = error_text;
    return;
end;

error_lines = strread(error_text, '%s', 'delimiter', '\n'); %regexp(error_text, '(.*)', 'tokens', 'dotexceptnewline');
if (isempty(error_lines))
    shortened_text = error_text;
    return;
end;

shortened_text = '';
for i = 1:length(error_lines)
    line = error_lines{i};
    [tokens] = regexp(char(line), ['^', get_line_expression(self)], 'tokens', 'once');
    if (length(tokens) == 2)
        token = tokens{1};
        token = strrep(token, '\', '/');
        str = cell(1);
        [str{1}, rem] = strtok(token, '\/'); %#ok
        j = 2;
        while (length(rem) > 0)
            [str{j}, rem] = strtok(rem, '\/'); %#ok
            j = j + 1;
        end;
        if (length(str) > 1)
            if (strcmp(token(1), '/'))
                line = ['/', str{1}, '/../', str{end - 1}, '/', str{end}];
            else
                line = [str{1}, '/../', str{end - 1}, '/', str{end}];
            end;
        else
            line = str(1);
        end;
        line = sprintf('  In %s at line %s', char(line), tokens{2});
    end;
    if (isempty(shortened_text))
        shortened_text = char(line);
    else
        shortened_text = sprintf('%s\n%s', char(shortened_text), char(line));
    end;
end;
