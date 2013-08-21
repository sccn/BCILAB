function suite = load_tests_from_mfile(self, name) %#ok
%test_loader/load_tests_from_mfile returns a test_suite with all
%test* methods from a .m-file.
%
%  Example
%  =======
%  load_tests_from_mfile is called from within the .m-file, that contains
%  the test* methods, e.g:
%         function test = test_example
%
%         test = load_tests_from_mfile(test_loader);
%
%             function test_method
%                 assert(0 == sin(0));
%             end
%         end
%
%  See also FUNCTION_TEST_CASE.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: load_tests_from_mfile.m 177 2007-01-06 14:02:49Z thomi $

if (nargin == 1)
    name = '';
end;

stack = dbstack;

str = textread(stack(2).file, '%s', 'whitespace', '', 'delimiter', '\n' );
idx = regexp(str, '^\s*function\s+\w*', 'start');
is_func = not(cellfun('isempty', idx));
is_func(1) = 0;
 
tokens = transpose(regexp(str(is_func),...
    '^\s*function\s+(\w*)\s*\%*.*',...
    'tokens'));

set_up_handle = 0;
tear_down_handle = 0;

for token = tokens
    fun = token{1}{:};
    if (strcmp(fun, 'set_up'))
        set_up_handle = evalin('caller', ['@() @', char(fun)]);
        set_up_handle = set_up_handle();
    end;
    if (strcmp(fun, 'tear_down'))
        tear_down_handle = evalin('caller', ['@() @', char(fun)]);
        tear_down_handle = tear_down_handle();
    end;
end;

suite = test_suite;
suite = set_name(suite, stack(2).name);
for token = tokens
    test = char(token{1}{:});
    pos = findstr('test', test);
    if (~isempty(name))
        if (strcmp(name, test))
            fun_handle = evalin('caller', ['@() @', test, '']);
            fun_handle = fun_handle();
            suite = function_test_case(fun_handle,...
                set_up_handle,...
                tear_down_handle);
            break;
        end;
    elseif (~isempty(pos) && (pos(1) == 1))
        fun_handle = evalin('caller', ['@() @', test, '']);
        fun_handle = fun_handle();
        suite = add_test(suite, function_test_case(fun_handle,...
            set_up_handle,...
            tear_down_handle));
    end;
end;


