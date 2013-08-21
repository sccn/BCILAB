function names = get_test_case_names(self, test_case_class) %#ok
%test_loader/get_test_case_names returns a list of string with all
%test* methods from the test_case_class.
%
%  Example
%  =======
%  get_test_case_names is usually called from
%  test_loader/load_tests_from_mfile:
%         names = get_test_case_names(self, test_case_class);
%
%  See also TEST_LOADER/LOAD_TESTS_FROM_MFILE.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: get_test_case_names.m 148 2007-01-02 20:08:23Z thomi $

t = reflect(test_case_class);
names = get_methods(t);
for i = size(names, 1):-1:1
    if (~strncmp(names(i), 'test', 4))
        names(i) = [];
    end;
end;
names = sortrows(names);
%try/catch mechanism removed as it hides errors in the constructor of the
%test_case_class.
%try
if (length(names) > 0)
    t = eval([test_case_class, '(''', char(names(1)), ''')']);
end;
%catch
%    names = [];
%    return;
%end;

if (~isa(t, 'test_case'))
    names = [];
end;