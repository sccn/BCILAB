function tests = map(self, test_case_class, test_names) %#ok
%test_loader/map returns a list of objects instantiated from the class
%test_case_class and the methods in test_names.
%
%  Example
%  =======
%  If you have for example a test_case my_test with two methods test_foo1
%  and test_foo2, then
%         map(test_loader, 'my_test', {'test_foo1' 'test_foo2'})
%  returns a list with two objects of my_tests, one instantiated with
%  test_foo1, the other with test_foo2.
%
%  See also TEST_LOADER/LOAD_TESTS_FROM_MFILE.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: map.m 35 2006-06-11 16:37:12Z thomi $

tests = {};
for i = 1:length(test_names)
    test = eval([test_case_class, '(''', char(test_names(i)), ''')']);
    tests{i} = test;
end;