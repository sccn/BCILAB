function suite = load_tests_from_test_case(self, test_case_class)
%test_loader/load_tests_from_test_case returns a test_suite with all
%test* methods from a test_case.
%  It returns an empty matrix, if the test is not found.
%
%  Example
%  =======
%         loader = test_loader;
%         suite = test_suite(load_tests_from_test_case(loader, 'my_test'));


%  This Software and all associated files are released unter the
%  GNU General Public License (GPL), see LICENSE for details.
%
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: load_tests_from_test_case.m 267 2007-03-10 12:38:34Z thomi $

suite = [];
names = get_test_case_names(self, test_case_class);
if (length(names) > 0)
    suite = test_suite(map(self, ...
        test_case_class, ...
        names));
end;
