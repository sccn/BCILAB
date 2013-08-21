function self = set_name(self, name)
%test_suite/set_name sets an optional name for the test suite.
%  The name is used by gui_test_runner to re-run a test_suite, which is
%  created by an .m-file.
%
%  Example:
%         function suite = all_tests
%
%         suite = test_suite;
%         suite = set_name(suite, 'all_tests');
%         suite = add_test(suite, my_test('test_foo'));
%         suite = add_test(suite, my_test('test_bar'));

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: set_name.m 63 2006-09-21 20:05:21Z thomi $

self.name = name;