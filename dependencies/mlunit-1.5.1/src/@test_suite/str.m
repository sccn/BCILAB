function s = str(self)
%test_suite/str returns a string with the name of the test suite.
%  The name has to be set with test_suite/set_name first.
%
%  Example: 
%         str(mlunit_all_tests)

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: str.m 63 2006-09-21 20:05:21Z thomi $

s = self.name;

