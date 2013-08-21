function s = str(self)
%gui_test_runner/str return a string with the class name.
%
%  Example
%  =======
%  For the class gui_test_runner, str will always return:
%           gui_test_runner
%
%  See also TEST_CASE/STR.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: str.m 155 2007-01-03 20:01:35Z thomi $

s = class(self);