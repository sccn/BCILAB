function self = get_object(self) %#ok
%gui_test_runner/get_object returns the singleton object of the gui window.
%
%  Example
%  =======
%         run(gui_test_runner, 'mlunit_all_tests');
%         runner = get_object(gui_test_runner)
%
%  See also GUI_TEST_RUNNER.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: get_object.m 158 2007-01-03 20:25:39Z thomi $

shh = get(0, 'ShowHiddenHandles');
set(0, 'ShowHiddenHandles', 'on');
handle = findall(0, 'Name', 'mlUnit');
self = get(handle, 'UserData');
set(0, 'ShowHiddenHandles', shh);
