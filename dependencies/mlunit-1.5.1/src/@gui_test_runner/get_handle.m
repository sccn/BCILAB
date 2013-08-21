function handle = get_handle(self)
%gui_test_runner/get_handle returns the graphics handle of the gui window,
%if the singleton object is passed as the input parameter.
%
%  Example
%  =======
%         run(gui_test_runner, 'mlunit_all_tests');
%         handle = get_handle(get_object(gui_test_runner))
%
%  See also GUI_TEST_RUNNER, GUI_TEST_RUNNER/GET_OBJECT.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: get_handle.m 158 2007-01-03 20:25:39Z thomi $

handle = self.handle;