function self = run(self, test, dock, shorten)
%gui_test_runner/run is an alternative method to execute the graphical 
%user interface of mlUnit. 
%
%  Example
%  =======
%  With the second parameter it is possible to specify a test method, test
%  case or test suite, that is executed immediately after the start of the
%  graphical user interface. 
%         Example: run(gui_test_runner, 'mlunit_all_tests');
%
%  See also GUI_TEST_RESULT, GUI_TEST_RUNNER.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: run.m 166 2007-01-04 21:19:31Z thomi $

if (nargin == 1)
    test = '';
    dock = 0;
    shorten = 0;
elseif (nargin == 2)
    dock = 0;
    shorten = 0;
elseif (nargin == 3)
    shorten = 0;
end;
self.test_case = test;
self.dock = dock;
self.shorten = shorten;
self.callback = 0;
gui(self);
self.handle = get_handle(get_object(self));