function self = gui_test_runner(callback)
%gui_test_runner contructor.
%  The constructor creates an object of the class gui_test_runner.
%
%  Class Info / Example
%  ====================
%  The class gui_test_runner runs a test_case or test_suite and displays
%  results in a graphical user interface (using gui_test_result). 
%         Example: run(gui_test_runner);
%
%  See also GUI_TEST_RESULT.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: gui_test_runner.m 166 2007-01-04 21:19:31Z thomi $

if (nargin == 0)
    callback = 0;
end;
self.callback = callback;
self.handle = 0;
self.handles = 0;
self.test_case = '';
self.dock = -1;
self.shorten = 0;
self = class(self, 'gui_test_runner');
