function self = close(self)
%gui_test_runner/close closes the graphical user interface of mlUnit.
%
%  Example
%  =======
%         run(gui_test_runner, 'mlunit_all_tests');
%         close(gui_test_runner);
%
%  See also GUI_TEST_RUNNER.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: close.m 158 2007-01-03 20:25:39Z thomi $

object = get_object(self);
if (~isempty(object) && (strcmp('gui_test_runner', class(object))))
    handle = get_handle(object);
    if (~isempty(handle))
        close(handle);
    end;
end;