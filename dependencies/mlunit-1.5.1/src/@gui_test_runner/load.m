function load(self) %#ok
%gui_test_runner/loads starts the graphical user interface of mlUnit with
%saved parameters from the .mat-file mlunit.tmp.
%
%  Example
%  =======
%         load(gui_test_runner);
%
%  See also GUI_TEST_RUNNER, GUI_TEST_RUNNER/SAVE.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: load.m 166 2007-01-04 21:19:31Z thomi $

try
    saved = load('mlunit.tmp', '-mat');
    run(gui_test_runner, '', saved.dock, saved.shorten);
    object = get_object(gui_test_runner);
    set(object.handles.gui_test_case, 'String', saved.test_case_name);
    delete('mlunit.tmp');
catch
end;

