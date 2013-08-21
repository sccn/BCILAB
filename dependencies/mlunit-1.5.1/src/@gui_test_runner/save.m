function save(self)
%gui_test_runner/save saves the parameters of the graphical user interface
%of mlUnit to a .mat-file with the name mlunit.tmp. 
%
%  Example
%  =======
%         save(gui_test_runner);
%
%  See also GUI_TEST_RUNNER, GUI_TEST_RUNNER/LOAD.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: save.m 166 2007-01-04 21:19:31Z thomi $

object = get_object(self);
if (strcmp(class(object), 'gui_test_runner')) 
    dock = object.dock; %#ok
    shorten = object.shorten; %#ok
    test_case_name = get(object.handles.gui_test_case, 'String'); %#ok
    save('mlunit.tmp', 'dock', 'shorten', 'test_case_name');
end;