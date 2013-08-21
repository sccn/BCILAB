function self = update(self)
%gui_test_result/update update the different gui object: the progress bar,
%the text area with the number of test ran etc., and the listbox with all
%errors and failures.
%
%  Example
%  =======
%  update is called after the adding of each error, failure or success in
%  add_error, etc:
%         self = update(self);
%
%  See also GUI_TEST_RESULT/ADD_ERROR, GUI_TEST_RESULT/ADD_FAILURE,
%           GUI_TEST_RESULT/ADD_SUCCESS.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: update.m 28 2006-06-11 15:11:59Z thomi $

progress_bar(self);
texts(self);
print_errors(self);
drawnow;
