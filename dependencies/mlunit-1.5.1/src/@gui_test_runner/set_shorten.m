function self = set_shorten(self, value)
%gui_test_runner/set_shorten activates the shortening of directory paths of
%an error message. This is helpful for a small width of the gui window.
%
%  Example
%  =======
%  The method is internal to the mlUnit framework and should not be called
%  directly.
%
%  See also GUI_TEST_RUNNER, GUI_TEST_RUNNER/GUI,
%           GUI_TEST_RUNNER/SHORTEN_ERROR_TEXT. 

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: set_shorten.m 164 2007-01-04 14:21:20Z thomi $

self.shorten = value;