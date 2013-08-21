function exists = method_exists(self, method_name)
%reflect/method_exists returns true, if a method with the name method_name 
%exists in the 'reflected' class.
%
%  Example
%  =======
%         r = reflect('test_case');
%         method_exists(r, 'run');  % Return true
%         method_exists(r, 'fail'); % Returns false
%
%  See also REFLECT.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: method_exists.m 23 2006-05-26 23:32:58Z thomi $

exists = ismember(method_name, get_methods(self));