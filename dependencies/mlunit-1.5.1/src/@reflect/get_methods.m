function methds = get_methods(self)
%reflect/get_methods returns the list of methods of the 'reflected' class.
%
%  Example
%  =======
%         r = reflect('test_case');
%         get_methods(r);           % Returns a cell array with all methods
%                                   % of the class test_case
%
%  See also REFLECT.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: get_methods.m 23 2006-05-26 23:32:58Z thomi $

methds = [];
if (length(self.class_name) > 0)
    d = methods(self.class_name);
    for i = 1:size(d, 1)
        method = cellstr(d(i));
        if (~strcmp(self.class_name, method))
            methds = [methds; method];
        end;
    end;
end;