function self = reflect(class_name)
%reflect constructer.
%  The constructor creates an object of the class reflect.
%
%  Class Info / Example
%  ====================
%  The class reflect helps to find out which methods to a class belong.
%  In fact it is simply a wrapper for the Matlab methods function,
%  providing a method checking whether a method within a class exists
%  or not, and a method returning all methods of a class as a cell array. 
%  Example:
%         r = reflect('test_case');
%         method_exists(r, 'run');  % Return true
%         method_exists(r, 'fail'); % Returns false
%         get_methods(r);           % Returns a cell array with all methods
%                                   % of the class test_case

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: reflect.m 23 2006-05-26 23:32:58Z thomi $

self.class_name = class_name;
self = class(self, 'reflect');
