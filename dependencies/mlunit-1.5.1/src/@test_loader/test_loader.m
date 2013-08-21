function self = test_loader
%test_loader constructor.
%  The constructor creates an object of the class test_loader.
%
%  Class Info / Example
%  ====================
%  The class test_loader is able to create a test_suite with all
%  test* methods from a test_case. 
%  Example:
%         loader = test_loader;
%         suite = test_suite(load_tests_from_test_case(loader, 'my_test'));

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_loader.m 35 2006-06-11 16:37:12Z thomi $

self.dummy = 0;
self = class(self, 'test_loader');