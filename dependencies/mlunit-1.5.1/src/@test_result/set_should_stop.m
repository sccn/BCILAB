function self = set_should_stop(self)
%test_result/set_should_stop indicates that the execution of tests
%should stop.
%
%  Example
%  =======
%         result = test_result;
%         % Do something, e.g. iterate through a number of tests, ...
%         result = set_should_stop(result);

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: set_should_stop.m 30 2006-06-11 15:53:00Z thomi $

self.should_stop = 1;