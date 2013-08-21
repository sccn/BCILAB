function self = test_sin_cos(self)
%test_null checks checks the correlation of cos and sin.
%
%  Example
%  =======
%  Use a test runner to run the test method:
%         Example: run(text_test_runner, test_sin('test_sin_cos'));

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_sin_cos.m 17 2006-05-26 16:23:58Z thomi $

assert_equals(cos(0), sin(pi/2));