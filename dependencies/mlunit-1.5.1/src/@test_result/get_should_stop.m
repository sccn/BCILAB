function should_stop = get_should_stop(self)
%test_result/get_should_stop returns whether the test should stop or not.
%
%  Example
%  =======
%  get_should_stop is called for example from test_suite/run:
%         get_should_stop(result)
%
%  See also TEST_SUITE/RUN.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: get_should_stop.m 30 2006-06-11 15:53:00Z thomi $

should_stop = self.should_stop;