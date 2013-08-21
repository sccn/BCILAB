function success = was_successful(self)
%test_result/was_successful returns whether the test was successful or not.
%
%  Example
%  =======
%  was_successful is called for example from text_test_result/run:
%         was_successful(result)
%
%  See also TEXT_TEST_RESULT/RUN.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: was_successful.m 30 2006-06-11 15:53:00Z thomi $

if (size(self.errors, 1) + size(self.failures, 1) == 0)
    success = 1;
else
    success = 0;
end;

