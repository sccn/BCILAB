function [self, result] = run(self, result)
%test_case/run executes the test case and saves the results in result.
%
%  Example
%  =======
%  There are two ways of calling run:
%
%  1) [test, result] = run(test) uses the default test result.
%
%  2) [test, result] = run(test, result) uses the result given as
%     paramater, which allows to collect the result of a number of tests
%     within one test result.
%
%  See also TEST_CASE.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: run.m 249 2007-01-26 22:59:59Z thomi $

if (nargin == 1)
    result = default_test_result(self);
end;

result = start_test(result, self);
try
    try
        self = set_up(self);
    catch
        result = add_error_with_stack(result, self, lasterror);
        return;
    end;
    
    ok = 0;
    try
        method = self.name;
        self = eval([method, '(self)']);
        ok = 1;
    catch
        err = lasterror;
        errmsg = err.message;
        failure = strfind(errmsg, 'MLUNIT FAILURE');
        if (size(failure) > 0)
            result = add_failure(result, ...
                self, ...
                errmsg(failure(1) + 15:length(errmsg)));
        else
            if (~isfield(err, 'stack'))
                err.stack(1).file = char(which(self.name));
                err.stack(1).line = '1';
                err.stack = vertcat(err.stack, dbstack('-completenames'));
            end;
            
            result = add_error_with_stack(result, self, err);
        end;
    end;
    
    try
        self = tear_down(self);    
    catch
        result = add_error_with_stack(result, self, lasterror);
        ok = 0;
    end;

    if (ok)
        result = add_success(result, self);
    end;
catch
end;
result = stop_test(result, self);

