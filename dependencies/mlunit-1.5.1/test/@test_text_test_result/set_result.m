function result = set_result(self, result) %#ok
%test_text_test_result/set_result sets up the contents of result with a
%number of sample results.
%
%  Example
%  =======
%         result = set_result(self, result);
%
%  See also TEST_TEXT_TEST_RESULT/TEST_RESULT.

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: set_result.m 48 2006-06-11 19:38:29Z thomi $

result = start_test(result, mock_test('test_method'));
result = add_success(result, mock_test('test_method'));
result = stop_test(result, mock_test('test_method'));
result = start_test(result, mock_test('test_method'));
result = add_error(result, mock_test('test_method'), 'foo error');
result = stop_test(result, mock_test('test_method'));
result = start_test(result, mock_test('test_method'));
result = add_failure(result, mock_test('test_method'), 'foo failure');
result = stop_test(result, mock_test('test_method'));
