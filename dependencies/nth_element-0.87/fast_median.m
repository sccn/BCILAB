%FAST_MEDIAN using C++ nth_element for better performance
%    MEDIANS = FAST_MEDIAN(ARR)
%        ARR is a 2D array of N columns of data
%        FAST_MEDIAN works with each column in turn
%        Returns MEDIANS, a row vector with the median for each column
%
%    As of R2010a, the MathWorks MEDIAN function relies on SORT.  Sorting
%    the entire dataset to get the median is wasteful.  C++ specifies a
%    more efficient selection algorithm called nth_element.  Typically this 
%    is implemented as a variation on "quickselect", AKA Hoare's Selection
%    Algorithm. 
%
%    (http://en.wikipedia.org/w/index.php?title=Selection_algorithm&oldid=397210793)
%
%    Theoretically, this should have average run time of O(n) instead of
%    the theoretical best case run time of O(n log n) for a SORT based 
%    algorithm.
%
%    FAST_MEDIAN simply wraps the C++ std::nth_element call, also using
%    some simple arithmetic to determine the proper rank for the median.
%    In order to emulate the MatLab MEDIAN function, FAST_MEDIAN treats
%    the median of an even number of elements as the arithmetic mean of
%    the "upper" and "lower" median.
%
%    Also emulating MEDIAN, FAST_MEDIAN will return medians in the same
%    datatype as the data passed in.  FAST_MEDIAN will work with any
%    numeric type except int64.  The code has lines to handle int64 but 
%    they are commented out as the C++ datatype for int64 is not standard
%    between different compilers.  (GCC uses "long long" while VC uses 
%    __int64.)  Note that for int types, FAST_MEDIAN does not round the
%    result; it is simply floored as is typical for integer arithmetic. 
%    This may result in FAST_MEDIAN differing in output from MEDIAN for
%    some integer inputs.
%
%    In contrast to MEDIAN, which attempts to return NaN for empty inputs,
%    FAST_MEDIAN returns empty output of the proper datatype for empty
%    inputs.
%
%    To compile FAST_MEDIAN, you must have MEX set up with a compiler. 
%    Then go to the directory containing fast_median.cpp and run:
%        > mex fast_median.cpp
%

% Version 0.87
% Peter H. Li 14-NOV-2013
% As required by MatLab Central FileExchange, licensed under the FreeBSD License
