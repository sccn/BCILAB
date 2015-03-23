%NTH_ELEMENT wrap of C++ nth_element, an efficient rank selection algorithm
%    OUTARR = NTH_ELEMENT(INARR, RANK)
%        INARR is a 2D array of data columns
%        RANK is an integer representing the selected rank to pivot around
%
%        NTH_ELEMENT works with each column in turn and calls C++ 
%            std::nth_element to iteratively pivot until the RANK element
%            is properly placed
%
%        OUTARR is a copy of INARR with the RANK element in the proper
%            position.  All elements before RANK will be less than RANK and
%            all elements after RANK will be greater, but no further sorting
%            is guaranteed.
%
%    See C++ documentation for std::nth_element for more information.
%
%    NTH_ELEMENT will work with any numeric data type except int64.  The 
%    code has lines to handle int64 but they are commented out as the C++ 
%    datatype for int64 is not standard between different compilers.  (GCC 
%    uses "long long" while VC uses __int64.)
%
%    To compile NTH_ELEMENT, you must have MEX set up with a compiler.
%    Then go to the directory that contains nth_element.cpp and run:
%        > mex nth_element.cpp
%
%    As of v0.86, OpenMP pragmas have been added to allow columns to be
%    processed in parallel.  To compile with multithread support, first
%    confirm that your compiler supports OpenMP and then run mex with the
%    appropriate flags.  For example, if you are using GCC >= 4.2, you
%    should be able to do:
%        > mex nth_element.cpp CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
%

% Changes (for whole package including fast_median)
%   Version 0.87 - Bugfix for handling empty inputs.
%   Version 0.86 - Bugfixes, both a compile time bug and a more insidious 
%                  multithreading runtime bug.
%   Version 0.85 - DO NOT USE; BUGS; use 0.86 or higher.
%                  Implemented requested feature, if you request two outputs it
%                  will return the indices of the partitioned elements.  Only
%                  for nth_element at the moment though, not fast_median.
%                  OpenMP for column parallel execution, tried to clean up
%                  includes while maintaining ease of compilation from mex.
%   Version 0.84 - Added in-place versions of nth_element and fast_median
%   Version 0.81 - Changed to BSD license, changed error IDs, added minor 
%                  documentation to CPP files
%   Version 0.8  - Initial release

% Version 0.87
% Peter H. Li 14-NOV-2013
% As required by MatLab Central FileExchange, licensed under the FreeBSD License
