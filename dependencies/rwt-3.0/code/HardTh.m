function  x = HardTh(y,thld)
%    x = HardTh(y,thld); 
%
%    HARDTH hard thresholds the input signal y with the threshold value
%    thld.
%
%    Input:  
%       y    : 1D or 2D signal to be thresholded
%       thld : threshold value
%
%    Output: 
%       x : Hard thresholded output (x = (abs(y)>thld).*y)
%
%  HERE'S AN EASY WAY TO RUN THE EXAMPLES:
%  Cut-and-paste the example you want to run to a new file 
%  called ex.m, for example. Delete out the % at the beginning 
%  of each line in ex.m (Can use search-and-replace in your editor
%  to replace it with a space). Type 'ex' in matlab and hit return.
%
%
%    Example:
%       y = makesig('WernerSorrows',8);
%       thld = 1;
%       x = HardTh(y,thld)
%       x = 1.5545 5.3175 0 1.6956  -1.2678 0 1.7332 0
%
%    See also: SoftTh
%
%Author: Haitao Guo  <harry@jazz.rice.edu>

x = (abs(y) > thld).*y; 
