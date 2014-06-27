function  x = SoftTh(y,thld)
%    x = SoftTh(y,thld); 
%
%    SOFTTH soft thresholds the input signal y with the threshold value
%    thld.
%
%    Input:  
%       y    : 1D or 2D signal to be thresholded
%       thld : Threshold value
%
%    Output: 
%       x : Soft thresholded output (x = sign(y)(|y|-thld)_+)
%
%  HERE'S AN EASY WAY TO RUN THE EXAMPLES:
%  Cut-and-paste the example you want to run to a new file 
%  called ex.m, for example. Delete out the % at the beginning 
%  of each line in ex.m (Can use search-and-replace in your editor
%  to replace it with a space). Type 'ex' in matlab and hit return.
%
%
%    Example:
%       y = makesig('Doppler',8);
%       thld = 0.2;
%       x = SoftTh(y,thld)
%       x = 0 0 0 -0.0703 0 0.2001 0.0483 0 
%
%    See also: HardTh
%
%    Reference: 
%       "De-noising via Soft-Thresholding" Tech. Rept. Statistics,
%       Stanford, 1992. D.L. Donoho.
%
%Author: Haitao Guo  <harry@jazz.rice.edu>

x = abs(y);
x = sign(y).*(x >= thld).*(x - thld); 
