%FAST_MEDIAN_IP In-place editing version of FAST_MEDIAN
%    MEDIANS = FAST_MEDIAN_IP(ARR)
%
% USE AT YOUR OWN RISK!
%
% See FAST_MEDIAN for general notes.
%
% This version edits the array in-place, so after you have run it your 
% original array will have changed.  This is counter to standard Matlab
% style.  It also requires calling undocumented Mathworks internals, so is
% not guaranteed to work at all on all Matlab versions and should be used
% at your own risk.
%

% Version 0.87
% Peter H. Li 26-JUL-2013
% As required by MatLab Central FileExchange, licensed under the FreeBSD License
