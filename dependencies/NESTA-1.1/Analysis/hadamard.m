function y = hadamard(x)
% y = hadamard(x)
%   x is a column vector or matrix
%   where the number of rows, m, must
%   be a power of 2
%
%   The transformation is symmetric and orthogonal.
%   To invert, apply again and then divide by m
%
%   This is the Help file for the mex file.
%   Mex file written by Peter Stobbe, stobbe@acm.caltech.edu, Aug 2008
%
%   This is NOT the builtin Matlab code "hadamard"
%
%   Note: in R2008b, Matlab added "fwht" and "ifwht" (the Fast Walsh-
%       Hadamart Transform and the inverse) to its Signal Processing
%       Toolbox.  With the default ordering and scaling, it's not
%       equivalent to this, but you can change this with the following:
%       y = length(x) * fwht( x, [], 'hadamard' );
%       Then y should be the same as hadamard(x) up to roundoff.
%       However, it appears that this code is faster than fwht.


error('You shouldn''t see this code -- mex file is not installed!');
