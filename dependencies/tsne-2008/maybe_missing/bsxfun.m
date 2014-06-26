% bsxfun does binary operation with singleton expansion (uses bsxarg)
%--------------------------------------------------------------------------
%  MATLAB (R) is a trademark of The Mathworks (R) Corporation
% 
%  Function:    bsxfun
%  Filename:    bsxfun.m
%  Programmer:  James Tursa
%  Version:     1.0
%  Date:        February 09, 2008
%  Copyright:   (c) 2008 by James Tursa, All Rights Reserved
%  Permission:  Permission is granted to freely distribute and use this code
%               as long as the header information is included.
%
% A poor-man's replacement for bsxfun for those users with earlier versions
% of MATLAB that do not have this function available to them. Uses the mex
% routine bsxarg to physically generate the singleton expanded arrays, and
% then calls the function. This method is not as fast as the MATLAB
% intrinsic bsxfun, but at least it gives you the same functionality.
%
% Apply element-by-element binary operation to two arrays with singleton
% expansion enabled.
%
% bsxfun is an m-file function that mimics the MATLAB intrinsic bsxfun.
% The m-file bsxfun differs from the MATLAB bsxfun in the following ways:
%  - Limited to 12 dimensions
%  - Allows the 1st input to be a char string as well as a function handle
%  - Physically does the singleton expansion of each input array
%
% This bsxfun m-file is mainly intended for users with older versions of
% MATLAB that do not have the MATLAB intrinsic bsxfun available.
%
% REQUIREMENT:  You must build bsxarg first (see bsxarg.c file)
%
% The following are the relavent excerpts from the MATLAB doc:
%
% Syntax
%
% C = bsxfun(fun,A,B)
%
% Description
%
% C = bsxfun(fun,A,B) applies an element-by-element binary operation to
% arrays A and B, with singleton expansion enabled. fun is a function handle or,
% a character string giving the name of a function, and can either be an M-file
% function or one of the following built-in functions:
%
%   @plus     Plus
%   @minus    Minus
%   @times    Array multiply
%   @rdivide  Right array divide
%   @ldivide  Left array divide
%   @power    Array power
%   @max      Binary maximum
%   @min      Binary minimum
%   @rem      Remainder after division
%   @mod      Modulus after division
%   @atan2    Four quadrant inverse tangent
%   @hypot    Square root of sum of squares
%   @eq       Equal
%   @ne       Not equal
%   @lt       Less than
%   @le       Less than or equal to
%   @gt       Greater than
%   @ge       Greater than or equal to
%   @and      Element-wise logical AND
%   @or       Element-wise logical OR
%   @xor      Logical exclusive OR
%
% If an M-file function is specified, it must be able to accept two input
% arrays of the same size as the result C.
%
% Each dimension of A and B must either be equal to each other, or equal to 1.
% Whenever a dimension of A or B is singleton (equal to 1), the array is 
% replicated along the dimension to match the other array. The array
% may be diminished if the corresponding dimension of the other array is 0.
%
% The size of the output array C is equal to:
%
% max(size(A),size(B)).*(size(A)>0 & size(B)>0).
%
% Examples
% 
% In this example, bsxfun is used to subtract the column means from the matrix A.
%
% A = magic(5);
% A = bsxfun(@minus, A, mean(A))
% A =
%    4    11   -12    -5     2
%   10    -8    -6     1     3
%   -9    -7     0     7     9
%   -3    -1     6     8   -10
%   -2     5    12   -11    -4
%
%--------------------------------------------------------------------------

function C = bsxfun(fun,A,B)
if( nargin ~= 3 )
    error('Need 3 arguments for bsxfun')
end
[AX BX] = bsxarg(A,B);
if( ischar(fun) )
    C = eval([fun '(AX,BX)']);
else
    C = fun(AX,BX);
end
return
end
