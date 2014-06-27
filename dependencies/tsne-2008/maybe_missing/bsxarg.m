% bsxarg does singleton expansion of two input arrays (related to bsxfun)
%******************************************************************************
% 
%  MATLAB (R) is a trademark of The Mathworks (R) Corporation
% 
%  Function:    bsxarg
%  Filename:    bsxarg.c
%  Programmer:  James Tursa
%  Version:     1.0
%  Date:        February 09, 2008
%  Copyright:   (c) 2008 by James Tursa, All Rights Reserved
%  Permission:  Permission is granted to freely distribute and use this code
%               as long as the header information is included.
% 
%  Returns the singleton expanded arrays that are virtually used in bsxfun.
%  Limited to a maximum of 12 dimensions.
% 
%  bsxfun performs binary operations on input arrays, where singleton
%  dimensions are virtually expanded to perform the operation. bsxarg will
%  actually do the singleton expansion and return the expanded arrays.
% 
%  Building:
% 
%  >> mex -setup
%    (then follow instructions to select a C / C++ compiler of your choice)
%  >> mex bsxarg.c
% 
%  Syntax
% 
%  [C D] = bsxarg(A,B)
% 
%  Description
% 
%  C = the expanded version of A.
%  D = the expanded version of B.
% 
%  Each dimension of A and B must either be equal to each other, or equal to 1.
%  Whenever a dimension of A or B is singleton (equal to 1), the array is 
%  virtually replicated along the dimension to match the other array. The array
%  may be diminished if the corresponding dimension of the other array is 0.
% 
%  The size of the output arrays C and D are equal to:
% 
%  max(size(A),size(B)).*(size(A)>0 & size(B)>0).
% 
%  Examples
%  
%   >> A = 10*rand(2,1,2)
%   A(:,:,1) =
%     9.3181
%     4.6599
%   A(:,:,2) =
%     4.1865
%     8.4622
% 
%   >> B = 10*rand(1,3)
%   B =
%     5.2515    2.0265    6.7214
% 
%   >> [C D] = bsxarg(A,B)
%   C(:,:,1) =
%    9.3181    9.3181    9.3181
%    4.6599    4.6599    4.6599
%   C(:,:,2) =
%    4.1865    4.1865    4.1865
%    8.4622    8.4622    8.4622
% 
%   D(:,:,1) =
%    5.2515    2.0265    6.7214
%    5.2515    2.0265    6.7214
%   D(:,:,2) =
%    5.2515    2.0265    6.7214
%    5.2515    2.0265    6.7214
% 
%  User's with older versions of MATLAB that do not have the bsxfun intrinsic
%  available to them can use this simple m-file to get that capability:
% 
%   function C = bsxfun(fun,A,B)
%   if( nargin ~= 3 )
%       error('Need 3 arguments for bsxfun')
%   end
%   [AX BX] = bsxarg(A,B);
%   if( ischar(fun) )
%       C = eval([fun '(AX,BX)']);
%   else
%       C = fun(AX,BX);
%   end
%   return
%   end
% 
% *******************************************************************************/

function [C D] = bsxarg(A,B)
disp('Error using bsxarg: You have not yet generated the bsxarg mex routine.');
disp('Do the following:');
disp(' ');
disp('>> mex -setup');
disp('  (Then follow instructins to select the C compiler of your choice. (e.g., lcc)');
disp('>> mex bsxarg.c');
disp(' ');
disp('That''s it! Now you are ready to use bsxarg.');
error(' ');
