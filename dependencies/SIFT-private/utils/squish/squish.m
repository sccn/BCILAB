function y = squish(x)
% y = squish(x);    SQUISH "x" to remove all singleton dimensions.
%   Since singleton dimensions can confuse many operations, this function removes
%   absolutely ALL singleton dimensions in "x". SQUEEZE is similar, however
%   it will not operate on 2D arrays, of which row vectors are included. Thus the
%   expected result may not always occur with SQUEEZE.
%   NOTE! This function will convert all row vectors to column vectors!
%   Example:        [1;2;3;4;5] = squish(shiftdim([1:5]',-5))

%   To see the differences between SQUISH and SQUEEZE compare the results
%   of the following for any positive or negative n:
%       size(squeeze(shiftdim([1:3]',n)))
%       size(squish(shiftdim([1:3]',n)))
%   created 08/16/2006 by Mirko Hrovat with Matlab ver.7.2
%   modified 02/25/2010 by MIH as per Jan Simon's suggestion (slight speed improvement)
%   contact: mhrovat@email.com

dims = size(x);
y = reshape(x,[dims(dims~=1),1,1]);    % the extra 1's help in case x is a scalar or vector.
%**************** END ***************