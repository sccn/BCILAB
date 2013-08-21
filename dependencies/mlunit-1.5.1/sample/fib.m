function y = fib(x)
%fib is a test-driven Fibonacci with MATLAB and mlUnit.
%
%  Example
%  =======
%  Running fib from the command line:
%
%         >> fib(0)
%         
%         ans =
%
%              0
%
%         >> fib(1)
%         
%         ans =
%
%              1
%
%         >> fib(16)
%         
%         ans =
%
%            987

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: fib.m 3 2006-05-17 17:24:38Z thomi $

if (x == 0)
    y = 0;
elseif (x <= 2) 
    y = 1;
else
    y = fib(x - 1) + fib(x - 2);
end;