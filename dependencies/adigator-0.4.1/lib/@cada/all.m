function y = all(x,varargin)
% CADA overloaded version of function ALL - calls cadaunarylogical  
if nargin == 2
  dim = varargin{1};
elseif x.func.size(1) == 1
  dim = 2;
else
  dim = 1;
end
y = cadaunarylogical(x,'all',dim);
