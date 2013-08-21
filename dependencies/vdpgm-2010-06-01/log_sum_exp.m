function y = log_sum_exp(x,dim);
% function y = log_sum_exp(x,dim);
%
% y = log( sum( exp(x), dim ) )
% y = log( sum( exp(x).*y, dim ) )
%
% x can be -inf but cannot be +inf.

% if has_inf(x)
%   warning(['x contains inf; x=' num2str(x)])
% end
% if has_nan(x)
%   x
%   error('x has NaN(s).')
% end

x_size = size(x);
if dim > 2 && dim == length(x_size)
  y_size = x_size(1:end-1);
else
  y_size = x_size;
  y_size(dim) = 1;
end
  
x_max = reshape(max(x, [], dim), y_size);
x_max(find(x_max==-inf)) = 0;
dims = ones(1, ndims(x));
dims(dim) = size(x, dim);
x = x - repmat(x_max, dims);
y = x_max + log(sum(exp(x), dim));

% Local Variables: ***
% mode: matlab ***
% End: ***
