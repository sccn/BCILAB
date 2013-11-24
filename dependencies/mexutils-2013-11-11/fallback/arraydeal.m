function varargout = arraydeal(x)
% MATLAB fallback for arraydeal function
varargout = cell(numel(x),1);
for k=1:numel(x)
    varargout{k} = x(k); end 
