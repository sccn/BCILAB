function varargout = chopdeal(data,sizes)
% Chop a vector into multiple vectors of varying sizes.
% Outputs... = chopdeal(Data,Sizes)
%
% This function operates simiarly to mat2cell, except that it accepts only vectors and returns
% the results as multiple outputs rather than in single cell array. This makes the function 
% equivalent to celldeal(mat2cell(Data,1,Sizes)).
%
% In:
%   Data : a vector to chop (single, double, or cell array).
%
%   Sizes : a double vector of sizes for the chop outputs; 
%           the sizes must sum to the length of the data.
%
% Out:
%   Outputs... : one output for each size, each of size [1xSizes(k)].
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-01

% the below implementation is the MATLAB fallback for the mex function of same name.
if nargout ~= length(sizes)
    error('The number of output arguments must match the number of elements in sizes.'); end
varargout = cell(1,length(sizes));
p=0;
for s=1:length(sizes)
    varargout{s} = data(p+(1:sizes(s)));
    p = p+sizes(s);
end
