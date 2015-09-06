function [idx value error] = getindex(A,els)
%
% function [idx value] = getindex(A,els)
% Returns the indices of values in vector A that are nearest to those in vector els
% 
% Inputs:   A   - lookup vector
%           els - values to search for in A
% Outputs:  idx - a vector the lengths of els containing the indices of
%                nearest values in A
%           value - a vector the length of els containing the values of
%                   nearest match for contents of els in A
%           error - vector the length of els containing the lookup error 
%                   (absolute value of the difference between contents of 
%                   els and their nearest values in A)
%
% Tim Mullen 2010, SCCN/INC, UCSD

if isempty(els)
    idx = [1 length(A)];
    return;
end

% if any(els<min(A)) || any(els>max(A))
%     fprintf('warning! some elements of els are not within range of A -- selecting nearest indices\n');
% end

L = length(els);
[value idx error] = deal(zeros(1,L));
for i=1:L
    [error(i) idx(i)] = min(abs(A-els(i)));
    value(i) = A(idx(i));
end
