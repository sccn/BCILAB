function [ nearestvalue, indexofnearest ] = utl_nearest( valuestosearch, targetvalue )
% For a group of values to search, find values nearest to each targetvalue
% [ nearestvalue, indexofnearest ] = utl_nearest( valuestosearch, targetvalue )
%
% Note: This function relies on the MATLAB built-in function min, so in the
% case of two valuestosearch which are equivalently near to a targetvalue,
% the first value found is returned
%
% In:
%   valuestosearch  : a vector or matrix of values to search
%   targetvalue : a vector of matrix of target values
%
% Out:
%   nearestvalue : the value from valuestosearch which is the nearest to
%   the targetvalue
%   indexofnearest : the index of nearestvalue within valuestosearch, if
%   valuestosearch is a matrix, indexofnearest returns the linear index
%
% Examples:
%   Find the value and index of desired frequency bins
%
%       freqs = 0:(250/256):250;
%       neededbins = [3 7; 8 10; 11 12];
%       [nearestbins, indexofnearest] = utl_nearest(freqs, neededbins);
%
%
%                               Dan Roberts, the MITRE Corporation
%                               2012-11-20


   [dummy, indexofnearest] = arrayfun(@(x) min(abs(valuestosearch(:) - x)), targetvalue); %#ok<ASGLU>
   nearestvalue = valuestosearch(indexofnearest);
    
end

