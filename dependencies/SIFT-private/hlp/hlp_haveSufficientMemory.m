function [res, bytesReq] = hlp_haveSufficientMemory(numelts,datatype)
% Check if we have sufficient memory to allocate an array.
% 
% Inputs:
%   numelts: The number of elements in the array. 
%       If the array is sparse, then numelts is the number of non-zero elements
%   datatype: The class of the array. 
%       This can be any of the following: logical,int*,uint*,single,double,char
% 
% Output: 
%   res: true if we have sufficient memory, else false
%   bytesReq: the number of bytesRequired to allocate the array
%
% Author: Tim Mullen, SCCN/INC/UCSD, Jan, 2014

switch datatype
    case {'logical','int8','uint8'} , N=1;
    case {'char','int16','uint16'}  , N=2;
    case {'double','int64','uint64'}, N=8;
    case {'single','int32','uint32'}, N=4;
    otherwise
        error('unsupported datatype %s',datatype);
end

bytesReq = numelts*N;
res = (bytesReq <= hlp_getAvailableMemory('Bytes'));