% vec - vectorize an array
%
% Copyright(c) 2009-2011 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt
function V=vec(M)
sz=size(M);
V=reshape(M, [prod(sz), 1]);
