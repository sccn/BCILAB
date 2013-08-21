% archive - pack variables into a struct
%
% Copyright(c) 2009 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt
function S=archive(varargin)

S = [];
for i=1:length(varargin)
  name =varargin{i};
  S = setfield(S, name, evalin('caller', name));
end
