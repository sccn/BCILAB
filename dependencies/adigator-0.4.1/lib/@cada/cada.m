function x = cada(varID,func,deriv)
% cada constructor function for the ADiGator algorithm. Creates the
% overloaded object given the .func and .deriv properties. If not all .func
% properties are defined, it sets them to empty.
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
if ~isfield(func,'name')
  func.name = [];
end
if ~isfield(func,'size')
  func.size = [0 0];
end
if ~isfield(func,'zerolocs')
  func.zerolocs = [];
end
if ~isfield(func,'value')
  func.value = [];
end
x.id = varID;
x.func = func;
x.deriv = deriv;
x = class(x,'cada');