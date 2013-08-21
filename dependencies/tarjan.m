function [components] = tarjan(graphMatrix)
%TARJAN finds Strongly Connected Component
%
% Synopsis
%   [COMPONENTS] = TARJAN(G)
%
% Description
%    COMPONENTS = TARJAN(G) searches for strongly connected components
%    using Tarjan's algorithm (it's actually depth first search).
%    G is an input directed graph. The function returns a vector COMPONENTS.
%    The value COMPONENTS(X) is number of component where the node X
%    belongs to.
%
% See also GRAPH/GRAPH, GRAPH/SPANNINGTREE.

%   Author(s): Panacek Martin (MartinPanacek@seznam.cz), P. Sucha
%   Copyright (c) 2007 CTU FEE
%   $Revision: 1896 $  $Date: 2007-10-12 08:13:54 +0200 (pá, 12 X 2007) $

% This file is part of Scheduling Toolbox.
% 
% Scheduling Toolbox is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License as
% published by the Free Software Foundation; either version 2 of the
% License, or (at your option) any later version.
% 
% Scheduling Toolbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with Scheduling Toolbox; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
% USA


global vertex_counter component_counter components dfs back stack stack_head;
% components(x) = component number

n = size(graphMatrix,1); % count of vertices
vertex_counter = 0;
component_counter = 0;
dfs = zeros(1,n);
back = zeros(1,n);
components = ones(1,n) * -1;
stack_head = 0;
stack = zeros(1,n);

for vertex = 1:n
    if components(vertex) == -1
        dfsearch(vertex, graphMatrix);
    end
end

return



% depth first search  
function [] = dfsearch(v , graphMatrix)
global vertex_counter component_counter components dfs back stack stack_head;
% stack + stack_head - stack of passed nodes

vertex_counter = vertex_counter + 1;
dfs(v) = vertex_counter;
back(v) = dfs(v);
components(v) = 0;
stack_head  = stack_head + 1;
stack_index = stack_head;
stack(stack_head) = v;


next = find(graphMatrix(v,:));
n = size(next,2);

for vertex = 1:n
    ver = next(vertex);
    if components(ver) == 0
        back(v) = min(back(ver), back(v));
    elseif components(ver) == -1
        dfsearch(ver, graphMatrix);
        back(v) = min(back(ver), back(v));
    end
end

if dfs(v) == back(v)
    component_counter = component_counter + 1;
    members = stack(stack_index : stack_head);
    components(members) = component_counter;
    stack(stack_index:stack_head) = 0;
    stack_head = stack_index -1;
end

return

%end of file
