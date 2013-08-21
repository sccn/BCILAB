function g = drawNetwork(varargin)
% Visualize a graph in a Matlab figure window by specifying an
% adjacency matrix and optionally node labels, descriptions, colors and the
% the layout algorithm. The best layout algorithms require that graphViz be
% installed, available free at <http://www.graphviz.org>.
%%
% Type doc graphViz4Matlab for more details.

    g = graphViz4Matlab(varargin{:}); 
  
end