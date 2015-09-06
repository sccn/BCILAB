
function NodeValue = hlp_computeGraphMeasure(varargin)
% NodeValue = hlp_computeGraphMeasure(causality,ch1,selectedvars,graphMetric)
%
% Compute a univariate graph-theoretic measure for a given node of a graph
%
% Inputs:
%
%   causality:      [num_vars_to x num_vars_from x <num_times> x <num_freqs>] causal matrix
%                   obtained from est_mvarConnectivity().
%   srcNodes:       indices of source nodes
%   selectedvars:   indices of target nodes
%   graphMetric:   which graph measure to compute (see below)
%
% Outputs:
%   
%   NodeValue:      [srcNodes x <num_times> x <num_freqs>] graph measures for the source nodes
%
%
% See Also: vis_causalBrainMovie3D(), est_mvarConnectivity(),
%           hlp_collapseFrequencies()
%
%
% References: 
% 
%   Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual. Chapter 6.
%   Available at: http://www.sccn.ucsd.edu/wiki/SIFT
%
% Author: Tim Mullen, 2010-2013, SCCN/INC/UCSD. 
% Email:  tim@sccn.ucsd.edu

% This function is part of the Source Information Flow Toolbox (SIFT)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

arg_define([0 Inf],varargin, ...
    arg_norep({'cmatrix','CausalMatrix'}, mandatory,[],'Causality array. Shape is [num_vars_to x num_vars_from x <num_times> x <num_freqs>]. Last two dims are optional','shape','matrix'), ...
    arg({'srcNodes','SourceNodes'},  [],[],'Indices of source nodes. Default is all nodes'), ...
    arg({'targNodes','TargetNodes'}, [],[],'Indices of target nodes. Default is all nodes'), ...
    arg({'graphMetric','GraphMetric'},'none',{'none','outflow','mag_outflow','inflow', ...
                                              'causalflow','outdegree','indegree',     ...
                                              'causaldegree','asymmetryratio'}, ...
                                              'Graph measure to compute'), ...
    arg({'ignoreSelfConn','IgnoreSelfConn'},true,[],'Ignore self connectivity in graph measure'));
           
% handle defaults
if isempty(srcNodes)
    srcNodes = 1:size(cmatrix,2);
end
if isempty(targNodes)
    targNodes = 1:size(cmatrix,1);
end

% remove source node from target nodes
if ignoreSelfConn
    cmatrix=hlp_setdiags(cmatrix,0);
%     targNodes = setdiff_bc(targNodes,srcNodes);
end

% compute graph metrics
switch lower(graphMetric)
    case 'none'
        NodeValue = zeros(1,size(cmatrix,3));
    case 'outflow'
        % Compute outflow from srcNodes in each freq
        NodeValue = squeeze(sum(cmatrix(targNodes,srcNodes,:,:),1));
    case 'mag_outflow'
        % Compute outflow from srcNodes in each freq, ignoring sign
        NodeValue = squeeze(sum(abs(cmatrix(targNodes,srcNodes,:,:)),1));
    case 'inflow'
        % Compute inflow to srcNodes in each freq
        NodeValue = squeeze(sum(cmatrix(srcNodes,targNodes,:,:),2));
    case 'causalflow'
        outflow =   squeeze(sum(cmatrix(targNodes,srcNodes,:,:),1));
        inflow  =   squeeze(sum(cmatrix(srcNodes,targNodes,:,:),2));
        NodeValue = outflow - inflow;
    case 'outdegree'
        % Compute number of outgoing edges from srcNodes in each freq
        NodeValue = squeeze(sum(logical(cmatrix(targNodes,srcNodes,:,:)),1));
    case 'indegree'
        % number of incoming edges to srcNodes in each freq
        NodeValue = squeeze(sum(logical(cmatrix(srcNodes,targNodes,:,:)),2));
    case 'causaldegree'
        % outdegree - indegree 
        outflow =   squeeze(sum(logical(cmatrix(targNodes,srcNodes,:,:)),1));
        inflow  =   squeeze(sum(logical(cmatrix(srcNodes,targNodes,:,:)),2));
        NodeValue = outflow - inflow;
    case 'asymmetryratio'
        % 1 if all edges are outgoing, -1 if all edges are incoming.
        % 0 if balanced
        outflow =   squeeze(sum(cmatrix(targNodes,srcNodes,:,:),1));
        inflow  =   squeeze(sum(cmatrix(srcNodes,targNodes,:,:),2));
        NodeValue = (outflow - inflow)./(outflow+inflow);
    otherwise
        % user wants to map a different Conn measure to this
        % (e.g., ERSP)
end

% enforce row vector output
if iscolumn(NodeValue)
    NodeValue = NodeValue';
end