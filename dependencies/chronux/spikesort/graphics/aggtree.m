function aggtree(spikes)

%  AGGTREE  Generates a graphical summary of an aggregation procedure.
%      AGGTREE(spikes)
% Takes a spike-sorting data structure that has undergone hierarchical
%   aggregation and visualizes the hierarchy.
% The function draws its output over the current axis.  The labels of the
%   initial clusters are shown at the bottom of the axis and bracket lines
%   indicate cluster merging.  The height of a bracket indicates the degree
%   of overlap between the two daughter clusters; short brackets imply high
%   overlap.  The bracket color evolves over the 'winter' colormap, with blue
%   representing the first aggregation step and green representing the last.
%   Finally, the number of spikes and isi statistic for each final cluster
%   are given at the top of the plot.

if (~isfield(spikes, 'hierarchy'))
	error('The input data structure to AGGTREE does not have a hierarchy defined.');
end
tree = spikes.hierarchy.tree;
assignments = spikes.hierarchy.assigns;
times = spikes.spiketimes;

% Get some numbers set up . . .
numorigclusts = length(unique(tree(:,[1,2])));  % note that maxorigclusts can be
maxorigclusts = max(unique(tree(:,[1,2])));     %    higher if there are empty clusters
aggsteps = size(tree,1);
numnodes = maxorigclusts + aggsteps;

% . . . and make some room to store tree info.
node = repmat(struct('parent', 0, 'lchild', 0, 'rchild', 0), [numnodes,1]);
xpos = zeros(numnodes,1);
ypos = zeros(numnodes,1);

% modify tree . . .
atree = zeros(aggsteps,4);
atree(:,[1,2]) = tree(:,[1,2]);   % original cluster names
atree(:,4) = 1.1*max(tree(:,3)) - tree(:,3);  % small when connect_strength is big

% Give every node a unique name; the 'aggregate' function output reuses label
%  names, but we'll want a separate identifier for each.
for step = 1:aggsteps
    oldname = atree(step,1);
    newname = maxorigclusts + step;
    atree(step,3) = newname;
    
    atree(find(atree(step+1:end, 1) == oldname) + step, 1) = newname;
    atree(find(atree(step+1:end, 2) == oldname) + step, 2) = newname;
end

% Build tree structure by assigning left/right child and parent indices
% to each aggregation node, determining node y locations along the way.
for step = 1:aggsteps
    lchild = atree(step,1);
    rchild = atree(step,2);
    parent = atree(step,3);
 
    % y position is the greater of the parents positions + an amount
    % derived from connection strength (computed above as 'atree(:,4)')
    ypos(parent) = max(ypos([lchild,rchild])) + atree(step,4);

    % store parent/child indices
    node(parent).lchild = lchild;
    node(parent).rchild = rchild;
    node(lchild).parent = parent;
    node(rchild).parent = parent;
end

% Final clusters are those with no parents (i.e., tree roots).
treeroots = find(cat(2,node.parent) == 0);

% Assign node x locations to leafs (i.e., original clusters labels)
% using a (LIFO) stack to walk the tree left to right.
nextleaf = 1;
stack = treeroots;
while (~isempty(stack))
    current = stack(1);
    stack = stack(2:end);
    if (node(current).lchild ~= 0)  % if not a leaf, push children (left first)
        stack = [node(current).lchild, node(current).rchild, stack];
    else   % if we're at a leaf, take the next available label
        xpos(current) = nextleaf;
        nextleaf = nextleaf + 1;
    end
end

% For interior (i.e., non-leaf nodes), x pos is just average of children's xpos
for this = 1:numnodes
    if (xpos(this) == 0)  % only interior nodes are unassigned
        xpos(this) = mean([xpos(node(this).lchild), xpos(node(this).rchild)]);
    end
end

% Make labels for the leafs (easy because they're ordered by their
% original assignment labels -- so we just write these down and resort
% based on y position to isolate leaves and x position to take the leaves
% from left to right.)
leaflabels = sortrows([(1:numnodes)' xpos ypos], [3,2]); % sort y then x
leaflabels(ypos ~= 0,:) = [];    % leafs have y pos == 0
leaflabels = num2str(leaflabels(:,1));

% Draw dots at the nodes . . .
cla
plot(xpos,ypos,'.','MarkerFaceColor',[0 0 1]);
hold on;

% Now the tree itself; lines with colors indicating aggregation order
cmap = winter(256);
cind = floor(linspace(1,256,aggsteps));  % make use of the entire colormap
for step = 1:aggsteps  % draw the brackets for each step
    x1 = xpos(atree(step,1));
    y1 = ypos(atree(step,1));
    x2 = xpos(atree(step,2));
    y2 = ypos(atree(step,2));
    y3 = ypos(atree(step,3));
    line([x1 x1 x2 x2], [y1 y3 y3 y2], 'Color', cmap(cind(step),:));
end

% Finally, draw extra lines to highlight the final clusters and
% display info for them.
tmin = size(spikes.waveforms,2)./spikes.Fs;
if (isfield(spikes, 'options') && isfield(spikes.options, 'refractory_period'))
    tref = spikes.options.refractory_period;
else   tref = max(0.002, tmin*1.5);
end
yscale = max(max(ypos)*1.5, 0.1);
for root = 1:length(treeroots)
    x1 = xpos(treeroots(root));
    y1 = ypos(treeroots(root));
    y2 = max(ypos) * 1.1;  % place final cluster node higher than the rest
    
    % figure out the original label that this root matches
    oldname = tree((atree(:,3) == treeroots(root)),1);
    if (isempty(oldname))
        oldname = treeroots(root);
    end

    % get cluster size and timing information for this final cluster
    members = find(assignments == oldname);
    clustsize = length(members);
    membertimes = times(members);
	[a, scores] = isiQuality(membertimes, membertimes, tmin, 0.010, tref, spikes.Fs);
    
    plot(x1,y2,'o','MarkerFaceColor',[0 0 0]);
    line([x1 x1],[y1 
		y2], 'LineWidth', 2, 'Color', [0 0 0]);
	stagger = (yscale/20) * rem(root-1,3);
	text(x1 - 0.61, y2*1.05 + stagger,    ['Clust #' num2str(oldname)], 'Color', [0 0.5 0.5], 'FontSize', 8);
    text(x1 - 0.59, y2*1.03 + stagger,    ['N=' num2str(clustsize)], 'Color', [1 0 0], 'FontSize', 8);
    text(x1 + 0.22, y2*0.90 + stagger/3,  ['isi=' num2str(scores(1),2)], 'Color', [0 0 1], 'FontSize', 8);
end

% Prettify the axes
hold off;
colormap(cmap);
%if (length(leaflabels) <= 32)
    set(gca, 'XTickLabel', leaflabels);
    set(gca, 'XTick', (1:size(leaflabels,1)), 'XTickLabelMode', 'manual');
%else
%    set(gca, 'XTickLabel', '');
%end
set(gca,'Xlim',[0 length(leaflabels)+1]);
set(gca,'Ylim',[0 yscale]);
