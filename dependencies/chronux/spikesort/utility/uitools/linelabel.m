function linelabel(vectors)
%LINELABEL         Labels plotted lines.
%   LINELABEL(LIBRARY), allows the user to select points and labels each
%   with the row number of its closest (Euclidean) match among the rows of
%   the matrix LIBRARY.  When the mouse is clicked in the axis, the
%   nearest point is selected and labeled.  The function stops when either
%   the escape or the enter key is pressed.  For this usage, LIBRARY must
%   by an (P x N) matrix, where P is the number of line objects in the
%   current plot and N is the number of points per line.
%
%   LINELABEL(QUERIES), where QUERIES is (M x N) for M less than the
%   number of lines P in the current plot, skips the interactive step.  It
%   instead searches through the lines on the current plot and, for each
%   row of QUERIES, finds the closest (i.e., Euclidean) match.  This match
%   is then labeled on the plot with the index of the corresponding row
%   from QUERIES.
%
%   LINELABEL('reset') deletes all text objects from the current axes and
%   sets the 'LineWidth' of all lines to 1.  USE WITH CARE, since these
%   effects are not restricted to those changes made by LINELABEL.

%%%%%%%%%% SPECIAL CASE
if(ischar(vectors) && strcmp(vectors, 'reset'))
    delete(findobj(gca, 'Type', 'Text'));
    set(findobj(gca, 'Type', 'Line'), 'LineWidth', 1);
    return;
end

%%%%%%%%%% ARGUMENT CHECKING
lines = findobj(gca, 'Type', 'Line'); 
if (isempty(lines)),  error('The plot does not contain any line objects.');  end

ydatalines = get(lines, 'YData');
xdatalines = get(lines, 'XData');
datacolors = get(lines, 'Color');

L = unique(cellfun('length', ydatalines));   % set of lengths of line objects
if (length(L) > 1)  % all lines not same length?
	error('LINELABEL requires all line objects in the current plot to have the same length.');
end
P = size(ydatalines, 1);   [M,N] = size(vectors);
xdatalines = cat(1, xdatalines{:});
ydatalines = cat(1, ydatalines{:});
datacolors = cat(1, datacolors{:});

if (N ~= L)
	error('The input matrix must have the same number of columns as the lines in the current plot.');
elseif (M > P)
	error('The input matrix can not have more rows than the number of lines in the current plot.');
else
	X = unique(xdatalines, 'rows');
	if (size(X, 1) > 1)
		error('LINELABEL requires all line objects in the current plot to share the same XData');
	end
end
if ((xdatalines(1) ~= 1) || (~all(all(diff(xdatalines, 1, 2) == 1))))
    warning(['This function is currently designed to work with XData that ' ... 
            'starts at 1 and is evenly spaced.  Behavior with current plot may be unexpected.']);
end


xlim = get(gca, 'XLim');  ylim = get(gca, 'YLim');
if (M == P)  %%%%%%%%%% INTERACTIVE CASE: label requested points with index of matches from 'vectors'
	while (true)
		[x,y,key] = ginput(1);
		if (isempty(x) || isequal(key, 27)),  break;   end
		nearestX = round(x);
		[howgood, index] = min(abs(ydatalines(:, nearestX) - y));
		if ((abs(((nearestX - x)./(xlim(2)-xlim(1)))) > 0.005) || ...  % too far from a valid x index
		    (abs(((  howgood   )./(ylim(2)-ylim(1)))) > 0.005))        % or too far from a valid y index
			continue;
		end
		[dist,ind] = min(pairdist(vectors, ydatalines(index,:),'nosqrt'), [], 1);
		label = num2str(ind);
        text(x, y, label, 'FontWeight', 'bold', 'FontSize', 14, 'Color', getcolor(datacolors, index));
	end
else         %%%%%%%%%% NONINTERACTIVE CASE: label closest matches to each row of 'vectors' with its index
	if (M > 100)
		areyousure = input(['Warning.  You are trying to label > 100 lines.\n' ...
			                'Enter y to continue or any other key to quit: '], 's');
		if (lower(areyousure(1)) ~= y)
			return;
		end
	end
	for test = 1:M
		match = sum((ydatalines - repmat(vectors(test,:), [size(ydatalines,1), 1])).^2, 2);
		[mn,index] = min(match);
		
		% Find the column that has the largest distance to the closest line (this'll
        % help if the line is an outlier for at least one coordinate).
        ydatacopy = ydatalines;
        ydatacopy(index,:) = Inf;        % don't consider self distance
		dist_to_lines = abs(ydatacopy - repmat(ydatalines(index,:),[P,1]));
		[junk,select] = max(min(dist_to_lines, [], 1), [], 2);    % largest dist to closest line
        text(X(select)*1.01, ydatalines(index,select)*0.99, num2str(test), ...  % this needs to be visible
               'FontWeight', 'bold', 'FontSize', 20, 'Color', getcolor(datacolors, index), ...
               'VerticalAlignment', 'baseline', 'HorizontalAlignment', 'right');
        set(lines(index), 'LineWidth', 3);
        uistack(lines(index), 'top');
	end
end


% Choose a color for the text label 
function color = getcolor(datacolors, index)
if (size(unique(datacolors, 'rows'), 1) > 1) % if the lines aren't all the same color,
    color = brighten(datacolors(index,:), 0.5);
else
    color = [1 1 1] - get(gca, 'Color');
end