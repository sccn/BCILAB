function activehull = compute_hull_withinhull(hull,active,type)
% compute the hull of active variables within the current hull
[ nhull p ] = size(hull);
activehull = 1:nhull;
switch type
	case {'grid + input_space', 'grid + kernels' }

		[ nhull p ] = size(hull);
		if nhull == length(active)
			activehull = 1:nhull;
		else
			toadd = [];
			for i=1:length(active);

				toadd = [ toadd; find( all(hull - repmat(hull(active(i),:),nhull,1)<=0,2) ) ];
			end

			activehull = unique(toadd);

		end
end

