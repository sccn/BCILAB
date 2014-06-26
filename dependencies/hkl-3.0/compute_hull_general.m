function hull = compute_hull_withinhull(hull,active,type)

switch type
    case 'grid'
        toadd = [];
        for i=1:size(hull,1)
            toadd = [ toadd; ind2subv(hull(i,:),1:prod(hull( (i),:)))];
        end
        hull = unique([toadd; hull],'rows');

end

