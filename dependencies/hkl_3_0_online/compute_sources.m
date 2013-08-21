function sources = compute_sources(hull,qs);
[nhull,p] = size(hull);
active = find(sum(hull-1,1)>0);
inactive = 1:p;
inactive(active) = [];

[a,b] = sort(sum(hull(:,active)-1,1),'descend');
active = active(b);

 
sources_loc = find_sources_complement_grid_fast_int_all_active_c(int32(hull(:,active)'),int32(qs(active)));

sources = ones(size(sources_loc,1),p);
sources(:,active) = sources_loc;
additional_sources = ones(length(inactive),p);
additional_sources(:,inactive) = eye(length(inactive))+1;
sources = [ sources; additional_sources];

