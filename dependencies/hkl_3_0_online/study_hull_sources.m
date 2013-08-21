load temphull

% OLD
tic;  sources = find_sources_complement_grid_fast_int_c(int32(hull'),int32(qs)); toc


% NEW (ACTIVE ONES REMOVED BY HAND)

tic;
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
additional_sources(:,inactive) = 2*eye(length(inactive));
sources = [ sources; additional_sources];


toc


% UPDATING
load temphull
tic
% put the sources to add at the end of the list of sources
sources_diff = 1:size(sources,1);
sources_diff(addhull_sources) = [];
sources_to_keep = sources([ sources_diff ],:);
sources_to_add = sources([ addhull_sources ],:);
ltoadd = size(sources_to_add,1);
ltokeep = size(sources_to_keep,1);

while ~isempty(sources_to_add)
    % at each iteration, take one of the sources to add, remove it and add
    % all parents, then check descendants
    toadd = [];
    for i=1:p

        temp  = sources_to_add(1,:) ;
        if (temp(i)<qs(i))   % has not reached the end of the dimension
            temp(i) = temp(i)+1;

            % check that temp is not a descendant of any sources
            if ~any(all([sources_to_keep; sources_to_add(2:end,:)] - repmat(temp,ltoadd+ltokeep-1,1)<=0,2))
                toadd = [ toadd;  temp];
            end
        end
    end
    sources_to_keep = [ sources_to_keep; toadd ];
    sources_to_add(1,:) = [];
    ltoadd = size(sources_to_add,1);
    ltokeep = size(sources_to_keep,1);


end
toc


