function sources_to_keep = update_sources(qs,sources,addhull_sources);

p = size(sources,2);
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
    rest_of_sources = [sources_to_keep; sources_to_add(2:end,:)];
    local_to_add = sources_to_add(1,:) ;
    toadd  = add_one_source_c(int32(rest_of_sources'),int32(local_to_add),int32(qs));
    
    
%     toadd_old = [];
%     for i=1:p
%         temp  = local_to_add;
%         if (temp(i)<qs(i))   % has not reached the end of the dimension
%             temp(i) = temp(i)+1;
%             % check that temp is not a descendant of any sources
%             if ~any(all( rest_of_sources - repmat(temp,ltoadd+ltokeep-1,1)<=0,2))
%                 toadd_old = [ toadd_old;  temp];
%             end
%         end
%     end
%     
%     if (size(toadd,1)~=size(toadd_old,1)), keyboard; end
%         toadd
%         toadd_old
%         pause
        sources_to_keep = [ sources_to_keep; toadd ];
        sources_to_add(1,:) = [];
        ltoadd = size(sources_to_add,1);
        ltokeep = size(sources_to_keep,1);
    end
    
