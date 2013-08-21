function groups = find_blocks(in,th_max)
if nargin<2
   th_max = 0.6;
end

groups = [];    
[xl yl] = size(in);
% maxth = 0.05;
% minth = 0.3;
maxi = max(max(in));
% if maxi<0.8
%     return
% end
maxth = maxi*th_max
if (max(max(in)) / (sum(sum(in))/(xl*xl-xl))) < 3
    fprintf('Clustering using MI may not be accurate since MI matrix has low variance')
    fprintf('\n')
end
i = 1;
found = 1;
while found
    maxconf = 0;
    for m_size = xl:-1:2
        
        found = 0;
        for tempi = 1:xl-m_size+1
            temp_matrix = in(tempi:tempi+m_size-1,tempi:tempi+m_size-1);
            [dependent conf] = check_dependency(temp_matrix,maxth);
            
            if dependent
                if conf>maxconf
                    maxconf = conf;
                    res = [tempi m_size];
                end
                found = 1;
            end
        end
        
    end
    if found
        groups = [groups res maxconf];
        in(res(1),:) = 0;
        in(:,res(1)) = 0;
        in(:,res(1)+res(2)-1)=0;
        in(res(1)+res(2)-1,:)=0;
        in(res(1):res(1)+res(2)-1,res(1):res(1)+res(2)-1) = 0;
    end
end




