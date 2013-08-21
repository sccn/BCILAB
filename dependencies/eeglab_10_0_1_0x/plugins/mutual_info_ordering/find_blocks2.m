function groups = find_blocks2(in,th_max,start_idx,growing)
start_idx;

if nargin<4
    starting = 1;
    growing = 1;
end
groups = [];
[xl yl] = size(in);
xl;
if xl == 0 || xl == 1
    return
end

% maxth = 0.05;
% minth = 0.3;
maxi = max(max(in));
if (xl == 2)
    [dep conf] = check_dependency(in,th_max);
    if dep
        groups = [start_idx start_idx+1 conf];
    end
else
    % if (max(max(in)) / (sum(sum(in))/(xl*xl-xl))) < 3
    %     fprintf('Clustering using MI may not be accurate since MI matrix has low variance')
    %     fprintf('\n')
    % end
    %i = 1;
    
    [y i] = max(in);
    [dummy idxj]= max(y);
    idxi = i(idxj);
    idxj;
    [dep conf] = check_dependency(in(idxj:idxi,idxj:idxi),th_max);
    if dep
        confn =conf;
        while dep
            conf = confn;
            idxi = idxi + 1;
            if idxi<=xl
                [dep confn] = check_dependency(in(idxj:idxi,idxj:idxi),th_max);
            else
                break;
            end
            
        end
        idxi = idxi -1;
        dep = 1;
        confn = conf;
        while dep
            conf = confn;
            idxj = idxj - 1;
            if idxj>=1
                [dep confn] = check_dependency(in(idxj:idxi,idxj:idxi),th_max);
                
            else
                break
            end
        end
        idxj = idxj + 1;
        groups = [start_idx+idxj-1 start_idx+idxi-1 conf find_blocks2(in(1:idxj-1,1:idxj-1),th_max,start_idx) ...
            find_blocks2(in(idxi+1:end,idxi+1:end),th_max,start_idx+idxi)];
    else
        if (idxi - idxj == 1)
            groups = [find_blocks2(in(1:idxj-1,1:idxj-1),th_max,start_idx) ...
            find_blocks2(in(idxi+1:end,idxi+1:end),th_max,start_idx+idxi)];
            
        else
            [dep conf] = check_dependency(in(idxj:idxi-1,idxj:idxi-1),th_max);
            if ~dep
                [dep conf] = check_dependency(in(idxj+1:idxi,idxj+1:idxi),th_max);
                if~dep
                    groups = find_blocks2(in(idxj+1:idxi-1,idxj+1:idxi-1),th_max,start_idx+idxj,0);
                else
                    confn = conf;
                    while dep
                        conf = confn;
                        idxi = idxi + 1;
                        if idxi<=xl
                            [dep confn] = check_dependency(in(idxj:idxi,idxj:idxi),th_max);
                        else
                            break;
                        end
                    end
                    idxi = idxi - 1;
                    groups = [start_idx+idxj start_idx+idxi-1 conf find_blocks2(in(1:idxj,1:idxj),th_max,start_idx) ...
                        find_blocks2(in(idxi+1:end,idxi+1:end),th_max,start_idx+idxi)];
                end
            else
                confn = conf;
                while dep
                    conf = confn;
                    idxj = idxj - 1;
                    if idxj>=1
                        [dep confn] = check_dependency(in(idxj:idxi,idxj:idxi),th_max);
                        
                    else
                        break
                    end
                end
                idxj = idxj + 1;
                groups = [start_idx+idxj-1 start_idx+idxi-2 conf find_blocks2(in(1:idxj-1,1:idxj-1),th_max,start_idx) ...
                    find_blocks2(in(idxi:end,idxi:end),th_max,start_idx+idxi-1)];
            end
            
        end
    end
end










% found = 1;
%
% while found
%     maxconf = 0;
%     for m_size = xl:-1:2
%
%         found = 0;
%         for tempi = 1:xl-m_size+1
%             temp_matrix = in(tempi:tempi+m_size-1,tempi:tempi+m_size-1);
%             [dependent conf] = check_dependency(temp_matrix,maxth);
%
%             if dependent
%                 if conf>maxconf
%                     maxconf = conf;
%                     res = [tempi m_size];
%                 end
%                 found = 1;
%             end
%         end
%
%     end
%     if found
%         groups = [groups res maxconf];
%         in(res(1),:) = 0;
%         in(:,res(1)) = 0;
%         in(:,res(1)+res(2)-1)=0;
%         in(res(1)+res(2)-1,:)=0;
%         in(res(1):res(1)+res(2)-1,res(1):res(1)+res(2)-1) = 0;
%     end
%     endonf] = check_dependency(in(idxj+1:idxi,idxj+1:idxi),th_max);
%     if~dep
%         groups = findblocks2(in(idxj+1:idxi-1,idxj+1:idxi-1),th_max,0,start_idx+idxj)
%     else
%
%
%
%     end








%     found = 1;
%
%     while found
%         maxconf = 0;
%         for m_size = xl:-1:2
%
%             found = 0;
%             for tempi = 1:xl-m_size+1
%                 temp_matrix = in(tempi:tempi+m_size-1,tempi:tempi+m_size-1);
%                 [dependent conf] = check_dependency(temp_matrix,maxth);
%
%                 if dependent
%                     if conf>maxconf
%                         maxconf = conf;
%                         res = [tempi m_size];
%                     end
%                     found = 1;
%                 end
%             end
%
%         end
%         if found
%             groups = [groups res maxconf];
%             in(res(1),:) = 0;
%             in(:,res(1)) = 0;
%             in(:,res(1)+res(2)-1)=0;
%             in(res(1)+res(2)-1,:)=0;
%             in(res(1):res(1)+res(2)-1,res(1):res(1)+res(2)-1) = 0;
%         end
%     end
