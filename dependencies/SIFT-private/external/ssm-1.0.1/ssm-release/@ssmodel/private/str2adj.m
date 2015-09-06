function A = str2adj(AC)
%% Adjacency matrix entries %%
adj_M   = eye(18);
adj_H   = adj_M(1, :);
adj_Z   = adj_M(2, :);
adj_T   = adj_M(3, :);
adj_R   = adj_M(4, :);
adj_Q   = adj_M(5, :);
adj_c   = adj_M(6, :);
adj_a1  = adj_M(7, :);
adj_P1  = adj_M(8, :);
adj_Hd  = adj_M(9, :);
adj_Zd  = adj_M(10, :);
adj_Td  = adj_M(11, :);
adj_Rd  = adj_M(12, :);
adj_Qd  = adj_M(13, :);
adj_cd  = adj_M(14, :);
adj_Hng = adj_M(15, :);
adj_Qng = adj_M(16, :);
adj_Znl = adj_M(17, :);
adj_Tnl = adj_M(18, :);
A       = zeros(length(AC), 18);
for i = 1 : length(AC)
    n   = length(AC{i});
    for j = 1 : n
        switch AC{i}(j)
            case 'H'
                if j < n && AC{i}(j+1) == 'd', A(i, :) = A(i, :) + adj_Hd;
                elseif j < n-1 && strcmp(AC{i}(j+1:j+2), 'ng'), A(i, :) = A(i, :) + adj_Hng;
                else A(i, :) = A(i, :) + adj_H;
                end
            case 'Z'
                if j < n && AC{i}(j+1) == 'd', A(i, :) = A(i, :) + adj_Zd;
                elseif j < n-1 && strcmp(AC{i}(j+1:j+2), 'nl'), A(i, :) = A(i, :) + adj_Znl;
                else A(i, :) = A(i, :) + adj_Z;
                end
            case 'T'
                if j < n && AC{i}(j+1) == 'd', A(i, :) = A(i, :) + adj_Td;
                elseif j < n-1 && strcmp(AC{i}(j+1:j+2), 'nl'), A(i, :) = A(i, :) + adj_Tnl;
                else A(i, :) = A(i, :) + adj_T;
                end
            case 'R'
                if j < n && AC{i}(j+1) == 'd', A(i, :) = A(i, :) + adj_Rd;
                else A(i, :) = A(i, :) + adj_R;
                end
            case 'Q'
                if j < n && AC{i}(j+1) == 'd', A(i, :) = A(i, :) + adj_Qd;
                elseif j < n-1 && strcmp(AC{i}(j+1:j+2), 'ng'), A(i, :) = A(i, :) + adj_Qng;
                else A(i, :) = A(i, :) + adj_Q;
                end
            case 'c'
                if j < n && AC{i}(j+1) == 'd', A(i, :) = A(i, :) + adj_cd;
                else A(i, :) = A(i, :) + adj_c;
                end
            case 'a'
                A(i, :) = A(i, :) + adj_a1;
            case 'P'
                A(i, :) = A(i, :) + adj_P1;
        end
    end
end
