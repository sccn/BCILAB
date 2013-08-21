function res = cropvector(in,startpoints,endpoints)

for i = 1:length(startpoints)+1
    if i ==1
        res = in(:,1:startpoints(i)-1);
    else
        if i ~=length(startpoints) + 1
            res = [res in(:,endpoints(i-1)+1:startpoints(i)-1)];
        else
            res = [res in(:,endpoints(i-1)+1:end)];
        end
    end
end