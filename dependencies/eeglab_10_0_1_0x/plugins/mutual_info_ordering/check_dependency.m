function [res conf] = check_dependency(in,maxth)

%maxth = 0.05;
minth = 0.5;
c1 = 0.5;
c2 = 10;
c3 = 1;
conf = [];

maxi = max(max(in));
[a b] = size(in);
avg = sum(sum(in))/(a*b-a);
for i = 1:a
    in(i,i) = Inf;
end
mini = min(min(in));
if maxi<maxth
    res = 0;
    return;
else
    if mini<minth*maxi
        res = 0;
        return
    else
        res = 1;
        conf = (c2*(maxi/maxth) + c3*avg);
    end
end

    
