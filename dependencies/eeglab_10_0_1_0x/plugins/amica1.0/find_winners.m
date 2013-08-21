function winners = find_winners(x);
winners = zeros(1,size(x,2),size(x,3));
for i = 1:size(x,3)
    for j = 1:size(x,2)
        [dmmy in] = max(x(:,j,i));
        winners(1,j,i) = in;
    end
end
